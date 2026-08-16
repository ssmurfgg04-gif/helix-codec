/**
 * Minimal htslib-compatible BAM/CRAM reader compiled to WASM.
 *
 * Provides the core htslib API surface:
 *   - hts_open() / hts_close()
 *   - sam_hdr_read() / sam_hdr_destroy()
 *   - sam_read1()
 *   - bam1_t structure access
 *
 * Uses zlib for BGZF decompression. The BAM binary format parser is
 * implemented directly (same format as htslib, but without the full
 * htslib dependency chain).
 *
 * Compile:
 *   emcc htslib-wasm-core.c -O3 -I../zlib-src -L../zlib-src -lz \
 *     -s EXPORTED_FUNCTIONS='[...]' -o htslib_wasm.js
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <emscripten.h>

/* Minimal zlib declarations (we link against compiled zlib) */
#include "zlib.h"

/* ---------------------------------------------------------------------------
 * BAM format constants
 * --------------------------------------------------------------------------- */

#define BAM_MAGIC 0x014D4142  /* "BAM\1" little-endian */

/* BAM flag bits */
#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024
#define BAM_FSUPPLEMENTARY 2048

/* CIGAR operations */
static const char CIGAR_OPS[] = "MIDNSHP=X";

/* ---------------------------------------------------------------------------
 * htsFile — opaque file handle
 * --------------------------------------------------------------------------- */

typedef enum { HTS_FMT_SAM, HTS_FMT_BAM, HTS_FMT_CRAM, HTS_FMT_UNKNOWN } htsFormat;

typedef struct {
    uint8_t *data;       /* File data (loaded into memory) */
    int32_t   data_len;  /* Total data length */
    int32_t   pos;       /* Current read position */
    htsFormat format;    /* Detected format */
    uint8_t *uncompressed; /* BGZF-decompressed data */
    int32_t   uncomp_len;   /* Decompressed data length */
    int32_t   uncomp_pos;   /* Position in decompressed data */
} htsFile;

/* ---------------------------------------------------------------------------
 * sam_hdr_t — BAM header
 * --------------------------------------------------------------------------- */

typedef struct {
    char    *text;       /* Header text */
    int32_t  l_text;     /* Length of header text */
    int32_t  n_ref;      /* Number of reference sequences */
    char   **ref_name;   /* Reference sequence names */
    int32_t *ref_len;    /* Reference sequence lengths */
    char    *format_version; /* SAM format version */
} sam_hdr_t;

/* ---------------------------------------------------------------------------
 * bam1_t — BAM alignment record
 * --------------------------------------------------------------------------- */

typedef struct {
    int32_t  core_tid;       /* Reference sequence ID (-1 if unmapped) */
    int32_t  core_pos;       /* 0-based leftmost coordinate */
    uint32_t core_bin;       /* BAI bin number */
    uint8_t  core_qual;      /* Mapping quality */
    uint16_t core_l_qname;   /* Length of read name (including NUL) */
    uint16_t core_n_cigar;   /* Number of CIGAR operations */
    uint16_t core_flag;      /* BAM flag */
    int32_t  core_l_qseq;    /* Length of read sequence */
    int32_t  core_mtid;      /* Mate reference ID */
    int32_t  core_mpos;      /* 0-based mate position */
    int32_t  core_isize;     /* Insert size */

    /* Variable-length data */
    uint32_t l_data;         /* Current length of variable-length data */
    uint32_t m_data;         /* Allocated length */
    uint8_t *data;           /* Read name + CIGAR + seq + qual + aux */
} bam1_t;

/* Accessor macros */
#define bam_get_qname(b)  ((char*)((b)->data))
#define bam_get_cigar(b)  ((uint32_t*)((b)->data + (b)->core_l_qname))
#define bam_get_seq(b)    ((uint8_t*)((b)->data + (b)->core_l_qname + (b)->core_n_cigar * 4))
#define bam_get_qual(b)   ((b)->data + (b)->core_l_qname + (b)->core_n_cigar * 4 + (((b)->core_l_qseq + 1) >> 1))
#define bam_get_aux(b)    ((b)->data + (b)->core_l_qname + (b)->core_n_cigar * 4 + (((b)->core_l_qseq + 1) >> 1) + (b)->core_l_qseq)
#define bam_seqi(b, i)    ((i) % 2 ? (bam_get_seq(b)[(i)/2] & 0xF) : (bam_get_seq(b)[(i)/2] >> 4 & 0xF))

/* ---------------------------------------------------------------------------
 * Helper: read little-endian integers from buffer
 * --------------------------------------------------------------------------- */

static inline int32_t read_i32(const uint8_t *p) {
    return (int32_t)(p[0] | (p[1] << 8) | (p[2] << 16) | (p[3] << 24));
}
static inline uint32_t read_u32(const uint8_t *p) {
    return (uint32_t)(p[0] | (p[1] << 8) | (p[2] << 16) | (p[3] << 24));
}
static inline int16_t read_i16(const uint8_t *p) {
    return (int16_t)(p[0] | (p[1] << 8));
}
static inline uint16_t read_u16(const uint8_t *p) {
    return (uint16_t)(p[0] | (p[1] << 8));
}

/* ---------------------------------------------------------------------------
 * BGZF decompression
 * --------------------------------------------------------------------------- */

/**
 * Decompress BGZF data using zlib.
 * BGZF is a series of gzip blocks, each containing a BGZF extra field.
 * We decompress block by block.
 */
static uint8_t *bgzf_decompress(const uint8_t *data, int32_t data_len, int32_t *out_len) {
    /* Allocate output buffer (start with 4x input size) */
    int32_t out_cap = data_len * 4;
    uint8_t *out = (uint8_t *)malloc(out_cap);
    int32_t out_pos = 0;

    int32_t in_pos = 0;
    while (in_pos < data_len) {
        /* Check gzip magic */
        if (in_pos + 10 > data_len || data[in_pos] != 0x1F || data[in_pos+1] != 0x8B) {
            /* Not gzip — treat as uncompressed */
            if (out_pos + data_len > out_cap) {
                out_cap = out_pos + data_len;
                out = (uint8_t *)realloc(out, out_cap);
            }
            memcpy(out + out_pos, data + in_pos, data_len - in_pos);
            out_pos += data_len - in_pos;
            break;
        }

        /* Read gzip header to find compressed block size */
        /* BGZF block structure:
         *   0-1: magic (0x1F 0x8B)
         *   2:   method (8 = deflate)
         *   3:   flags
         *   4-7: mtime
         *   8:   xfl
         *   9:   os
         *   10-11: xlen (extra field length)
         *   12-xlen: extra field (must contain BGZF: 0x42 0x43 0x02 0x01)
         *   After extra: compressed data
         *   Then: CRC32 (4 bytes) + ISIZE (4 bytes)
         */

        /* Read the BGZF extra field to get block size */
        if (in_pos + 12 > data_len) break;

        uint8_t flags = data[in_pos + 3];
        int32_t pos = in_pos + 10;

        /* Skip extra field if FEXTRA flag is set */
        int32_t bsize = 0;
        if (flags & 0x04) { /* FEXTRA */
            if (pos + 2 > data_len) break;
            uint16_t xlen = read_u16(data + pos);
            /* Look for BGZF extra subfield (SI1=66 'B', SI2=67 'C', SLEN=2) */
            int32_t xp = pos + 2;
            while (xp + 4 <= pos + 2 + xlen) {
                if (data[xp] == 0x42 && data[xp+1] == 0x43 && data[xp+2] == 0x02) {
                    bsize = read_u16(data + xp + 3) + 1; /* bsize = total block size */
                    break;
                }
                /* Skip this subfield */
                xp += 4 + data[xp + 3];
            }
            pos += 2 + xlen;
        }

        if (bsize == 0) {
            /* Not BGZF — skip this block (shouldn't happen for valid BGZF) */
            bsize = 65536; /* default guess */
        }

        /* Decompress this block using zlib */
        if (in_pos + bsize > data_len) bsize = data_len - in_pos;

        /* The compressed data starts after the gzip header.
         * For simplicity, use zlib's inflate with raw deflate (wbits=-15). */
        z_stream strm;
        memset(&strm, 0, sizeof(strm));
        strm.zalloc = Z_NULL;
        strm.zfree = Z_NULL;
        strm.opaque = Z_NULL;

        /* Initialize for raw deflate (gzip header already parsed) */
        int ret = inflateInit2(&strm, -15);
        if (ret != Z_OK) {
            in_pos += bsize;
            continue;
        }

        strm.avail_in = bsize - (pos - in_pos); /* compressed data length */
        strm.next_in = (Bytef *)(data + pos);

        /* Allocate space for decompressed output */
        int32_t block_out_cap = 65536; /* BGZF blocks decompress to at most 65536 bytes */
        uint8_t *block_out = (uint8_t *)malloc(block_out_cap);

        strm.avail_out = block_out_cap;
        strm.next_out = block_out;

        ret = inflate(&strm, Z_FINISH);

        if (ret == Z_STREAM_END || ret == Z_OK) {
            int32_t block_out_len = block_out_cap - strm.avail_out;
            /* Grow output buffer if needed */
            if (out_pos + block_out_len > out_cap) {
                out_cap = (out_pos + block_out_len) * 2;
                out = (uint8_t *)realloc(out, out_cap);
            }
            memcpy(out + out_pos, block_out, block_out_len);
            out_pos += block_out_len;
        }

        free(block_out);
        inflateEnd(&strm);
        in_pos += bsize;
    }

    *out_len = out_pos;
    return out;
}

/* ---------------------------------------------------------------------------
 * hts_open / hts_close — file opening
 * --------------------------------------------------------------------------- */

EMSCRIPTEN_KEEPALIVE
htsFile *hts_open_mem(const uint8_t *data, int32_t data_len) {
    htsFile *fp = (htsFile *)calloc(1, sizeof(htsFile));
    if (!fp) return NULL;

    fp->data = (uint8_t *)malloc(data_len);
    memcpy(fp->data, data, data_len);
    fp->data_len = data_len;
    fp->pos = 0;

    /* Detect format */
    if (data_len >= 2 && data[0] == 0x1F && data[1] == 0x8B) {
        /* Gzip/BGZF — could be BAM */
        fp->format = HTS_FMT_BAM;
        /* Decompress BGZF */
        fp->uncompressed = bgzf_decompress(data, data_len, &fp->uncomp_len);
        fp->uncomp_pos = 0;
    } else if (data_len >= 4 && read_u32(data) == BAM_MAGIC) {
        fp->format = HTS_FMT_BAM;
        fp->uncompressed = fp->data;
        fp->uncomp_len = data_len;
        fp->uncomp_pos = 0;
    } else if (data_len > 0 && data[0] == '@') {
        fp->format = HTS_FMT_SAM;
        fp->uncompressed = fp->data;
        fp->uncomp_len = data_len;
        fp->uncomp_pos = 0;
    } else {
        fp->format = HTS_FMT_UNKNOWN;
    }

    return fp;
}

EMSCRIPTEN_KEEPALIVE
void hts_close(htsFile *fp) {
    if (!fp) return;
    if (fp->uncompressed && fp->uncompressed != fp->data) free(fp->uncompressed);
    if (fp->data) free(fp->data);
    free(fp);
}

/* ---------------------------------------------------------------------------
 * sam_hdr_read — read BAM header
 * --------------------------------------------------------------------------- */

EMSCRIPTEN_KEEPALIVE
sam_hdr_t *sam_hdr_read(htsFile *fp) {
    if (!fp || fp->format != HTS_FMT_BAM) return NULL;

    const uint8_t *d = fp->uncompressed;
    int32_t len = fp->uncomp_len;
    int32_t pos = 0;

    /* Verify BAM magic */
    if (len < 4 || read_u32(d) != BAM_MAGIC) return NULL;
    pos = 4;

    /* Read header text */
    if (pos + 4 > len) return NULL;
    int32_t l_text = read_i32(d + pos); pos += 4;
    if (pos + l_text > len) return NULL;

    sam_hdr_t *hdr = (sam_hdr_t *)calloc(1, sizeof(sam_hdr_t));
    hdr->l_text = l_text;
    hdr->text = (char *)malloc(l_text + 1);
    memcpy(hdr->text, d + pos, l_text);
    hdr->text[l_text] = '\0';
    pos += l_text;

    /* Read reference sequences */
    if (pos + 4 > len) { free(hdr->text); free(hdr); return NULL; }
    hdr->n_ref = read_i32(d + pos); pos += 4;
    hdr->ref_name = (char **)calloc(hdr->n_ref, sizeof(char *));
    hdr->ref_len = (int32_t *)calloc(hdr->n_ref, sizeof(int32_t));

    for (int32_t i = 0; i < hdr->n_ref; i++) {
        if (pos + 4 > len) break;
        int32_t l_name = read_i32(d + pos); pos += 4;
        if (pos + l_name > len) break;
        hdr->ref_name[i] = (char *)malloc(l_name);
        memcpy(hdr->ref_name[i], d + pos, l_name - 1); /* exclude NUL */
        hdr->ref_name[i][l_name - 1] = '\0';
        pos += l_name;
        if (pos + 4 > len) break;
        hdr->ref_len[i] = read_i32(d + pos); pos += 4;
    }

    fp->uncomp_pos = pos; /* Save position for sam_read1 */
    return hdr;
}

EMSCRIPTEN_KEEPALIVE
void sam_hdr_destroy(sam_hdr_t *hdr) {
    if (!hdr) return;
    if (hdr->text) free(hdr->text);
    if (hdr->ref_name) {
        for (int32_t i = 0; i < hdr->n_ref; i++) free(hdr->ref_name[i]);
        free(hdr->ref_name);
    }
    if (hdr->ref_len) free(hdr->ref_len);
    free(hdr);
}

/* Accessors for the header */
EMSCRIPTEN_KEEPALIVE
int32_t hdr_n_ref(sam_hdr_t *hdr) { return hdr ? hdr->n_ref : 0; }

EMSCRIPTEN_KEEPALIVE
const char *hdr_ref_name(sam_hdr_t *hdr, int32_t i) {
    return (hdr && i >= 0 && i < hdr->n_ref) ? hdr->ref_name[i] : NULL;
}

EMSCRIPTEN_KEEPALIVE
int32_t hdr_ref_len(sam_hdr_t *hdr, int32_t i) {
    return (hdr && i >= 0 && i < hdr->n_ref) ? hdr->ref_len[i] : 0;
}

EMSCRIPTEN_KEEPALIVE
const char *hdr_text(sam_hdr_t *hdr) { return hdr ? hdr->text : NULL; }

EMSCRIPTEN_KEEPALIVE
int32_t hdr_l_text(sam_hdr_t *hdr) { return hdr ? hdr->l_text : 0; }

/* ---------------------------------------------------------------------------
 * sam_read1 — read next BAM alignment
 * --------------------------------------------------------------------------- */

EMSCRIPTEN_KEEPALIVE
int32_t sam_read1(htsFile *fp, sam_hdr_t *hdr, bam1_t *b) {
    if (!fp || fp->format != HTS_FMT_BAM || !b) return -1;

    const uint8_t *d = fp->uncompressed;
    int32_t len = fp->uncomp_len;
    int32_t pos = fp->uncomp_pos;

    if (pos + 4 > len) return -1; /* EOF */

    /* Read block_size */
    int32_t block_size = read_i32(d + pos); pos += 4;
    if (pos + block_size > len) return -1; /* truncated */

    int32_t block_end = pos + block_size;

    /* Read core fields */
    b->core_tid = read_i32(d + pos); pos += 4;
    b->core_pos = read_i32(d + pos); pos += 4;
    b->core_l_qname = d[pos]; pos += 1;
    b->core_qual = d[pos]; pos += 1;
    b->core_bin = read_u16(d + pos); pos += 2;
    b->core_n_cigar = read_u16(d + pos); pos += 2;
    b->core_flag = read_u16(d + pos); pos += 2;
    b->core_l_qseq = read_i32(d + pos); pos += 4;
    b->core_mtid = read_i32(d + pos); pos += 4;
    b->core_mpos = read_i32(d + pos); pos += 4;
    b->core_isize = read_i32(d + pos); pos += 4;

    /* Read variable-length data */
    int32_t var_len = block_size - 32; /* core is 32 bytes */
    if (var_len < 0) var_len = 0;

    if (b->m_data < (uint32_t)var_len) {
        b->data = (uint8_t *)realloc(b->data, var_len);
        b->m_data = var_len;
    }
    memcpy(b->data, d + pos, var_len);
    b->l_data = var_len;

    fp->uncomp_pos = block_end;
    return 0; /* success */
}

/* ---------------------------------------------------------------------------
 * bam1_t lifecycle
 * --------------------------------------------------------------------------- */

EMSCRIPTEN_KEEPALIVE
bam1_t *bam_init1(void) {
    bam1_t *b = (bam1_t *)calloc(1, sizeof(bam1_t));
    return b;
}

EMSCRIPTEN_KEEPALIVE
void bam_destroy1(bam1_t *b) {
    if (!b) return;
    if (b->data) free(b->data);
    free(b);
}

/* ---------------------------------------------------------------------------
 * bam1_t accessors for JS wrapper
 * --------------------------------------------------------------------------- */

EMSCRIPTEN_KEEPALIVE
int32_t bam_core_tid(bam1_t *b) { return b ? b->core_tid : -1; }
EMSCRIPTEN_KEEPALIVE
int32_t bam_core_pos(bam1_t *b) { return b ? b->core_pos : -1; }
EMSCRIPTEN_KEEPALIVE
uint8_t bam_core_qual(bam1_t *b) { return b ? b->core_qual : 0; }
EMSCRIPTEN_KEEPALIVE
uint16_t bam_core_flag(bam1_t *b) { return b ? b->core_flag : 0; }
EMSCRIPTEN_KEEPALIVE
int32_t bam_core_l_qseq(bam1_t *b) { return b ? b->core_l_qseq : 0; }
EMSCRIPTEN_KEEPALIVE
int32_t bam_core_mtid(bam1_t *b) { return b ? b->core_mtid : -1; }
EMSCRIPTEN_KEEPALIVE
int32_t bam_core_mpos(bam1_t *b) { return b ? b->core_mpos : -1; }
EMSCRIPTEN_KEEPALIVE
int32_t bam_core_isize(bam1_t *b) { return b ? b->core_isize : 0; }

/* Get read name as C string (pointer into bam1_t data) */
EMSCRIPTEN_KEEPALIVE
const char *bam_qname(bam1_t *b) { return b ? bam_get_qname(b) : NULL; }

/* Get CIGAR as uint32 array (pointer into bam1_t data) */
EMSCRIPTEN_KEEPALIVE
const uint32_t *bam_cigar(bam1_t *b) { return b ? bam_get_cigar(b) : NULL; }
EMSCRIPTEN_KEEPALIVE
uint16_t bam_n_cigar(bam1_t *b) { return b ? b->core_n_cigar : 0; }

/* Decode a sequence to ASCII (caller must free) */
EMSCRIPTEN_KEEPALIVE
char *bam_seq_str(bam1_t *b) {
    if (!b || b->core_l_qseq <= 0) return NULL;
    static const char SEQ_LUT[] = "=ACMGRSVTWYHKDBN";
    char *seq = (char *)malloc(b->core_l_qseq + 1);
    uint8_t *s = bam_get_seq(b);
    for (int32_t i = 0; i < b->core_l_qseq; i++) {
        seq[i] = SEQ_LUT[bam_seqi(b, i)];
    }
    seq[b->core_l_qseq] = '\0';
    return seq;
}

/* Get quality string (Phred+33, caller must free) */
EMSCRIPTEN_KEEPALIVE
char *bam_qual_str(bam1_t *b) {
    if (!b || b->core_l_qseq <= 0) return NULL;
    char *qual = (char *)malloc(b->core_l_qseq + 1);
    uint8_t *q = bam_get_qual(b);
    for (int32_t i = 0; i < b->core_l_qseq; i++) {
        qual[i] = q[i] == 0xFF ? '*' : (char)(q[i] + 33);
    }
    qual[b->core_l_qseq] = '\0';
    return qual;
}
