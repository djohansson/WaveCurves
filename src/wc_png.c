/* wc_png.c -- dependency-free PNG writer.
 *
 * Uses stored (uncompressed) deflate blocks, so the files are large but the
 * code is ~80 lines and needs no zlib. It exists so that the headless mode can
 * produce a reference image without dragging in an image library.
 */

#include "wc_png.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t crc_table[256];
static bool     crc_ready = false;

static void crc_init(void)
{
    if (crc_ready) return;
    for (uint32_t n = 0; n < 256; ++n) {
        uint32_t c = n;
        for (int k = 0; k < 8; ++k)
            c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : (c >> 1);
        crc_table[n] = c;
    }
    crc_ready = true;
}

static uint32_t crc_update(uint32_t crc, const uint8_t *buf, size_t len)
{
    for (size_t i = 0; i < len; ++i) crc = crc_table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    return crc;
}

static void put_u32(uint8_t *p, uint32_t v)
{
    p[0] = (uint8_t)(v >> 24); p[1] = (uint8_t)(v >> 16);
    p[2] = (uint8_t)(v >> 8);  p[3] = (uint8_t)v;
}

static void write_chunk(FILE *f, const char *type, const uint8_t *data, uint32_t len)
{
    uint8_t hdr[4];
    put_u32(hdr, len);
    fwrite(hdr, 1, 4, f);
    fwrite(type, 1, 4, f);
    if (len) fwrite(data, 1, len, f);

    uint32_t crc = crc_update(0xFFFFFFFFu, (const uint8_t *)type, 4);
    if (len) crc = crc_update(crc, data, len);
    crc ^= 0xFFFFFFFFu;
    put_u32(hdr, crc);
    fwrite(hdr, 1, 4, f);
}

bool wc_png_write(const char *path, const uint32_t *pixels, int w, int h)
{
    crc_init();

    FILE *f = fopen(path, "wb");
    if (!f) return false;

    static const uint8_t sig[8] = { 137, 'P', 'N', 'G', 13, 10, 26, 10 };
    fwrite(sig, 1, 8, f);

    uint8_t ihdr[13];
    put_u32(ihdr, (uint32_t)w);
    put_u32(ihdr + 4, (uint32_t)h);
    ihdr[8] = 8;      /* bit depth   */
    ihdr[9] = 6;      /* colour type: RGBA */
    ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
    write_chunk(f, "IHDR", ihdr, 13);

    /* Raw scanlines with filter byte 0. The source words are RGBA8888 in
     * little-endian order, which is already R,G,B,A byte-wise. */
    size_t stride = (size_t)w * 4 + 1;
    size_t raw_len = stride * (size_t)h;
    uint8_t *raw = (uint8_t *)malloc(raw_len);
    if (!raw) { fclose(f); return false; }
    for (int y = 0; y < h; ++y) {
        uint8_t *row = raw + stride * (size_t)y;
        row[0] = 0;
        memcpy(row + 1, pixels + (size_t)y * w, (size_t)w * 4);
    }

    /* zlib stream: 2-byte header, stored deflate blocks, adler32 trailer. */
    const size_t BLOCK = 65535;
    size_t nblocks = (raw_len + BLOCK - 1) / BLOCK;
    if (nblocks == 0) nblocks = 1;
    size_t z_len = 2 + nblocks * 5 + raw_len + 4;
    uint8_t *z = (uint8_t *)malloc(z_len);
    if (!z) { free(raw); fclose(f); return false; }

    size_t zi = 0;
    z[zi++] = 0x78; z[zi++] = 0x01;
    size_t off = 0;
    for (size_t bi = 0; bi < nblocks; ++bi) {
        size_t n = raw_len - off;
        if (n > BLOCK) n = BLOCK;
        z[zi++] = (bi == nblocks - 1) ? 1 : 0;
        z[zi++] = (uint8_t)(n & 0xFF);
        z[zi++] = (uint8_t)(n >> 8);
        z[zi++] = (uint8_t)(~n & 0xFF);
        z[zi++] = (uint8_t)((~n >> 8) & 0xFF);
        memcpy(z + zi, raw + off, n);
        zi += n;
        off += n;
    }

    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < raw_len; ++i) {
        a = (a + raw[i]) % 65521u;
        b = (b + a) % 65521u;
    }
    put_u32(z + zi, (b << 16) | a);
    zi += 4;

    write_chunk(f, "IDAT", z, (uint32_t)zi);
    write_chunk(f, "IEND", NULL, 0);

    free(z);
    free(raw);
    fclose(f);
    return true;
}
