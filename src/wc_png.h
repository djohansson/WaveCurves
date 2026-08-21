#ifndef WC_PNG_H
#define WC_PNG_H

#include <stdbool.h>
#include <stdint.h>

/* Write an RGBA8888 image (bytes in memory: R,G,B,A) to `path`. */
bool wc_png_write(const char *path, const uint32_t *pixels, int w, int h);

#endif /* WC_PNG_H */
