/* s4p2taps.c — standalone: .s4p -> channel_taps.txt (+ impulse dump)
 *
 *   gcc -O2 -o s4p2taps s4p2taps.c s4p_channel.c -lm
 *   ./s4p2taps channel.s4p [out_taps.txt]
 */
#include <stdio.h>
#include <stdlib.h>
#include "s4p_channel.h"

#define OSF        16
#define DATA_RATE  60e9
#define FS         (OSF * DATA_RATE)
#define MAX_TAPS   4096

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.s4p> [out_taps.txt]\n", argv[0]);
        return 1;
    }
    const char *out = (argc >= 3) ? argv[2] : "channel_taps.txt";

    static double h[MAX_TAPS];
    S4pOpts opt; s4p_default_opts(&opt);

    int L = channel_from_s4p(argv[1], h, MAX_TAPS, FS, &opt);
    if (L <= 0) return 1;

    FILE *fp = fopen(out, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot write '%s'\n", out); return 1; }
    for (int i = 0; i < L; i++) fprintf(fp, "%.12e\n", h[i]);
    fclose(fp);
    printf("  wrote %d taps -> %s\n", L, out);
    return 0;
}