// Reads every file named on the command line, writes <file>.l1 and <file>.l2: raw deflate from
// isal_deflate_stateless with GKL's configuration.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <isa-l/igzip_lib.h>
static unsigned char *slurp(const char *p, size_t *n) {
  FILE *f = fopen(p, "rb"); fseek(f, 0, SEEK_END); *n = ftell(f); fseek(f, 0, SEEK_SET);
  unsigned char *b = malloc(*n + 1); fread(b, 1, *n, f); fclose(f); return b;
}
int main(int argc, char **argv) {
  for (int a = 1; a < argc; a++) {
    size_t n; unsigned char *in = slurp(argv[a], &n);
    for (int level = 1; level <= 2; level++) {
      size_t worst = n + 5 * ((n + 65534) / 65535 + 1) + 64;
      unsigned char *out = malloc(worst);
      size_t lb = level == 2 ? ISAL_DEF_LVL2_DEFAULT : ISAL_DEF_LVL1_DEFAULT;
      unsigned char *lvl = calloc(1, lb);
      struct isal_zstream s;
      isal_deflate_stateless_init(&s);
      s.level = level; s.level_buf = lvl; s.level_buf_size = lb;
      s.next_in = in; s.avail_in = n; s.next_out = out; s.avail_out = worst;
      s.end_of_stream = 1; s.flush = NO_FLUSH;
      int rc = isal_deflate_stateless(&s);
      if (rc) { fprintf(stderr, "%s level %d rc %d\n", argv[a], level, rc); return 1; }
      char path[4096]; snprintf(path, sizeof path, "%s.l%d", argv[a], level);
      FILE *f = fopen(path, "wb"); fwrite(out, 1, worst - s.avail_out, f); fclose(f);
      free(out); free(lvl);
    }
    free(in);
  }
  return 0;
}
