/*
 * Which ISA-L configuration is GKL's level 1?
 *
 * Decision 0029 established that levels 1 and 2 of `IntelDeflater` reach `isal_deflate_stateless`
 * rather than zlib. That names a function, not a configuration, and igzip's output depends on
 * three more things the call site chooses: the ISA-L level, the size of the token buffer, and
 * whether the stream is closed in one call.
 *
 * They are found by trying rather than by reading, because the disassembly is misleading here. The
 * branch contains `calloc(1, 0x141D0)`, which looks exactly like a level buffer and is not one:
 * feeding ISA-L that size reproduces GKL on three fixtures and misses the fourth by 62 bytes,
 * which is the shape of a near miss that a smaller corpus would have called a match.
 *
 * The answer, and the only combination that matches all four:
 *
 *     level           1                          (2 gives 19141 where GKL gives 19044)
 *     level_buf_size  ISAL_DEF_LVL1_DEFAULT      (0x141D0 gives 63373 where GKL gives 63311)
 *     end_of_stream   1
 *
 * and it explains why Java levels 1 and 2 produce identical bytes: GKL does not pass the level
 * through, so both land on ISA-L level 1.
 *
 * Build and run (linux/amd64, against ISA-L 2.30.0):
 *
 *     gcc -O2 -Iisa-l/include igzip_config.c isa-l/.libs/libisal.a -o probe && ./probe
 *
 * with the four fixtures in /fx and the output compared against the `gkl` rows of
 * `real-x86-64.txt`. See docs/decisions/0031-gkl-levels-one-and-two-are-isal-level-one.md.
 */
#include <stdio.h>
#include <stdlib.h>
#include "igzip_lib.h"
static const char *names[] = {"acgt","runs","random","acgt-2blocks"};
int main(void){
  for(int f=0;f<4;f++){
    char p[256]; snprintf(p,sizeof p,"/fx/%s.bin",names[f]);
    FILE*fp=fopen(p,"rb"); fseek(fp,0,SEEK_END); long n=ftell(fp); fseek(fp,0,SEEK_SET);
    unsigned char*in=malloc(n); if(fread(in,1,n,fp)!=(size_t)n) return 1; fclose(fp);
    struct isal_zstream s; isal_deflate_stateless_init(&s);
    unsigned char*lb=calloc(1,ISAL_DEF_LVL1_DEFAULT);
    s.level=1; s.level_buf=lb; s.level_buf_size=ISAL_DEF_LVL1_DEFAULT;
    s.next_in=in; s.avail_in=n;
    unsigned long cap=2*n+4096; unsigned char*o=malloc(cap);
    s.next_out=o; s.avail_out=cap; s.end_of_stream=1;
    isal_deflate_stateless(&s);
    unsigned long len=cap-s.avail_out;
    snprintf(p,sizeof p,"/out/%s.final.raw",names[f]);
    FILE*w=fopen(p,"wb"); fwrite(o,1,len,w); fclose(w);
    printf("%s %lu\n",names[f],len);
  }
  return 0;
}
