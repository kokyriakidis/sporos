#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
static char *seq; static long n; static int K, W, T;
typedef struct { int id; long cnt; } job_t;
static void *worker(void *arg){
  job_t *j=(job_t*)arg; long chunk=(n+T-1)/T; long st=j->id*chunk, en=st+chunk+K+W-2; if(en>n) en=n;
  mm128_v p={0,0,0}; void *km=km_init();
  // mm_sketch takes int len; process in 64 Mbp pieces with overlap
  long cnt=0; for(long s=st; s<en; s+=(1L<<26)){ long e=s+(1L<<26)+K+W-2; if(e>en) e=en; p.n=0; mm_sketch(km, seq+s, (int)(e-s), W, K, 0, 0, &p); cnt+=p.n; if(e==en) break; }
  j->cnt=cnt; km_destroy(km); return 0;
}
int main(int argc, char **argv){
  T = argc>2? atoi(argv[2]) : 1;
  if(!strcmp(argv[1],"random")){ n=100000000; seq=malloc(n); srand(42); for(long i=0;i<n;i++) seq[i]="ACGT"[rand()&3]; }
  else { FILE *f=fopen(argv[1],"rb"); fseek(f,0,SEEK_END); long sz=ftell(f); fseek(f,0,SEEK_SET); char *b=malloc(sz); fread(b,1,sz,f); fclose(f);
    seq=malloc(sz); n=0; int hdr=0; for(long i=0;i<sz;i++){ char c=b[i]; if(c=='>') hdr=1; if(c=='\n'){hdr=0; continue;} if(!hdr) seq[n++]=c; } free(b); }
  int ks[2]={21,31}, ws[2]={11,19};
  for(int t=0;t<2;t++){ K=ks[t]; W=ws[t]; if(K>28){ printf("k=%d > 28 unsupported by mm_sketch; using k=28 w=%d\n",K,W); K=28; }
    double best=1e9; long cnt=0;
    for(int r=0;r<5;r++){ pthread_t th[64]; job_t jb[64]; double t0=now();
      for(int i=0;i<T;i++){ jb[i].id=i; pthread_create(&th[i],0,worker,&jb[i]); }
      cnt=0; for(int i=0;i<T;i++){ pthread_join(th[i],0); cnt+=jb[i].cnt; }
      double el=now()-t0; if(el<best) best=el; }
    printf("mm_sketch %s k=%d w=%d T=%d n=%ld out=%ld time=%.4fs %.3f ns/bp %.3f Gbp/s\n", argv[1], K, W, T, n, cnt, best, best*1e9/n, n/best/1e9);
  }
  return 0;
}
