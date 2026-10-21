// minimap2 seed lookups: index a reference, then for every read compute its
// minimizers (mm_sketch) and fetch every reference hit of each (mm_idx_get),
// skipping minimizers with more than max_occ hits. Times the read loop.
//   mm2_lookup <ref.fa> <reads.fa> k w max_occ
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
#include "minimap.h"
#include "mmpriv.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct { uint32_t qpos, ref, rpos, strand, occ; } hit_t;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }

int main(int argc, char *argv[]) {
	int k = atoi(argv[3]), w = atoi(argv[4]), max_occ = atoi(argv[5]);
	mm_idxopt_t iopt; mm_mapopt_t mopt;
	mm_verbose = 1;
	mm_set_opt(0, &iopt, &mopt);
	iopt.k = k, iopt.w = w, iopt.flag = 0, iopt.batch_size = 0x7fffffffffffffffL;
	for (int threads = 1; threads <= 8; threads += 7) {
		double t0 = now();
		mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, 0);
		mm_idx_t *mi = mm_idx_reader_read(r, threads);
		printf("build threads=%d: %.3f s\n", threads, now() - t0);
		mm_idx_reader_close(r);
		if (threads == 8) { mm_idx_destroy(mi); break; }
		// reads into memory
		gzFile fp = gzopen(argv[2], "r"); kseq_t *ks = kseq_init(fp);
		int n = 0, m = 0; char **seqs = 0; int *lens = 0; int64_t bases = 0;
		while (kseq_read(ks) >= 0) {
			if (n == m) { m = m ? m * 2 : 1024; seqs = realloc(seqs, m * sizeof(char*)); lens = realloc(lens, m * sizeof(int)); }
			seqs[n] = malloc(ks->seq.l + 1); memcpy(seqs[n], ks->seq.s, ks->seq.l + 1); lens[n] = ks->seq.l; bases += ks->seq.l; n++;
		}
		kseq_destroy(ks); gzclose(fp);
		hit_t *hits = malloc(sizeof(hit_t) * 1 << 20); size_t cap = 1 << 20;
		double best = 1e30; size_t total = 0;
		for (int rep = 0; rep < 3; rep++) {
			mm128_v a = {0, 0, 0}; total = 0;
			double t = now();
			for (int i = 0; i < n; i++) {
				a.n = 0;
				mm_sketch(0, seqs[i], lens[i], w, k, 0, 0, &a);
				size_t nh = 0;
				for (size_t j = 0; j < a.n; j++) {
					int cnt;
					const uint64_t *cr = mm_idx_get(mi, a.a[j].x >> 8, &cnt);
					if (cnt == 0 || (max_occ > 0 && cnt > max_occ)) continue;
					if (nh + cnt > cap) { cap = (nh + cnt) * 2; hits = realloc(hits, sizeof(hit_t) * cap); }
					uint32_t qpos = (uint32_t)a.a[j].y >> 1;
					for (int c = 0; c < cnt; c++) {
						hit_t *h = &hits[nh++];
						h->qpos = qpos, h->ref = cr[c] >> 32, h->rpos = (uint32_t)cr[c] >> 1;
						h->strand = (cr[c] ^ a.a[j].y) & 1, h->occ = cnt;
					}
				}
				total += nh;
			}
			double e = now() - t;
			if (e < best) best = e;
			free(a.a);
		}
		printf("lookup minimizers, per read: %.1f Mbp/s, %.2f hits/read (%.3f s)\n", bases / best / 1e6, (double)total / n, best);
		mm_idx_destroy(mi);
	}
	return 0;
}
