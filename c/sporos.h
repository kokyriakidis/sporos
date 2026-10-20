/*
 * sporos from C and C++: (w, k)-minimizers, closed and open syncmers, and the maximal exact matches
 * (MEMs) and maximal unique matches (MUMs) of a query against an indexed reference. Build the library
 * with `pixi run build-c [target-cpu]`, which leaves it in build/c beside the Mojo runtime libraries it
 * loads and this header; link with `-Lbuild/c -lsporos` and put build/c on the program's library
 * search path (an rpath, say). Built for an explicit `target-cpu`, the library runs only on CPUs with
 * that instruction set; the default is the platform's oldest (x86-64, apple-m1).
 *
 * Sequences are ASCII bytes: A, C, G, T and U in either case are bases, anything else (N, say) is a
 * base no seed covers and no match crosses. Positions are 0-based. `threads` of zero or less uses every
 * physical core; sporos runs its own threads and keeps no state between calls, so any function may be
 * called from many threads at once, an index queried by many at once.
 */
#ifndef SPOROS_H
#define SPOROS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* A negative length, or seed parameters out of range: k in [1, 64] and w >= 1 for minimizers,
 * 1 <= s <= k and an open offset in [0, k - s] for syncmers. */
#define SPOROS_INVALID_ARGUMENT (-1)

/*
 * The start positions of the (w, k)-minimizers of `seq`, in increasing order, each once per run of
 * windows sharing it; canonical (strand-independent) when `canonical` is nonzero. Writes the first
 * `out_capacity` into `out` and returns how many there are, which may be more (`len` is always
 * enough room), or SPOROS_INVALID_ARGUMENT.
 */
int64_t sporos_minimizers(const uint8_t *seq, int64_t len, int32_t k, int32_t w, int32_t canonical,
                          int32_t threads, uint32_t *out, int64_t out_capacity);

/*
 * The start positions of the syncmers of `seq`: k-mers whose smallest s-mer starts at offset 0 or
 * k - s (closed, `open_offset` negative), or at `open_offset` (open). Otherwise as sporos_minimizers.
 */
int64_t sporos_syncmers(const uint8_t *seq, int64_t len, int32_t k, int32_t s, int32_t open_offset,
                        int32_t canonical, int32_t threads, uint32_t *out, int64_t out_capacity);

/*
 * Open syncmers with several allowed offsets: k-mers whose smallest s-mer starts at any of the
 * `offset_count` (>= 1) offsets in [0, k - s]. Otherwise as sporos_syncmers.
 */
int64_t sporos_syncmers_open(const uint8_t *seq, int64_t len, int32_t k, int32_t s, const int32_t *offsets,
                             int64_t offset_count, int32_t canonical, int32_t threads, uint32_t *out,
                             int64_t out_capacity);

/* A reference index for matches of at least a minimum length. */
typedef struct sporos_index sporos_index;

/*
 * Indexes `n_seqs` reference sequences, `seqs[i]` of `lens[i]` bytes, for matches of at least
 * `min_len` bases; the sequences are copied, so they may be freed afterwards. Null for no sequences,
 * a negative length, `min_len` below 1, or 2^32 bases or more in all. Free with sporos_index_free.
 */
sporos_index *sporos_index_new(const uint8_t *const *seqs, const int64_t *lens, int64_t n_seqs, int32_t min_len,
                               int32_t threads);

/* Frees an index; null is nothing to free. */
void sporos_index_free(sporos_index *index);

/* How many reference sequences an index holds. */
int64_t sporos_index_sequences(const sporos_index *index);

/* Which matches sporos_index_matches reports. */
#define SPOROS_MAXMATCH 0     /* Every MEM, MUMmer's -maxmatch. */
#define SPOROS_MUM 1          /* MEMs unique in the reference and the query, MUMmer's -mum. */
#define SPOROS_MUMREFERENCE 2 /* MEMs unique in the reference, MUMmer's -mumreference. */

/* One match: `length` bases of reference sequence `ref_seq` from `ref_pos` equal the query's from
 * `query_pos`. */
typedef struct {
    uint32_t ref_seq;   /* Index of the reference sequence, in the order given to sporos_index_new. */
    uint32_t ref_pos;   /* 0-based start within that sequence. */
    uint32_t query_pos; /* 0-based start within the query, or within its reverse complement. */
    uint32_t length;    /* Bases matched, at least the index's `min_len`. */
} sporos_match;

/* One query's matches. */
typedef struct sporos_matches sporos_matches;

/*
 * The matches of `query` against `index` under `mode` (SPOROS_MAXMATCH, SPOROS_MUM or
 * SPOROS_MUMREFERENCE), sorted by query position, then reference position. With
 * `reverse_complement` nonzero the query's reverse complement is matched instead, and query
 * positions count along it, as MUMmer reports reverse matches. Null for a bad mode, a negative
 * length or a query of 2^32 bases or more. Free with sporos_matches_free.
 */
sporos_matches *sporos_index_matches(const sporos_index *index, const uint8_t *query, int64_t qlen, int32_t mode,
                                     int32_t reverse_complement, int32_t threads);

/* How many matches there are. */
int64_t sporos_matches_count(const sporos_matches *matches);

/* The matches, sporos_matches_count of them, valid until sporos_matches_free. */
const sporos_match *sporos_matches_data(const sporos_matches *matches);

/* Frees matches; null is nothing to free. */
void sporos_matches_free(sporos_matches *matches);

#ifdef __cplusplus
}
#endif

#endif /* SPOROS_H */
