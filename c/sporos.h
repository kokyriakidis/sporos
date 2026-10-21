/*
 * sporos from C and C++: (w, k)-minimizers, closed and open syncmers, the maximal exact matches
 * (MEMs), maximal unique matches (MUMs) and frequency-capped MEMs of a query against an indexed
 * reference, the super-maximal exact matches (SMEMs) of reads against an FM-index, and seed-index
 * lookups. Build the library
 * with `pixi run build-c [target-cpu]`, which leaves it in build/c beside the Mojo runtime libraries it
 * loads and this header; link with `-Lbuild/c -lsporos` and put build/c on the program's library
 * search path (an rpath, say). Built for an explicit `target-cpu`, the library runs only on CPUs with
 * that instruction set; the default is the platform's oldest (x86-64, apple-m1).
 *
 * Sequences are ASCII bytes: A, C, G, T and U in either case are bases, anything else (N, say) is a
 * base no seed covers and no match crosses. Positions are 0-based. `threads` of zero or less uses every
 * physical core; sporos runs its own threads and keeps no state between calls, so any function may be
 * called from many threads at once, an index queried by many at once. `max_occurrences` and
 * `min_occurrences` count a matched string's occurrences in the reference throughout, as in the Mojo
 * library.
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

/* A file that could not be written. */
#define SPOROS_IO_ERROR (-2)

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
 * positions count along it, as MUMmer reports reverse matches.
 *
 * `max_occurrences` N > 0 relaxes the modes' "once" to "at most N times" (Mumemto's -F; for
 * SPOROS_MAXMATCH it caps the reference count), and `min_occurrences` (1 for none) keeps matches whose
 * string occurs at least that often in the reference. Null for a bad mode, a negative length, a
 * query of 2^32 bases or more, or limits that admit nothing. Free with sporos_matches_free.
 */
sporos_matches *sporos_index_matches(const sporos_index *index, const uint8_t *query, int64_t qlen, int32_t mode,
                                     int32_t reverse_complement, int32_t max_occurrences, int32_t min_occurrences,
                                     int32_t threads);

/* How many matches there are. */
int64_t sporos_matches_count(const sporos_matches *matches);

/* The matches, sporos_matches_count of them, valid until sporos_matches_free. */
const sporos_match *sporos_matches_data(const sporos_matches *matches);

/* Frees matches; null is nothing to free. */
void sporos_matches_free(sporos_matches *matches);

/* An FM-index of reference sequences, for the SMEMs of many reads. */
typedef struct sporos_fm_index sporos_fm_index;

/*
 * Indexes `n_seqs` reference sequences, `seqs[i]` of `lens[i]` bytes, named `names[i]` (NUL-terminated;
 * `names` null names them seq0, seq1, ...); the sequences are copied. Null for no sequences or a
 * negative length. Free with sporos_fm_index_free.
 */
sporos_fm_index *sporos_fm_index_new(const uint8_t *const *seqs, const int64_t *lens, const char *const *names,
                                     int64_t n_seqs, int32_t threads);

/* Reads an index sporos_fm_index_save (or `sporos index`) wrote; null when it cannot. */
sporos_fm_index *sporos_fm_index_load(const char *path, int32_t threads);

/* Writes an index to `path`: 0, or SPOROS_IO_ERROR. */
int32_t sporos_fm_index_save(const sporos_fm_index *index, const char *path);

/* Frees an FM-index; null is nothing to free. */
void sporos_fm_index_free(sporos_fm_index *index);

/* How many reference sequences an FM-index holds. */
int64_t sporos_fm_index_sequences(const sporos_fm_index *index);

/* Sequence `i`'s name, valid while the index is; null for no such sequence. */
const char *sporos_fm_index_name(const sporos_fm_index *index, int64_t i);

/* One SMEM: query bases [start, end) occur `count` times in the reference, on either strand. */
typedef struct {
    uint32_t start;    /* 0-based query start. */
    uint32_t end;      /* One past the last base. */
    int64_t count;     /* Reference occurrences, both strands. */
    int64_t first_hit; /* Its hits in sporos_smems_hits: `hit_count` from here. */
    int64_t hit_count; /* `count`, or 0 when that is above `max_hits`. */
} sporos_smem;

/* Where an SMEM occurs. */
typedef struct {
    uint32_t ref_seq; /* Index of the reference sequence. */
    uint32_t ref_pos; /* 0-based position of the match's leftmost base in it. */
    int32_t reverse;  /* Nonzero when the query matches the reverse strand. */
} sporos_smem_hit;

/* Many queries' SMEMs. */
typedef struct sporos_smems sporos_smems;

/*
 * The SMEMs of length >= `min_len` of `n_queries` queries, both strands: maximal query intervals that
 * occur in the reference and lie in no longer one that does, in query order, with up to `max_hits`
 * hits each (an SMEM with more lists none). With `min_occurrences` f > 1 they are f-SMEMs, the
 * maximal intervals occurring at least f times (minibwa's -s); `max_occurrences` > 0 drops SMEMs
 * occurring more often. Null for bad arguments. Free with sporos_smems_free.
 */
sporos_smems *sporos_fm_index_smems(const sporos_fm_index *index, const uint8_t *const *queries,
                                    const int64_t *lens, int64_t n_queries, int32_t min_len, int32_t max_hits,
                                    int32_t min_occurrences, int32_t max_occurrences, int32_t threads);

/* How many SMEMs query `query` has, and the SMEMs, valid until sporos_smems_free. */
int64_t sporos_smems_count(const sporos_smems *smems, int64_t query);
const sporos_smem *sporos_smems_data(const sporos_smems *smems, int64_t query);

/* Every SMEM's hits (see sporos_smem's `first_hit`). */
const sporos_smem_hit *sporos_smems_hits(const sporos_smems *smems);

/* Frees SMEMs; null is nothing to free. */
void sporos_smems_free(sporos_smems *smems);

/* Which reference k-mers a seed index holds. */
#define SPOROS_MINIMIZER_SEEDS 0 /* (w, k)-minimizers. */
#define SPOROS_SYNCMER_SEEDS 1   /* Syncmers, with `w` the s-mer length. */
#define SPOROS_ALL_KMERS 2       /* Every k-mer. */

/* A hash index of reference k-mers (k <= 32), probed with reads' seeds. */
typedef struct sporos_seed_index sporos_seed_index;

/*
 * Indexes the seeds of `n_seqs` reference sequences under `scheme`: syncmers are closed for a negative
 * `open_offset`, else open at it. Null for no sequences, a negative length, bad seed parameters or
 * 2^32 bases or more in all. Free with sporos_seed_index_free.
 */
sporos_seed_index *sporos_seed_index_new(const uint8_t *const *seqs, const int64_t *lens, int64_t n_seqs, int32_t k,
                                         int32_t w, int32_t scheme, int32_t open_offset, int32_t threads);

/* Frees a seed index; null is nothing to free. */
void sporos_seed_index_free(sporos_seed_index *index);

/* The smallest occurrence count that `keep_fraction` of the indexed k-mers stay within (a
 * `max_occurrences` for lookups, like GraphAligner's quantile cutoff). */
int64_t sporos_seed_index_occurrence_cutoff(const sporos_seed_index *index, double keep_fraction);

/* A read k-mer found in the reference. */
typedef struct {
    uint32_t query_pos;   /* The k-mer's 0-based start in the read. */
    uint32_t ref_seq;     /* Index of the reference sequence. */
    uint32_t ref_pos;     /* The reference k-mer's 0-based start in it. */
    uint32_t occurrences; /* The k-mer's reference occurrences, both strands. */
    int32_t reverse;      /* Nonzero when it matches the reference's reverse strand. */
} sporos_seed_hit;

/* Many reads' seed hits. */
typedef struct sporos_seed_hits sporos_seed_hits;

/*
 * The hits of `n_reads` reads' seeds (of every k-mer, with `all_kmers` nonzero), skipping k-mers that
 * occur more than `max_occurrences` times (0 for no limit). Null for bad arguments. Free with
 * sporos_seed_hits_free.
 */
sporos_seed_hits *sporos_seed_index_lookup(const sporos_seed_index *index, const uint8_t *const *reads,
                                           const int64_t *lens, int64_t n_reads, int32_t all_kmers,
                                           int64_t max_occurrences, int32_t threads);

/* How many hits read `read` has, and the hits in read order, valid until sporos_seed_hits_free. */
int64_t sporos_seed_hits_count(const sporos_seed_hits *hits, int64_t read);
const sporos_seed_hit *sporos_seed_hits_data(const sporos_seed_hits *hits, int64_t read);

/* Frees seed hits; null is nothing to free. */
void sporos_seed_hits_free(sporos_seed_hits *hits);

#ifdef __cplusplus
}
#endif

#endif /* SPOROS_H */
