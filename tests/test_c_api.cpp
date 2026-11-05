// The C API (see c/sporos.h): seeds, matches, SMEMs and seed-index hits held to what the Mojo library
// itself answers on the same inputs (tests/api_expected.mojo writes those), matches held to a
// brute-force oracle sharing no code with the library, known small answers, refused arguments,
// truncated output, saved and loaded FM-indexes and four threads querying one index at once.
//
//     pixi run test-c
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "sporos.h"

// The checks failed so far; any one fails the run, after every check has been reported.
static int failures = 0;

// Reports `condition` with its place when it fails, and counts it, without stopping the run.
#define CHECK(condition)                                                                    \
    do {                                                                                    \
        if (!(condition)) {                                                                 \
            std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition);    \
            ++failures;                                                                     \
        }                                                                                   \
    } while (0)

using Match = std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>;  // ref_seq, ref_pos, query_pos, length

static const uint8_t *bytes(const std::string &text) { return reinterpret_cast<const uint8_t *>(text.data()); }

static std::vector<uint32_t> minimizers(const std::string &seq, int k, int w, int canonical, int threads) {
    std::vector<uint32_t> out(seq.size());
    int64_t count = sporos_minimizers(bytes(seq), seq.size(), k, w, canonical, threads, out.data(), out.size());
    CHECK(count >= 0 && count <= int64_t(seq.size()));
    out.resize(std::max<int64_t>(count, 0));
    return out;
}

static std::vector<uint32_t> syncmers(const std::string &seq, int k, int s, int open, int canonical, int threads) {
    std::vector<uint32_t> out(seq.size());
    int64_t count = sporos_syncmers(bytes(seq), seq.size(), k, s, open, canonical, threads, out.data(), out.size());
    CHECK(count >= 0 && count <= int64_t(seq.size()));
    out.resize(std::max<int64_t>(count, 0));
    return out;
}

// Open syncmers at any of several offsets.
static std::vector<uint32_t> syncmers_open(const std::string &seq, int k, int s, const std::vector<int32_t> &offsets,
                                           int canonical, int threads) {
    std::vector<uint32_t> out(seq.size());
    int64_t count = sporos_syncmers_open(bytes(seq), seq.size(), k, s, offsets.data(), offsets.size(), canonical,
                                         threads, out.data(), out.size());
    CHECK(count >= 0 && count <= int64_t(seq.size()));
    out.resize(std::max<int64_t>(count, 0));
    return out;
}

static sporos_index *index_of(const std::vector<std::string> &refs, int min_len, int threads) {
    std::vector<const uint8_t *> seqs;
    std::vector<int64_t> lens;
    for (const std::string &ref : refs) seqs.push_back(bytes(ref)), lens.push_back(ref.size());
    return sporos_index_new(seqs.data(), lens.data(), refs.size(), min_len, threads);
}

// Pointers and lengths of many sequences, for the calls taking several.
struct Many {
    std::vector<const uint8_t *> seqs;
    std::vector<int64_t> lens;
    explicit Many(const std::vector<std::string> &texts) {
        for (const std::string &text : texts) seqs.push_back(bytes(text)), lens.push_back(text.size());
    }
};

using Hit = std::tuple<uint32_t, uint32_t, int>;                         // ref_seq, ref_pos, reverse
using Smem = std::tuple<uint32_t, uint32_t, int64_t, std::vector<Hit>>;  // start, end, count, hits

static std::vector<std::vector<Smem>> smems(const sporos_fm_index *index, const std::vector<std::string> &reads,
                                            int min_len, int max_hits, int min_occ, int max_occ, int threads) {
    Many many(reads);
    sporos_smems *found = sporos_fm_index_smems(index, many.seqs.data(), many.lens.data(), reads.size(), min_len,
                                                max_hits, min_occ, max_occ, threads);
    CHECK(found != nullptr);
    std::vector<std::vector<Smem>> out(reads.size());
    if (!found) return out;
    const sporos_smem_hit *hits = sporos_smems_hits(found);
    for (size_t q = 0; q < reads.size(); ++q) {
        const sporos_smem *data = sporos_smems_data(found, q);
        for (int64_t i = 0; i < sporos_smems_count(found, q); ++i) {
            std::vector<Hit> listed;
            for (int64_t h = data[i].first_hit; h < data[i].first_hit + data[i].hit_count; ++h)
                listed.emplace_back(hits[h].ref_seq, hits[h].ref_pos, hits[h].reverse != 0);
            out[q].emplace_back(data[i].start, data[i].end, data[i].count, listed);
        }
    }
    sporos_smems_free(found);
    return out;
}

using SeedHit = std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>;  // query_pos, ref_seq, ref_pos, occ, reverse

static std::vector<std::vector<SeedHit>> lookup(const sporos_seed_index *index, const std::vector<std::string> &reads,
                                                int all_kmers, int64_t max_occ, int threads) {
    Many many(reads);
    sporos_seed_hits *found =
        sporos_seed_index_lookup(index, many.seqs.data(), many.lens.data(), reads.size(), all_kmers, max_occ, threads);
    CHECK(found != nullptr);
    std::vector<std::vector<SeedHit>> out(reads.size());
    if (!found) return out;
    for (size_t r = 0; r < reads.size(); ++r) {
        const sporos_seed_hit *data = sporos_seed_hits_data(found, r);
        for (int64_t i = 0; i < sporos_seed_hits_count(found, r); ++i)
            out[r].emplace_back(data[i].query_pos, data[i].ref_seq, data[i].ref_pos, data[i].occurrences,
                                data[i].reverse != 0);
    }
    sporos_seed_hits_free(found);
    return out;
}

static std::vector<Match> matches(const sporos_index *index, const std::string &query, int mode, int reverse,
                                  int threads, int max_occurrences = 0, int min_occurrences = 1) {
    sporos_matches *found = sporos_index_matches(index, bytes(query), query.size(), mode, reverse, max_occurrences,
                                                 min_occurrences, threads);
    CHECK(found != nullptr);
    std::vector<Match> out;
    if (!found) return out;
    const sporos_match *data = sporos_matches_data(found);
    for (int64_t i = 0; i < sporos_matches_count(found); ++i)
        out.emplace_back(data[i].ref_seq, data[i].ref_pos, data[i].query_pos, data[i].length);
    sporos_matches_free(found);
    return out;
}

// A base's 2-bit identity, -1 for anything no match covers.
static int code(char c) {
    switch (c | 0x20) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': case 'u': return 3;
        default: return -1;
    }
}

static bool same(char a, char b) { return code(a) >= 0 && code(a) == code(b); }

static std::string reverse_complement(const std::string &seq) {
    std::string out;
    for (size_t i = seq.size(); i-- > 0;) out += code(seq[i]) < 0 ? 'N' : "TGCA"[code(seq[i])];
    return out;
}

// How often `query[start, start + length)` occurs in `text`.
static int occurrences(const std::string &text, const std::string &query, size_t start, size_t length) {
    int count = 0;
    for (size_t i = 0; i + length <= text.size(); ++i) {
        size_t j = 0;
        while (j < length && same(text[i + j], query[start + j])) ++j;
        count += j == length;
    }
    return count;
}

// Every MEM of at least `min_len` by trying every pair of starts, kept by `mode` and the occurrence
// limits from how often its string occurs in the references and the query, sorted as the library sorts them.
static std::vector<Match> oracle(const std::vector<std::string> &refs, std::string query, size_t min_len, int mode,
                                 bool reverse, int max_occ = 0, int min_occ = 1) {
    int cap = max_occ > 0 ? max_occ : mode == SPOROS_MAXMATCH ? 0 : 1;
    if (reverse) query = reverse_complement(query);
    std::vector<Match> out;
    for (size_t r = 0; r < refs.size(); ++r) {
        const std::string &ref = refs[r];
        for (size_t i = 0; i < ref.size(); ++i)
            for (size_t j = 0; j < query.size(); ++j) {
                if (i > 0 && j > 0 && same(ref[i - 1], query[j - 1])) continue;
                size_t length = 0;
                while (i + length < ref.size() && j + length < query.size() && same(ref[i + length], query[j + length]))
                    ++length;
                if (length < min_len) continue;
                if (cap > 0 || min_occ > 1) {
                    int in_refs = 0;
                    for (const std::string &other : refs) in_refs += occurrences(other, query, j, length);
                    if ((cap > 0 && in_refs > cap) || in_refs < min_occ) continue;
                    if (mode == SPOROS_MUM && occurrences(query, query, j, length) > cap) continue;
                }
                out.emplace_back(r, i, j, length);
            }
    }
    std::sort(out.begin(), out.end(), [](const Match &a, const Match &b) {
        return std::make_tuple(std::get<2>(a), std::get<0>(a), std::get<1>(a), std::get<3>(a)) <
               std::make_tuple(std::get<2>(b), std::get<0>(b), std::get<1>(b), std::get<3>(b));
    });
    return out;
}

static std::string random_dna(size_t length, std::mt19937 &random) {
    std::uniform_int_distribution<int> base(0, 7);
    std::string out;
    for (size_t i = 0; i < length; ++i) out += "ACGTacgt"[base(random)];
    return out;
}

// The known answers: every k-mer is its own window's minimizer at w = 1 and a closed syncmer at s = k.
static void known() {
    std::vector<uint32_t> all = {0, 1, 2, 3, 4, 5};
    CHECK(minimizers("ACGTACGTAC", 5, 1, 1, 1) == all);
    CHECK(syncmers("ACGTACGTAC", 5, 5, -1, 1, 1) == all);
    // No k-mer covers an N.
    CHECK(minimizers("ACGTNACGTAC", 4, 1, 0, 1) == std::vector<uint32_t>({0, 5, 6, 7}));
    CHECK(syncmers("ACGTNACGTAC", 4, 4, -1, 0, 1) == std::vector<uint32_t>({0, 5, 6, 7}));
    CHECK(minimizers("ACG", 5, 3, 1, 1).empty());
    CHECK(sporos_minimizers(nullptr, 0, 5, 3, 1, 1, nullptr, 0) == 0);
    CHECK(sporos_syncmers(nullptr, 0, 5, 3, -1, 1, 1, nullptr, 0) == 0);

    // Refused: k outside [1, 64], w below 1, s outside [1, k], an open offset past k - s, a negative length.
    uint32_t slot[16];
    const std::string text = "ACGTACGTACGTACGT";
    const uint8_t *seq = bytes(text);
    CHECK(sporos_minimizers(seq, 16, 0, 3, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_minimizers(seq, 16, 65, 3, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_minimizers(seq, 16, 5, 0, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_minimizers(seq, -1, 5, 3, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_syncmers(seq, 16, 5, 6, -1, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_syncmers(seq, 16, 5, 0, -1, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_syncmers(seq, 16, 5, 3, 3, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_syncmers(seq, -1, 5, 3, -1, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    const int32_t offsets[] = {0, 3};
    CHECK(sporos_syncmers_open(seq, 16, 5, 3, offsets, 2, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);
    CHECK(sporos_syncmers_open(seq, 16, 5, 3, offsets, 0, 1, 1, slot, 16) == SPOROS_INVALID_ARGUMENT);

    // Too little room: the full count, and nothing written past the room given.
    std::vector<uint32_t> room(8, UINT32_MAX);
    CHECK(sporos_minimizers(seq, 16, 4, 1, 0, 1, room.data(), 3) == 13);
    CHECK(room[0] == 0 && room[1] == 1 && room[2] == 2 && room[3] == UINT32_MAX);

    // An exact repeat in the query is a MEM but no MUM; reverse matches count along the reverse complement.
    std::vector<std::string> refs = {"TTTTTTTTACGTACCGTAGGCATCAATTTTTTTT", "GGGGGGGCCGATTACAGGATCAGGGGGGGG"};
    sporos_index *index = index_of(refs, 10, 1);
    CHECK(index != nullptr);
    CHECK(sporos_index_sequences(index) == 2);
    std::string query = "ACGTACCGTAGGCATCAA";
    CHECK(matches(index, query, SPOROS_MUM, 0, 1) == std::vector<Match>({Match(0, 8, 0, 18)}));
    CHECK(matches(index, query + "CC" + query, SPOROS_MUM, 0, 1).empty());
    CHECK(matches(index, query + "CC" + query, SPOROS_MUMREFERENCE, 0, 1).size() == 2);
    CHECK(matches(index, reverse_complement("CCGATTACAGGATC"), SPOROS_MUM, 1, 1) ==
          std::vector<Match>({Match(1, 7, 0, 14)}));
    CHECK(matches(index, "", SPOROS_MAXMATCH, 0, 1).empty());
    CHECK(sporos_index_matches(index, bytes(query), query.size(), 3, 0, 0, 1, 1) == nullptr);
    CHECK(sporos_index_matches(index, bytes(query), -1, 0, 0, 0, 1, 1) == nullptr);
    // Limits: a repeat twice in the query passes a query cap of 2; a floor above the cap, or below 1, is refused.
    CHECK(matches(index, query + "CC" + query, SPOROS_MUM, 0, 1, 2).size() == 2);
    CHECK(matches(index, query, SPOROS_MAXMATCH, 0, 1, 0, 2).empty());
    CHECK(sporos_index_matches(index, bytes(query), query.size(), SPOROS_MUM, 0, 0, 2, 1) == nullptr);
    CHECK(sporos_index_matches(index, bytes(query), query.size(), SPOROS_MAXMATCH, 0, 0, 0, 1) == nullptr);
    sporos_index_free(index);

    // SMEMs: the query sits once in the first reference, named as given.
    const char *names[] = {"chrA", "chrB"};
    Many many(refs);
    sporos_fm_index *fm = sporos_fm_index_new(many.seqs.data(), many.lens.data(), names, 2, 1);
    CHECK(fm != nullptr && sporos_fm_index_sequences(fm) == 2);
    CHECK(std::string(sporos_fm_index_name(fm, 1)) == "chrB");
    CHECK(sporos_fm_index_name(fm, 2) == nullptr);
    std::vector<std::vector<Smem>> found = smems(fm, {"GG" + query + "GG"}, 10, 20, 1, 0, 1);
    CHECK(found.size() == 1 && found[0].size() == 1);
    CHECK(found[0][0] == Smem(2, 20, 1, {Hit(0, 8, 0)}));
    CHECK(smems(fm, {reverse_complement(query)}, 10, 20, 1, 0, 1)[0] == std::vector<Smem>({Smem(0, 18, 1, {Hit(0, 8, 1)})}));
    CHECK(smems(fm, {query}, 10, 20, 2, 0, 1)[0].empty());  // nothing occurs twice
    CHECK(smems(fm, {}, 10, 20, 1, 0, 1).empty());
    const uint8_t *q = bytes(query);
    int64_t qlen = query.size(), negative_length = -1;
    CHECK(sporos_fm_index_smems(fm, &q, &qlen, 1, 0, 20, 1, 0, 1) == nullptr);
    CHECK(sporos_fm_index_smems(fm, &q, &negative_length, 1, 10, 20, 1, 0, 1) == nullptr);
    CHECK(sporos_fm_index_save(fm, "/nonexistent/dir/x.sfm") == SPOROS_IO_ERROR);
    CHECK(sporos_fm_index_load("/nonexistent/dir/x.sfm", 1) == nullptr);
    sporos_fm_index_free(fm);
    sporos_fm_index *unnamed = sporos_fm_index_new(many.seqs.data(), many.lens.data(), nullptr, 2, 1);
    CHECK(unnamed != nullptr && std::string(sporos_fm_index_name(unnamed, 0)) == "seq0");
    sporos_fm_index_free(unnamed);
    CHECK(sporos_fm_index_new(many.seqs.data(), many.lens.data(), nullptr, 0, 1) == nullptr);

    // Seed index: every 5-mer of the first reference (ACGTA also occurs as TACGT at 7), refused parameters.
    sporos_seed_index *seeds = sporos_seed_index_new(many.seqs.data(), many.lens.data(), 1, 5, 1, SPOROS_ALL_KMERS, -1, 1);
    CHECK(seeds != nullptr);
    std::vector<std::vector<SeedHit>> hits = lookup(seeds, {"ACGTACC"}, 0, 0, 1);
    CHECK(hits.size() == 1 &&
          hits[0] == std::vector<SeedHit>({SeedHit(0, 0, 8, 2, 0), SeedHit(0, 0, 7, 2, 1), SeedHit(1, 0, 9, 1, 0),
                                           SeedHit(2, 0, 10, 1, 0)}));
    sporos_seed_index_free(seeds);
    CHECK(sporos_seed_index_new(many.seqs.data(), many.lens.data(), 1, 65, 1, SPOROS_ALL_KMERS, -1, 1) == nullptr);
    CHECK(sporos_seed_index_new(many.seqs.data(), many.lens.data(), 1, 5, 1, 3, -1, 1) == nullptr);
    sporos_seed_index_free(nullptr);
    sporos_smems_free(nullptr);
    sporos_seed_hits_free(nullptr);
    sporos_fm_index_free(nullptr);

    CHECK(index_of({}, 10, 1) == nullptr);
    CHECK(index_of(refs, 0, 1) == nullptr);
    int64_t negative = -1;
    const uint8_t *one = bytes(refs[0]);
    CHECK(sporos_index_new(&one, &negative, 1, 10, 1) == nullptr);
    sporos_index_free(nullptr);
    sporos_matches_free(nullptr);
}

// Matches on random references and a query built from their pieces, against the oracle.
static void random_matches() {
    std::mt19937 random(7);
    for (int round = 0; round < 6; ++round) {
        std::vector<std::string> refs = {random_dna(300, random), random_dna(250, random)};
        refs[1][100] = 'N';
        refs[1].replace(200, 40, refs[0].substr(30, 40));
        std::string query = refs[0].substr(20, 60) + random_dna(30, random) +
                            reverse_complement(refs[1].substr(50, 80)) + refs[0].substr(40, 30) +
                            refs[1].substr(150, 70);
        query[90 + round] = query[90 + round] == 'A' ? 'C' : 'A';
        size_t min_len = 10 + 3 * round;
        sporos_index *index = index_of(refs, min_len, round % 3);
        for (int mode = 0; mode < 3; ++mode)
            for (int reverse = 0; reverse < 2; ++reverse) {
                std::vector<Match> got = matches(index, query, mode, reverse, 1 + round % 2);
                std::vector<Match> expected = oracle(refs, query, min_len, mode, reverse);
                CHECK(got == expected);
                if (got != expected)
                    std::fprintf(stderr, "  round %d mode %d reverse %d: %zu matches, expected %zu\n", round, mode,
                                 reverse, got.size(), expected.size());
                // Occurrence limits on a repeat planted in both references.
                CHECK(matches(index, query, mode, reverse, 1, 2) == oracle(refs, query, min_len, mode, reverse, 2));
                CHECK(matches(index, query, mode, reverse, 1, 3, 2) ==
                      oracle(refs, query, min_len, mode, reverse, 3, 2));
            }
        sporos_index_free(index);
    }
}

// Every answer the Mojo library wrote, through the C API.
static void against_mojo(const char *path) {
    std::ifstream file(path);
    CHECK(file.good());
    std::map<std::string, std::string> seqs;
    std::map<int, sporos_index *> indexes;
    sporos_fm_index *fm = nullptr;
    std::string line;
    int records = 0;
    while (std::getline(file, line)) {
        std::istringstream fields(line);
        std::string kind, name;
        fields >> kind >> name;
        if (kind == "seq") {
            fields >> seqs[name];
            continue;
        }
        ++records;
        const std::string &seq = seqs[name];
        if (kind == "minimizers" || kind == "syncmers") {
            int k, second, canonical, threads;
            std::vector<int32_t> open;
            fields >> k >> second;
            if (kind == "syncmers") {
                std::string list, offset;
                fields >> list;
                std::istringstream offsets(list);
                while (std::getline(offsets, offset, ',')) open.push_back(std::stoi(offset));
            }
            fields >> canonical >> threads;
            size_t count;
            fields >> count;
            std::vector<uint32_t> expected(count);
            for (uint32_t &position : expected) fields >> position;
            std::vector<uint32_t> got = kind == "minimizers" ? minimizers(seq, k, second, canonical, threads)
                                     : open.size() > 1 ? syncmers_open(seq, k, second, open, canonical, threads)
                                                       : syncmers(seq, k, second, open[0], canonical, threads);
            CHECK(got == expected);
            if (got != expected) std::fprintf(stderr, "  %s\n", line.substr(0, 60).c_str());
            continue;
        }
        std::vector<std::string> refs = {seqs["ref0"], seqs["ref1"]};
        if (kind == "smems" || kind == "seeds") {
            size_t n;
            fields >> n;
            std::vector<std::string> reads;
            for (size_t r = 0; r < n; ++r) reads.push_back(seqs[name + std::to_string(r)]);
            bool same_answers = true;
            if (kind == "smems") {
                int min_len, max_hits, min_occ, max_occ, threads;
                fields >> min_len >> max_hits >> min_occ >> max_occ >> threads;
                if (!fm) {
                    Many many(refs);
                    const char *names[] = {"ref0", "ref1"};
                    fm = sporos_fm_index_new(many.seqs.data(), many.lens.data(), names, 2, 2);
                }
                std::vector<std::vector<Smem>> got = smems(fm, reads, min_len, max_hits, min_occ, max_occ, threads);
                for (size_t r = 0; r < n; ++r) {
                    size_t count;
                    fields >> count;
                    std::vector<Smem> expected(count);
                    for (Smem &smem : expected) {
                        size_t listed;
                        fields >> std::get<0>(smem) >> std::get<1>(smem) >> std::get<2>(smem) >> listed;
                        std::get<3>(smem).resize(listed);
                        for (Hit &hit : std::get<3>(smem)) fields >> std::get<0>(hit) >> std::get<1>(hit) >> std::get<2>(hit);
                    }
                    same_answers &= got[r] == expected;
                }
            } else {
                int k, w, scheme, open, all_kmers, threads;
                int64_t max_occ;
                fields >> k >> w >> scheme >> open >> max_occ >> all_kmers >> threads;
                Many many(refs);
                sporos_seed_index *index =
                    sporos_seed_index_new(many.seqs.data(), many.lens.data(), 2, k, w, scheme, open, threads);
                CHECK(index != nullptr);
                std::vector<std::vector<SeedHit>> got = lookup(index, reads, all_kmers, max_occ, threads);
                sporos_seed_index_free(index);
                for (size_t r = 0; r < n; ++r) {
                    size_t count;
                    fields >> count;
                    std::vector<SeedHit> expected(count);
                    for (SeedHit &hit : expected)
                        fields >> std::get<0>(hit) >> std::get<1>(hit) >> std::get<2>(hit) >> std::get<3>(hit) >>
                            std::get<4>(hit);
                    same_answers &= got[r] == expected;
                }
            }
            CHECK(same_answers);
            if (!same_answers) std::fprintf(stderr, "  %s\n", line.substr(0, 60).c_str());
            continue;
        }
        if (kind == "cutoff") {
            // The name field holds k.
            int k = std::stoi(name), w, scheme, open;
            double keep;
            int64_t expected;
            fields >> w >> scheme >> open >> keep >> expected;
            Many many(refs);
            sporos_seed_index *index = sporos_seed_index_new(many.seqs.data(), many.lens.data(), 2, k, w, scheme, open, 1);
            CHECK(index != nullptr && sporos_seed_index_occurrence_cutoff(index, keep) == expected);
            sporos_seed_index_free(index);
            continue;
        }
        int min_len, mode, reverse, max_occ, min_occ, threads;
        size_t count;
        fields >> min_len >> mode >> reverse >> max_occ >> min_occ >> threads >> count;
        if (!indexes.count(min_len)) indexes[min_len] = index_of(refs, min_len, threads);
        std::vector<Match> expected(count);
        for (Match &match : expected)
            fields >> std::get<0>(match) >> std::get<1>(match) >> std::get<2>(match) >> std::get<3>(match);
        std::vector<Match> got = matches(indexes[min_len], seq, mode, reverse, threads, max_occ, min_occ);
        CHECK(got == expected);
        if (got != expected) std::fprintf(stderr, "  %s\n", line.substr(0, 60).c_str());
    }
    CHECK(records > 0);

    // Four threads at once on one index, each answer the same as alone.
    sporos_index *index = indexes.begin()->second;
    const std::string &query = seqs["query"];
    std::vector<Match> alone = matches(index, query, SPOROS_MUM, 0, 1);
    std::vector<std::vector<Match>> each(4);
    std::vector<std::thread> threads;
    for (int t = 0; t < 4; ++t)
        threads.emplace_back([&, t] {
            for (int round = 0; round < 5; ++round) each[t] = matches(index, query, SPOROS_MUM, 0, 2);
        });
    for (std::thread &thread : threads) thread.join();
    for (const std::vector<Match> &found : each) CHECK(found == alone);
    for (auto &entry : indexes) sporos_index_free(entry.second);

    // A saved FM-index loads with its names and the same SMEMs.
    std::vector<std::string> reads;
    for (int r = 0; seqs.count("read" + std::to_string(r)); ++r) reads.push_back(seqs["read" + std::to_string(r)]);
    CHECK(fm != nullptr && !reads.empty());
    if (!fm) return;
    const char *saved = "build/test_c_api.sfm";
    CHECK(sporos_fm_index_save(fm, saved) == 0);
    sporos_fm_index *loaded = sporos_fm_index_load(saved, 2);
    CHECK(loaded != nullptr && sporos_fm_index_sequences(loaded) == 2);
    if (loaded) {
        CHECK(std::string(sporos_fm_index_name(loaded, 1)) == "ref1");
        CHECK(smems(loaded, reads, 19, 20, 1, 0, 3) == smems(fm, reads, 19, 20, 1, 0, 1));
        sporos_fm_index_free(loaded);
    }
    sporos_fm_index_free(fm);
}

// Runs every check, and fails when any did.
int main(int argc, char **argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: test_c_api <expected.txt from tests/api_expected.mojo>\n");
        return 2;
    }
    known();
    random_matches();
    against_mojo(argv[1]);
    if (failures) {
        std::fprintf(stderr, "test_c_api: %d checks failed\n", failures);
        return 1;
    }
    std::printf("test_c_api: all checks passed\n");
    return 0;
}
