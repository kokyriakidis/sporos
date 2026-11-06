"""The k-mer count histogram's peaks (`ha_analyze_count`, hist.cpp:74): the homozygous coverage peak, which sets
the high-occurrence cut-off and the candidate cap, and the heterozygous one."""


@fieldwise_init
struct Peaks(ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """`ha_analyze_count`'s answer: the homozygous peak (its return) and the heterozygous one (-1 for none)."""

    var hom: Int
    var het: Int


def analyze_count(cnt: List[Int], start_cnt: Int) -> Peaks:
    """hist.cpp:74 `ha_analyze_count` without the printing and without a given genome size (`m_peak_hom < 0`)."""
    var n_cnt = len(cnt)
    var start = 1 if cnt[1] > 0 else 2
    var low_i = start if start > start_cnt else start_cnt
    var i = low_i + 1
    while i < n_cnt:
        if cnt[i] > cnt[i - 1]:
            break
        i += 1
    low_i = i - 1
    if low_i == n_cnt - 1:
        return Peaks(-1, -1)
    var max_i = low_i + 1
    var max_v = cnt[max_i]
    for j in range(low_i + 1, n_cnt):
        if cnt[j] > max_v:
            max_v = cnt[j]
            max_i = j
    # A smaller peak on the low end.
    var max2: Int = -1
    var max2_i = -1
    var j = max_i - 1
    while j > low_i:
        if cnt[j] >= cnt[j - 1] and cnt[j] >= cnt[j + 1]:
            if cnt[j] > max2:
                max2 = cnt[j]
                max2_i = j
        j -= 1
    if max2_i > low_i and max2_i < max_i:
        var min_v = max_v
        for m in range(max2_i + 1, max_i):
            if cnt[m] < min_v:
                min_v = cnt[m]
        if Float64(max2) < Float64(max_v) * 0.05 or Float64(min_v) > Float64(max2) * 0.95:
            max2_i = -1
    # A smaller peak on the high end.
    var max3: Int = -1
    var max3_i = -1
    for m in range(max_i + 1, n_cnt - 1):
        if cnt[m] >= cnt[m - 1] and cnt[m] >= cnt[m + 1]:
            if cnt[m] > max3:
                max3 = cnt[m]
                max3_i = m
    if max3_i > max_i:
        var min_v = max_v
        for m in range(max_i + 1, max3_i):
            if cnt[m] < min_v:
                min_v = cnt[m]
        if (
            Float64(max3) < Float64(max_v) * 0.05
            or Float64(min_v) > Float64(max3) * 0.95
            or Float64(max3_i) > Float64(max_i) * 2.5
        ):
            max3_i = -1
    if max3_i > 0:
        return Peaks(max3_i, max_i)
    return Peaks(max_i, max2_i if max2_i > 0 else -1)
