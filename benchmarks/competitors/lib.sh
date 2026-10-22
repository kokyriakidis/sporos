# Shared helpers for run_mems.sh and run_seeds.sh (sourced, not run).
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
BENCH=${BENCH:-$HOME/sporos-bench}
D=$BENCH/data
OUT=$BENCH/out
ENVBIN=$BENCH/env/.pixi/envs/default/bin
SBIN=${SBIN:-$BENCH/bin}  # sporos binaries snapshotted by setup_xaros.sh
RES=${RES:-$HERE/results}
REPS=${REPS:-3}
# A first run slower than this many seconds is not repeated (keeps the large pair affordable).
REPEAT_UNDER=${REPEAT_UNDER:-60}
mkdir -p "$OUT" "$RES"

# SPOROS_ONLY=1 measures sporos alone: every other tool's rows stay as the last full run left them, and the
# output checks run again against the outputs that run kept in $OUT.
SPOROS_ONLY=${SPOROS_ONLY:-0}
# True when other tools are not measured again: `cached || measure ...`.
cached() { [ "$SPOROS_ONLY" = 1 ]; }
# drop_rows <tsv> <column> <regex>: removes the rows whose column matches (the header stays), so they can be
# measured again.
drop_rows() {
  [ -f "$1" ] || return 0
  awk -F'\t' -v c="$2" -v re="$3" 'NR == 1 || $c !~ re' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}

# Busy CPUs over a short window, from /proc/stat (nothing of ours runs while it samples).
busy_cpus() {
  local a b
  a=($(head -1 /proc/stat)); sleep 0.5; b=($(head -1 /proc/stat))
  local tot=0 idle i
  for i in 1 2 3 4 5 6 7 8; do tot=$((tot + b[i] - a[i])); done
  idle=$(( (b[4] - a[4]) + (b[5] - a[5]) ))
  awk -v t="$tot" -v i="$idle" -v n="$(nproc)" 'BEGIN{printf "%.2f", (t>0 ? (t-i)/t*n : 0)}'
}

# Wait (up to ~10 min) while other jobs keep more than $BUSY_MAX CPUs busy.
BUSY_MAX=${BUSY_MAX:-1.5}
wait_idle() {
  local tries=0 b
  while :; do
    b=$(busy_cpus)
    awk -v b="$b" -v m="$BUSY_MAX" 'BEGIN{exit !(b <= m)}' && return 0
    tries=$((tries + 1))
    [ $tries -gt 60 ] && { echo "warning: machine still busy ($b CPUs), running anyway" >&2; return 0; }
    echo "machine busy ($b CPUs), waiting" >&2; sleep 5
  done
}

load1() { cut -d' ' -f1 /proc/loadavg; }

# measure <tsv> <label...tab-separated fields> -- <stdout file|-> cmd...
# Appends: fields, reps, best wall s, user s (of best run), sys s, max RSS MB, busy CPUs before, loadavg.
measure() {
  local tsv=$1; shift
  local fields=()
  while [ "$1" != "--" ]; do fields+=("$1"); shift; done
  shift
  local out=$1; shift
  local best=; local bu=; local bs=; local rss=0 r=0 tf
  tf=$(mktemp)
  wait_idle
  local busy; busy=$(busy_cpus)
  local ld; ld=$(load1)
  while [ $r -lt "$REPS" ]; do
    if [ "$out" = "-" ]; then
      /usr/bin/time -v "$@" > /dev/null 2> "$tf" || { echo "FAILED: $*" >&2; tail -5 "$tf" >&2; }
    else
      /usr/bin/time -v "$@" > "$out" 2> "$tf" || { echo "FAILED: $*" >&2; tail -5 "$tf" >&2; }
    fi
    r=$((r + 1))
    local wall user sys mrss
    wall=$(awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; print s}' "$tf")
    user=$(awk -F': ' '/User time \(seconds\)/{print $2}' "$tf")
    sys=$(awk -F': ' '/System time \(seconds\)/{print $2}' "$tf")
    mrss=$(awk -F': ' '/Maximum resident set size/{printf "%.0f", $2/1024}' "$tf")
    [ "$mrss" -gt "$rss" ] && rss=$mrss
    if [ -z "$best" ] || awk -v a="$wall" -v b="$best" 'BEGIN{exit !(a<b)}'; then best=$wall; bu=$user; bs=$sys; fi
    awk -v a="$wall" -v m="$REPEAT_UNDER" 'BEGIN{exit !(a>m)}' && break
  done
  rm -f "$tf"
  local line
  line=$(printf "%s\t" "${fields[@]}")
  printf "%s%s\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\n" "$line" "$r" "$best" "$bu" "$bs" "$rss" "$busy" "$ld" | tee -a "$tsv"
}

machine_info() {
  echo "host: $(hostname)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') ($(nproc) threads)"
  echo "mem: $(free -g | awk '/Mem:/{print $2}') GB"
  echo "kernel: $(uname -r)"
  echo "date: $(date -Iseconds)"
  echo "uptime: $(uptime)"
}
