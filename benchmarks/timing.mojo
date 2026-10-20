"""Timing helpers that hold up on a busy machine.

`thread_cpu_ns` reads the calling thread's CPU time, which excludes time the
thread spends descheduled, so single-thread measurements stay comparable
while other jobs run. Report the minimum of several runs.
"""

from std.ffi import external_call
from std.sys.info import CompilationTarget


def thread_cpu_ns() -> Int:
    # CLOCK_THREAD_CPUTIME_ID is 16 on macOS and 3 on Linux.
    comptime clock_id = 16 if CompilationTarget.is_macos() else 3
    var spec = List[Int64](length=2, fill=0)
    _ = external_call["clock_gettime", Int32](Int32(clock_id), spec.unsafe_ptr())
    return Int(spec[0]) * 1_000_000_000 + Int(spec[1])
