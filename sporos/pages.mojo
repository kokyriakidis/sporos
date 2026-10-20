"""Transparent huge pages for large buffers that are filled at random or by
many threads at once."""

from std.ffi import external_call
from std.sys.info import CompilationTarget


def advise_huge_pages(address: Int, bytes: Int):
    """Backs an uninitialized buffer with transparent huge pages on Linux
    (needed where THP is `madvise`-only); a no-op elsewhere. The buffer's
    whole pages are released first, since the allocator may hand out memory
    already faulted in as small pages: their contents are lost."""
    comptime if CompilationTarget.is_linux():
        comptime MADV_DONTNEED: Int32 = 4
        comptime MADV_HUGEPAGE: Int32 = 14
        var start = (address + 4095) & ~4095
        var end = (address + bytes) & ~4095
        if end - start >= 1 << 21:
            _ = external_call["madvise", Int32](start, end - start, MADV_HUGEPAGE)
            _ = external_call["madvise", Int32](start, end - start, MADV_DONTNEED)
