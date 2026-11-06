//! rapidgzip-core behind a C interface, for sporos's `ByteStream` (sporos/stream.mojo): one gzip file (one member
//! or many, plain gzip or BGZF) decoded by many threads with rapidgzip's marker/window method, and read in order.
//! A thin layer: sporos does everything with the bytes.
//!
//! Conventions:
//! - A call that can fail returns 0 (or a count), or -1 with a message `srgz_error` reads (per thread, until its
//!   next failure).
//! - Paths go in as a pointer and a length, bytes as Unix takes them, no terminator.
//! - A stream is used by one thread at a time; its decoding runs on its own worker threads.

use rapidgzip_core::{Decoder, DecoderReader};
use std::cell::RefCell;
use std::io::Read;
use std::panic::{AssertUnwindSafe, catch_unwind};

/// The interface's version: sporos checks it, so a library left from another build is not used by mistake.
pub const ABI: u32 = 1;

#[unsafe(no_mangle)]
pub extern "C" fn srgz_abi() -> u32 {
    ABI
}

thread_local! {
    static ERROR: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

fn fail(message: String) -> i64 {
    ERROR.with(|e| *e.borrow_mut() = message.into_bytes());
    -1
}

/// The last failure's message on this thread (its length into `len`).
///
/// # Safety
/// `len` must be writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn srgz_error(len: *mut usize) -> *const u8 {
    ERROR.with(|e| {
        let e = e.borrow();
        // SAFETY: the caller passes a writable `len`.
        unsafe { *len = e.len() };
        e.as_ptr()
    })
}

/// An open file being decoded; `None` once its end was read and verified, or after a failure.
pub struct Stream {
    reader: Option<DecoderReader>,
}

fn guarded<F: FnOnce() -> Result<i64, String>>(body: F) -> i64 {
    match catch_unwind(AssertUnwindSafe(body)) {
        Ok(Ok(value)) => value,
        Ok(Err(message)) => fail(message),
        Err(_) => fail(String::from("rapidgzip panicked")),
    }
}

/// Opens `path` for decoding on up to `threads` worker threads (0: the machine's default); the stream into `out`.
///
/// # Safety
/// `path` must point to `len` readable bytes and `out` must be writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn srgz_open(path: *const u8, len: usize, threads: u32, out: *mut *mut Stream) -> i64 {
    guarded(|| {
        // SAFETY: the caller passes `len` readable bytes.
        let bytes = unsafe { std::slice::from_raw_parts(path, len) };
        let path = std::path::Path::new(<std::ffi::OsStr as std::os::unix::ffi::OsStrExt>::from_bytes(bytes));
        let decoder = Decoder::builder()
            .decoder_threads(threads as usize)
            .build()
            .map_err(|e| format!("rapidgzip configuration: {e}"))?;
        let reader = decoder.open(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let stream = Box::new(Stream { reader: Some(reader) });
        // SAFETY: the caller passes a writable `out`.
        unsafe { *out = Box::into_raw(stream) };
        Ok(0)
    })
}

/// Fills `buf` (`len` bytes) with the next decoded bytes; how many, fewer only at the end, 0 once the end was
/// reached and every member's checksum and size verified.
///
/// # Safety
/// `stream` must come from `srgz_open` and not be closed; `buf` must be `len` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn srgz_read(stream: *mut Stream, buf: *mut u8, len: usize) -> i64 {
    guarded(|| {
        // SAFETY: the caller passes an open stream and `len` writable bytes.
        let stream = unsafe { &mut *stream };
        let out = unsafe { std::slice::from_raw_parts_mut(buf, len) };
        let Some(reader) = stream.reader.as_mut() else {
            return Ok(0);
        };
        let mut done = 0;
        while done < len {
            let got = reader.read(&mut out[done..]).map_err(|e| format!("rapidgzip: {e}"))?;
            if got == 0 {
                // The end: verify what remains (every member's CRC32 and size) before saying so.
                let reader = stream.reader.take().expect("the reader is present");
                reader.finish().map_err(|e| format!("rapidgzip: {e}"))?;
                break;
            }
            done += got;
        }
        Ok(done as i64)
    })
}

/// Closes a stream, cancelling any decoding still running.
///
/// # Safety
/// `stream` must come from `srgz_open` and not be closed already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn srgz_close(stream: *mut Stream) {
    if !stream.is_null() {
        // SAFETY: the caller passes a stream from `srgz_open`, closed only here.
        drop(unsafe { Box::from_raw(stream) });
    }
}
