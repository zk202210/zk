// https://docs.rs/psx/latest/src/psx/macros.rs.html
//
// Gets a file size in bytes at compile-time.
#[macro_export]
macro_rules! file_size {
    ($file:expr) => {{
        const N: usize = include_bytes!($file).len();
        N
    }};
}

#[macro_export]
macro_rules! include_qwords {
    ($file:expr) => {{
        const N: usize = $crate::file_size!($file);
        const _: () = {
            if N % 8 != 0 {
                panic!("File size is not a multiple of 8 bytes");
            }
        };
        const M: usize = N / 8;
        const RET: [u64; M] = unsafe {
            std::mem::transmute(*include_bytes!($file))
        };
        &RET
    }};
}

// Includes the specified file as a `&mut [u32; N]`. This is exactly like
// `include_bytes`, expect that it checks the size is a multiple of 4 bytes at
// compile-time to allow creating a mutable `u32` slice rather than a `u8`
// slice.
#[macro_export]
macro_rules! include_dwords {
    ($file:expr) => {{
        const N: usize = $crate::file_size!($file);
        const _: () = {
            if N % 4 != 0 {
                panic!("File size is not a multiple of 4 bytes");
            }
        };
        const M: usize = N / 4;
        const RET: [u32; M] = unsafe {
            std::mem::transmute(*include_bytes!($file))
        };
        &RET
    }};
}

#[macro_export]
macro_rules! include_words {
    ($file:expr) => {{
        const N: usize = $crate::file_size!($file);
        const _: () = {
            if N % 2 != 0 {
                panic!("File size is not a multiple of 2 bytes");
            }
        };
        const M: usize = N / 2;
        const RET: [u16; M] = unsafe {
            std::mem::transmute(*include_bytes!($file))
        };
        &RET
    }};
}

pub(crate) use include_words;
pub(crate) use include_dwords;
#[allow(unused_imports)]
pub(crate) use include_qwords;
