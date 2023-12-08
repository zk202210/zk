#[link(name = "c")]
extern "C" {
    pub fn eventfd(initval: u32, flag: i32) -> i32;
    pub fn eventfd_read(fd: i32, value: *mut u64) -> i32;
    pub fn eventfd_write(fd: i32, value: u64) -> i32;
    pub fn close(fd: i32) -> i32;
}
