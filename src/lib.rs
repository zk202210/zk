#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(stdsimd)]
#![feature(stmt_expr_attributes)]

mod macros;
mod eventfd;
pub mod field;
pub mod merkle;
pub mod hash;
pub mod ldt;
pub mod vpd;
pub mod vpdv2;
pub mod mpcvpd;
pub mod plainvpd;
pub mod plainvpd_slow;
pub mod circuit;
pub mod tape;
pub mod poly;
pub mod util;

pub mod gkr;

pub mod zkp;
pub mod app;
pub mod matmul;

#[allow(non_snake_case)]
pub mod channel;

extern crate queues;

// #[cfg(debug_assertions)]
pub mod statistic;
