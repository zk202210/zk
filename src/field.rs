#[cfg(all(target_arch = "x86", target_feature = "avx2"))]
use std::arch::x86::*;
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
use std::arch::x86_64::*;
use std::ops::*;
use std::fmt::Debug;
use rand::Rng;

use primitive_types::U256 as u256;
// use primitive_types::U512 as u512;

use super::macros::include_words;
use super::macros::include_dwords;

#[cfg(debug_assertions)]
pub use super::statistic::*;


pub mod common;
pub use self::common::Field;
pub use self::common::FieldPack;
pub use self::common::Pack;
pub use self::common::PField;
pub use self::common::VField;
pub use self::common::Size;
pub use self::common::Random;
pub use self::common::Pow;
pub use self::common::Characteristic;

mod coset;
pub use self::coset::Coset;
pub use self::coset::*;

mod gf2p8;
pub use self::gf2p8::GF2p8;

mod gf2p16;
pub use self::gf2p16::{GF2p16, P4GF2p16, P8GF2p16, P16GF2p16, VGF2p16};

mod gf2p32;
pub use self::gf2p32::{GF2p32, P4GF2p32, P8GF2p32, P16GF2p32, VGF2p32};

mod gf2p64;
pub use self::gf2p64::{GF2p64, P4GF2p64, P8GF2p64, P16GF2p64, VGF2p64};

mod gf2p128;
pub use self::gf2p128::{GF2p128, P4GF2p128, P8GF2p128, P16GF2p128, VGF2p128};

mod gf2p192;
pub use self::gf2p192::{GF2p192, P4GF2p192, P8GF2p192, P16GF2p192, VGF2p192};

mod gf2p256;
pub use self::gf2p256::{GF2p256, P4GF2p256, P8GF2p256, P16GF2p256, VGF2p256};

mod m61p2;
pub use self::m61p2::{M61p2, P4M61p2, P8M61p2, P16M61p2, VM61p2};

mod m127p2;
pub use self::m127p2::{M127p2, P4M127p2, P8M127p2, P16M127p2, VM127p2};

mod p192;
pub use self::p192::{P192, P4P192, P8P192, P16P192, VP192};

mod p256;
pub use self::p256::{P256, P4P256, P8P256, P16P256, VP256};

mod vgf2p16;
mod vgf2p32;
mod p4gf2p16;
mod p8gf2p16;
mod p16gf2p16;
mod p4m61p2;


#[macro_export]
macro_rules! field_env {
    ($fn:tt ::< $F:tt $(, $Ts:tt)* >($($arg:expr),*)) => {
        $fn::<
            $F,
            <$F as FieldPack>::P4,
            <$F as FieldPack>::P8,
            <$F as FieldPack>::P16
            $(, $Ts)*
            >($($arg),*)
    };
    // (fn $name:ident
    //  <$F:tt $(, $Ts:tt $(: $Traits:tt)*)*>
    //  ($($arg:ident : $type:ty),*)
    //  (-> $ret:ty)?
    //  (where ())?
    //  $body) => {
    // }
}

// pub(crate) use field_env;
