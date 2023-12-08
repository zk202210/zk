#!/usr/bin/env python3

from pyfinite import ffield


# fn euclidean_div_quotient(a: u64, b: u64) -> u64 {
#     let mut q: u64 = 0;
#     let mut a: u64 = a;
#     let mut diff: i32 = b.leading_zeros() as i32 - a.leading_zeros() as i32;
#     while diff >= 0 && a != 0 {
#         q |= 1 << diff;
#         a ^= b << diff;
#         diff = b.leading_zeros() as i32 - a.leading_zeros() as i32;
#     }
#     q
# }
# pub fn init() {
#     unsafe {
#         // init inv table
#         // modolus: x^16 + x^5 + x^3 + x + 1
#         let modulus: u64 = 0b10000000000101011;
#         for idx in 0..65536 {
#             let mut t0 = 0;
#             let mut t1 = 1;
#             let mut r0 = modulus;
#             let mut r1 = idx;
#             while r1 != 0 {
#                 let q = euclidean_div_quotient(r0, r1);
#                 (r0, r1) = (r1, r0 ^ clmul(q, r1));
#                 (t0, t1) = (t1, t0 ^ clmul(q, t1));
#             }
#             GF2P16_INV_LUT[idx as usize] = t0 as u16;
#         }
#         // init lift table
#         // ring morphism: x |--> y^14 + y^8 + y^7 + y^6 + y^3 + y^2 + 1
#         GF2P16_LFT_LUT[0] = 0;
#         GF2P16_LFT_LUT[1] = 1;
#         let gen = Self(0b100000111001101);
#         let mut pow = gen;
#         for bit in 1..8 {
#             let s = 1 << bit;
#             for idx in 0..s {
#                 GF2P16_LFT_LUT[s + idx] = GF2P16_LFT_LUT[idx] ^ pow.0;
#             }
#             pow *= gen;
#         }
#     }
# }


def gf2p8_mut_lut_str():
    F = ffield.FField(8, gen=0b100011011, useLUT=0)
    with open("../resources/mut.lut", "w") as f:
        for i in range(256):
            f.write("    ")
            for j in range(256):
                f.write("0x{:02x},".format(F.Multiply(i, j)))
                f.write("\n")


def gf2p8_mut_lut_bin():
    F = ffield.FField(8, gen=0b100011011, useLUT=0)
    with open("../resources/gf2p8.mut.lut", "wb") as f:
        for i in range(256):
            for j in range(256):
                f.write(F.Multiply(i, j).to_bytes(1, 'little'))


def gf2p8_inv_lut_bin():
    F = ffield.FField(8, gen=0b100011011, useLUT=0)
    with open("../resources/gf2p8.inv.lut", "wb") as f:
        f.write(b'\x00')
        for i in range(1, 256):
            f.write(F.Inverse(i).to_bytes(1, 'little'))


def packgf2p16_mut_lut_bin():
    F = ffield.FField(16, gen=0b10000000000101011, useLUT=0)
    with open("../resources/packgf2p16.mut.lut", "wb") as f:
        for i in range(0, 65536):
            # a0lo, a2hi, a2lo, a0hi, a1lo, a3hi, a3lo, a1hi
            for a in range(0, 16):
                a0 = F.Multiply(i, a)
                f.write((a0 & 0xff).to_bytes(1, 'little'))
            for a in range(0, 16):
                a2 = F.Multiply(i, a << 8)
                f.write((a2 >> 8).to_bytes(1, 'little'))

            for a in range(0, 16):
                a2 = F.Multiply(i, a << 8)
                f.write((a2 & 0xff).to_bytes(1, 'little'))
            for a in range(0, 16):
                a0 = F.Multiply(i, a)
                f.write((a0 >> 8).to_bytes(1, 'little'))

            for a in range(0, 16):
                a1 = F.Multiply(i, a << 4)
                f.write((a1 & 0xff).to_bytes(1, 'little'))
            for a in range(0, 16):
                a3 = F.Multiply(i, a << 12)
                f.write((a3 >> 8).to_bytes(1, 'little'))

            for a in range(0, 16):
                a3 = F.Multiply(i, a << 12)
                f.write((a3 & 0xff).to_bytes(1, 'little'))
            for a in range(0, 16):
                a1 = F.Multiply(i, a << 4)
                f.write((a1 >> 8).to_bytes(1, 'little'))

            # # a0lo, a2hi, a1lo, a3hi, a2lo, a0hi, a3lo, a1hi
            # for a in range(0, 16):
            #     a0 = F.Multiply(i, a)
            #     f.write((a0 & 0xff).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a2 = F.Multiply(i, a << 8)
            #     f.write((a2 >> 8).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a1 = F.Multiply(i, a << 4)
            #     f.write((a1 & 0xff).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a3 = F.Multiply(i, a << 12)
            #     f.write((a3 >> 8).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a2 = F.Multiply(i, a << 8)
            #     f.write((a2 & 0xff).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a0 = F.Multiply(i, a)
            #     f.write((a0 >> 8).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a3 = F.Multiply(i, a << 12)
            #     f.write((a3 & 0xff).to_bytes(1, 'little'))
            # for a in range(0, 16):
            #     a1 = F.Multiply(i, a << 4)
            #     f.write((a1 >> 8).to_bytes(1, 'little'))


def gf2p16_mut_lut_bin():
    F = ffield.FField(16, gen=0b10000000000101011, useLUT=0)
    with open("../resources/gf2p16.mut.lut", "wb") as f:
        for i in range(0, 256):
            for a in range(0, 256):
                a0 = F.Multiply(i, a)
                f.write(a0.to_bytes(2, 'little'))
            for a in range(0, 256):
                a1 = F.Multiply(i, a << 8)
                f.write(a1.to_bytes(2, 'little'))
        for i in range(0, 256):
            for a in range(0, 256):
                a0 = F.Multiply(i << 8, a)
                f.write(a0.to_bytes(2, 'little'))
            for a in range(0, 256):
                a1 = F.Multiply(i << 8, a << 8)
                f.write(a1.to_bytes(2, 'little'))



def gf2p16_inv_lut_bin():
    F = ffield.FField(16, gen=0b10000000000101011, useLUT=0)
    with open("../resources/gf2p16.inv.lut", "wb") as f:
        f.write(b'\x00\x00')
        for i in range(1, 65536):
            f.write(F.Inverse(i).to_bytes(2, 'little'))


def gf2p16_lift_lut_bin():
    F = ffield.FField(16, gen=0b10000000000101011, useLUT=0)
    gen = 0b100000111001101
    power = gen
    LUT = [0] * 256
    with open("../resources/gf2p16.lft.lut", "wb") as f:
        f.write(b'\x00\x00')
        f.write(b'\x01\x00')
        LUT[0] = 0
        LUT[1] = 1
        for bit in range(1, 8):
            s = 1 << bit
            for idx in range(0, s):
                LUT[s + idx] = LUT[idx] ^ power
                f.write(LUT[s + idx].to_bytes(2, 'little'))
            power = F.Multiply(power, gen)


def gf2p32_mut_lut_bin():
    F = ffield.FField(32, gen=0b100000000000000000000000010001101, useLUT=0)
    with open("../resources/gf2p32.mut.lut", "wb") as f:
        for k in range(0, 2):
            for i in range(0, 65536):
                j = i << (k * 16)
                # lo
                # a0lo, a2hi, a1lo, a3hi, a2lo, a0hi, a3lo, a1hi
                for a in range(0, 16):
                    a0 = F.Multiply(j, a) & 0xffff
                    f.write((a0 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 8) & 0xffff
                    f.write((a2 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 4) & 0xffff
                    f.write((a1 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 12) & 0xffff
                    f.write((a3 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 8) & 0xffff
                    f.write((a2 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a0 = F.Multiply(j, a) & 0xffff
                    f.write((a0 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 12) & 0xffff
                    f.write((a3 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 4) & 0xffff
                    f.write((a1 >> 8).to_bytes(1, 'little'))
                # hi
                for a in range(0, 16):
                    a0 = F.Multiply(j, a) >> 16
                    f.write((a0 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 8) >> 16
                    f.write((a2 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 4) >> 16
                    f.write((a1 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 12) >> 16
                    f.write((a3 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 8) >> 16
                    f.write((a2 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a0 = F.Multiply(j, a) >> 16
                    f.write((a0 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 12) >> 16
                    f.write((a3 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 4) >> 16
                    f.write((a1 >> 8).to_bytes(1, 'little'))

                # a4lo, a6hi, a5lo, a7hi, a6lo, a4hi, a7lo, a5hi
                # lo
                for a in range(0, 16):
                    a0 = F.Multiply(j, a << 16) & 0xffff
                    f.write((a0 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 24) & 0xffff
                    f.write((a2 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 20) & 0xffff
                    f.write((a1 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 28) & 0xffff
                    f.write((a3 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 24) & 0xffff
                    f.write((a2 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a0 = F.Multiply(j, a << 16) & 0xffff
                    f.write((a0 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 28) & 0xffff
                    f.write((a3 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 20) & 0xffff
                    f.write((a1 >> 8).to_bytes(1, 'little'))
                # hi
                for a in range(0, 16):
                    a0 = F.Multiply(j, a << 16) >> 16
                    f.write((a0 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 24) >> 16
                    f.write((a2 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 20) >> 16
                    f.write((a1 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 28) >> 16
                    f.write((a3 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a2 = F.Multiply(j, a << 24) >> 16
                    f.write((a2 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a0 = F.Multiply(j, a << 16) >> 16
                    f.write((a0 >> 8).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a3 = F.Multiply(j, a << 28) >> 16
                    f.write((a3 & 0xff).to_bytes(1, 'little'))
                for a in range(0, 16):
                    a1 = F.Multiply(j, a << 20) >> 16
                    f.write((a1 >> 8).to_bytes(1, 'little'))


# def gf2p32_inv_lut_bin():
#     F = ffield.FField(32, gen=0b100000000000000000000000010001101, useLUT=0)
#     with open("../resources/gf2p32.inv.lut", "wb") as f:
#         f.write(b'\x00\x00')
#         for i in range(1, 2**32):
#             f.write(F.Inverse(i).to_bytes(4, 'little'))


def gf2p32_lift_lut_bin():
    F = ffield.FField(32, gen=0b100000000000000000000000010001101, useLUT=0)
    # TODO: compute gen
    gen = 0b100000111001101
    power = gen
    LUT = [0] * 256
    with open("../resources/gf2p32.lft.lut", "wb") as f:
        f.write(b'\x00\x00\x00\x00')
        f.write(b'\x01\x00\x00\x00')
        LUT[0] = 0
        LUT[1] = 1
        for bit in range(1, 8):
            s = 1 << bit
            for idx in range(0, s):
                LUT[s + idx] = LUT[idx] ^ power
                f.write(LUT[s + idx].to_bytes(4, 'little'))
            power = F.Multiply(power, gen)


if __name__ == '__main__':
    # gf2p8_mut_lut_bin()
    # gf2p8_inv_lut_bin()
    # gf2p16_inv_lut_bin()
    # gf2p16_lift_lut_bin()
    # gf2p16_mut_lut_bin()
    packgf2p16_mut_lut_bin()
    # gf2p32_lift_lut_bin()
    # gf2p32_mut_lut_bin()
