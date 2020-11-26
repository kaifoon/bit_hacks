pub struct BitTwiddlingHacks {}

pub const INT_BITS: i32 = 32; //number of bits in an integer
pub const INT_MAX: i32 = 0x7fffffff;
pub const INT_MIN: i32 = -1 << (INT_BITS - 1);

//bool to i32
#[inline]
pub fn bi32(b: bool) -> i32 {
    b as i32
}

///bool to u32
#[inline]
pub fn bu32(b: bool) -> u32 {
    b as u32
}

///u32 to usize
#[inline]
pub fn usz(v: u32) -> usize {
    v as usize
}

impl BitTwiddlingHacks {
    // NOT
    pub fn not(v: i32) -> i32 {
        -(v + 1)
    }
    // 计算整数的符号，返回 -1, 0, 1
    // 负数 => -1
    // 0 => 0
    // 正数 => 1
    pub fn sign(v: i32) -> i32 {
        bi32(v > 0) - bi32(v < 0)
    }

    // 检测两个整数是否具有相反的符号
    // 相反返回true
    pub fn diff_sign(x: i32, y: i32) -> bool {
        (x ^ y) < 0
    }

    // 无分支的求绝对值
    pub fn abs(v: i32) -> i32 {
        let mask = v >> (INT_BITS - 1);

        (v ^ mask) - mask
    }

    // 无分支求最小最大值
    pub fn min(x: i32, y: i32) -> i32 {
        y ^ ((x ^ y) & -(bi32(x < y)))
    }

    pub fn max(x: i32, y: i32) -> i32 {
        x ^ ((x ^ y) & -(bi32(x < y)))
    }

    // 确定整数是否为 2 的幂
    pub fn is_pow2(v: i32) -> bool {
        v != 0 && ((v & (v - 1)) == 0)
    }

    // 条件设定位，没有分支
    // 对超标量CPU
    // ! is rust bitwise not
    pub fn condi_set(f: bool, mask: i32, w: i32) -> i32 {
        (w & !mask) | (-bi32(f) & mask)
    }

    // 条件设定负数，没有分支
    // f = true, 设为负数
    pub fn condi_neg(f: bool, v: i32) -> i32 {
        (v ^ -bi32(f)) + bi32(f)
    }

    // 根据 mask 来合并位
    pub fn merge_bits(a: i32, b: i32, mask: i32) -> i32 {
        a ^ ((a ^ b) & mask)
    }

    // 计算位为1的个数
    // Brian Kernighan's way
    pub fn count_bits(v: i32) -> i32 {
        let mut c = 0;
        let mut v = v as u32;
        while v != 0 {
            v &= v - 1;
            c += 1;
        }
        c
    }

    // 从 MSB到k位计算位为1的个数
    // 使用并行计算
    // 10101110，k=4，得到1010，所以个数为2
    pub fn count_bits_k(v: i32, k: i32) -> i32 {
        let mut r = v as u32 >> (INT_BITS - k);
        r = r - ((r >> 1) & !0u32 / 3);
        r = (r & !0u32 / 5) + ((r >> 2) & !0u32 / 5);
        r = (r + (r >> 4)) & !0u32 / 17;
        r = r % 255;

        r as i32
    }

    // 计算几偶性
    pub fn parity_native(v: i32) -> bool {
        let mut v = v as u32;
        let mut parity = false;
        while v != 0 {
            parity = !parity;
            v &= v - 1;
        }

        parity
    }

    pub fn parity_parallel(v: i32) -> bool {
        let mut v = v as u32;
        v ^= v >> 16;
        v ^= v >> 8;
        v ^= v >> 4;
        v &= 0xf;

        ((0x6996 >> v) & 0x1) == 1
    }

    // 交换值
    pub fn swap(mut x: i32, mut y: i32) -> (i32, i32) {
        x ^= y;
        y ^= x;
        x ^= y;

        (x, y)
    }

    // 交换特定的位
    // 如 交换从i = 1 开始的3个位 到从j = 5开始的3个位
    // 00101111 -> 11100011
    pub fn swap_bits(v: i32, i: i32, j: i32, n: i32) -> i32 {
        let v = v as u32;
        let x = ((v >> i) ^ (v >> j)) & ((1 << n) - 1);

        (v ^ ((x << i) | (x << j))) as i32
    }

    // 逆转位
    // 0b10110010  -> 0b01001101
    pub fn reverse_bits_native(v: i32) -> i32 {
        let mut s = INT_BITS - 1;
        let mut v = v as u32;
        let mut r = v;

        v >>= 1;
        while v != 0 {
            r <<= 1;
            r |= v & 1;
            s -= 1;
            v >>= 1;
        }

        r <<= s;

        r as i32
    }

    // 适用大数的位逆转 N-bit log(N)条指令
    pub fn reverse_bits_logn(mut v: i128) -> i128 {
        use std::mem;
        // size_of_val == 16
        let mut s = mem::size_of_val(&v) * 8;
        let mut mask = !0;

        s >>= 1;
        while s > 0 {
            mask ^= mask << s;
            v = ((v >> s) & mask) | ((v << s) & !mask);
            s >>= 1;
        }

        v
    }

    // 模除法
    // 当 n % d , n是正整数， d 是2 的幂适用
    pub fn modiv(n: i32, d: i32) -> i32 {
        n & (d - 1)
    }

    // log2 求 n位数的的log2 对数值， 时间复杂度O(logN)
    pub fn log2(v: u32) -> u32 {
        let mut v: u32 = v;
        let mut r: u32;
        let mut shift: u32;
        r = bu32(v > 0xFFFF) << 4;
        v >>= r;
        shift = bu32(v > 0xFF) << 3;
        v >>= shift;
        r |= shift;
        shift = bu32(v > 0xF) << 2;
        v >>= shift;
        r |= shift;
        shift = bu32(v > 0x3) << 1;
        v >>= shift;
        r |= shift;
        r | (v >> 1)
    }

    // log10对数值
    pub fn log10(v: i32) -> i32 {
        if v >= 1000000000 {
            9
        } else if v >= 100000000 {
            8
        } else if v >= 10000000 {
            7
        } else if v >= 1000000 {
            6
        } else if v >= 100000 {
            5
        } else if v >= 10000 {
            4
        } else if v >= 1000 {
            3
        } else if v >= 100 {
            2
        } else if v >= 10 {
            1
        } else {
            0
        }
    }

    // 从右边开始连续计算0的个数
    pub fn count_trailing_zeros(v: u32) -> u32 {
        let mut v = v;
        let mut c = 32u32;
        // NOTE: if 0 == v, then c = 31.
        //short circuit v == 0 as 32
        if v == 0 {
            return c;
        }
        if v & 0x1 != 0 {
            c = 0; // special case for odd v (assumed to happen half of the time)
        } else {
            c = 1;
            if (v & 0xffff) == 0 {
                v >>= 16;
                c += 16;
            }
            if (v & 0xff) == 0 {
                v >>= 8;
                c += 8;
            }
            if (v & 0xf) == 0 {
                v >>= 4;
                c += 4;
            }
            if (v & 0x3) == 0 {
                v >>= 2;
                c += 2;
            }
            c -= v & 0x1;
        }
        c
    }

    // 向上取pow2的值
    pub fn next_pow2(v: i32) -> i32 {
        let mut v = v;
        v += bi32(v == 0);
        v -= 1;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v + 1
    }

    // 向下取pow2的值
    pub fn prev_pow2(v: i32) -> i32 {
        let mut v = v;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v - (v >> 1)
    }

    // 2个坐标的交错位
    pub fn interleave2(x: u32, y: u32) -> u32 {
        let mut x = x;
        let mut y = y;
        x &= 0xFFFF;
        x = (x | (x << 8)) & 0x00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F;
        x = (x | (x << 2)) & 0x33333333;
        x = (x | (x << 1)) & 0x55555555;

        y &= 0xFFFF;
        y = (y | (y << 8)) & 0x00FF00FF;
        y = (y | (y << 4)) & 0x0F0F0F0F;
        y = (y | (y << 2)) & 0x33333333;
        y = (y | (y << 1)) & 0x55555555;

        x | (y << 1)
    }

    // 提取交错位的元素
    pub fn deinterleave2(v: u32, n: u32) -> u32 {
        let mut v = v;
        v = (v >> n) & 0x55555555;
        v = (v | (v >> 1)) & 0x33333333;
        v = (v | (v >> 2)) & 0x0F0F0F0F;
        v = (v | (v >> 4)) & 0x00FF00FF;
        v = (v | (v >> 16)) & 0x000FFFF;
        (v << 16) >> 16
    }

    // 3个坐标的交错位
    pub fn interleave3(x: u32, y: u32, z: u32) -> u32 {
        let mut x = x;
        let mut y = y;
        let mut z = z;
        x &= 0x3FF;
        x = (x | (x << 16)) & 4278190335;
        x = (x | (x << 8)) & 251719695;
        x = (x | (x << 4)) & 3272356035;
        x = (x | (x << 2)) & 1227133513;

        y &= 0x3FF;
        y = (y | (y << 16)) & 4278190335;
        y = (y | (y << 8)) & 251719695;
        y = (y | (y << 4)) & 3272356035;
        y = (y | (y << 2)) & 1227133513;
        x |= y << 1;

        z &= 0x3FF;
        z = (z | (z << 16)) & 4278190335;
        z = (z | (z << 8)) & 251719695;
        z = (z | (z << 4)) & 3272356035;
        z = (z | (z << 2)) & 1227133513;

        x | (z << 2)
    }

    // 提取交错位的元素
    pub fn deinterleave3(v: u32, n: u32) -> u32 {
        let mut v = v;
        v = (v >> n) & 1227133513;
        v = (v | (v >> 2)) & 3272356035;
        v = (v | (v >> 4)) & 251719695;
        v = (v | (v >> 8)) & 4278190335;
        v = (v | (v >> 16)) & 0x3FF;
        (v << 22) >> 22
    }

    // 按字典顺序计算下一位排列
    pub fn next_combination(v: u32) -> u32 {
        let t = v | (v - 1);
        let c = (Self::not(t as i32) & -Self::not(t as i32)) as u32 - 1;
        (t + 1) | (c >> (Self::count_trailing_zeros(v) + 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_not() {
        assert_eq!(BitTwiddlingHacks::not(170), -171);
        assert_eq!(BitTwiddlingHacks::not(0), -1);
        assert_eq!(BitTwiddlingHacks::not(-3), 2);
    }

    #[test]
    fn test_sign() {
        assert_eq!(BitTwiddlingHacks::sign(-100), -1);
        assert_eq!(BitTwiddlingHacks::sign(100), 1);
        assert_eq!(BitTwiddlingHacks::sign(0), 0);
        assert_eq!(BitTwiddlingHacks::sign(INT_MAX), 1);
        assert_eq!(BitTwiddlingHacks::sign(INT_MIN), -1);
    }

    #[test]
    fn test_abs() {
        assert_eq!(BitTwiddlingHacks::abs(0), 0);
        assert_eq!(BitTwiddlingHacks::abs(1), 1);
        assert_eq!(BitTwiddlingHacks::abs(-1), 1);
        assert_eq!(BitTwiddlingHacks::abs(INT_MAX), INT_MAX);
        assert_eq!(BitTwiddlingHacks::abs(-INT_MAX), INT_MAX);
        //BitTwiddlingHacks::abs(-INT_MIN) -- overflow
    }

    #[test]
    fn test_min() {
        assert_eq!(BitTwiddlingHacks::min(0, 0), 0);
        assert_eq!(BitTwiddlingHacks::min(-1, 1), -1);
        assert_eq!(BitTwiddlingHacks::min(INT_MAX, INT_MAX), INT_MAX);
        assert_eq!(BitTwiddlingHacks::min(INT_MIN, INT_MIN), INT_MIN);
        assert_eq!(BitTwiddlingHacks::min(INT_MAX, INT_MIN), INT_MIN);
    }

    #[test]
    fn test_max() {
        assert_eq!(BitTwiddlingHacks::max(0, 0), 0);
        assert_eq!(BitTwiddlingHacks::max(-1, 1), 1);
        assert_eq!(BitTwiddlingHacks::max(INT_MAX, INT_MAX), INT_MAX);
        assert_eq!(BitTwiddlingHacks::max(INT_MIN, INT_MIN), INT_MIN);
        assert_eq!(BitTwiddlingHacks::max(INT_MAX, INT_MIN), INT_MAX);
    }

    #[test]
    fn test_is_pow2() {
        assert!(!BitTwiddlingHacks::is_pow2(0));
        for i in 0..31 {
            assert!(BitTwiddlingHacks::is_pow2(1 << i));
        }
        assert!(!BitTwiddlingHacks::is_pow2(100));
        assert!(!BitTwiddlingHacks::is_pow2(0x7fffffff));
        assert!(!BitTwiddlingHacks::is_pow2(-1000000));
    }
    #[test]
    fn test_count_bits() {
        assert_eq!(BitTwiddlingHacks::count_bits(-10), 30);
    }
    #[test]
    fn test_count_bits_k() {
        assert_eq!(BitTwiddlingHacks::count_bits_k(-10, 2), 2);
        assert_eq!(BitTwiddlingHacks::count_bits_k(INT_MAX, 2), 1);
    }

    #[test]
    fn test_parity() {
        assert_eq!(BitTwiddlingHacks::parity_native(10), false);
        assert_eq!(BitTwiddlingHacks::parity_native(-10), false);
        assert_eq!(BitTwiddlingHacks::parity_native(8), true);
        assert_eq!(BitTwiddlingHacks::parity_parallel(8), true);
        assert_eq!(BitTwiddlingHacks::parity_parallel(-10), false);
        assert_eq!(BitTwiddlingHacks::parity_parallel(10), false);
    }

    #[test]
    fn test_swap_bits() {
        assert_eq!(
            BitTwiddlingHacks::swap_bits(0b00101111, 1, 5, 3),
            0b11100011
        );
    }

    #[test]
    fn test_reverse_bits() {
        assert_eq!(
            BitTwiddlingHacks::reverse_bits_native(0b00101111),
            -201326592
        );
        // 字面数(0x, 0b,0o)没有sign表示，需要转成unsigned， 然后转为 sign
        assert_eq!(
            BitTwiddlingHacks::reverse_bits_logn(0x7f123456789abcdedfedaba987d54321u128 as i128),
            0x84c2abe195d5b7fb7b3d591e6a2c48feu128 as i128
        );
    }

    #[test]
    fn test_count_trailing_zeros() {
        assert_eq!(BitTwiddlingHacks::count_trailing_zeros(0), 32);
        assert_eq!(BitTwiddlingHacks::count_trailing_zeros(1), 0);
        for i in 0..31 {
            assert_eq!(BitTwiddlingHacks::count_trailing_zeros(1 << i), i);
            if i > 0 {
                assert_eq!(BitTwiddlingHacks::count_trailing_zeros((1 << i) - 1), 0)
            }
        }
        assert_eq!(BitTwiddlingHacks::count_trailing_zeros(0xf81700), 8);
    }

    #[test]
    fn test_next_pow2() {
        for i in 0..31 {
            if i != 1 {
                assert_eq!(BitTwiddlingHacks::next_pow2((1 << i) - 1), 1 << i);
            }
            assert_eq!(BitTwiddlingHacks::next_pow2(1 << i), 1 << i);
            if i < 30 {
                assert_eq!(BitTwiddlingHacks::next_pow2((1 << i) + 1), 1 << (i + 1));
            }
        }
    }

    #[test]
    fn test_prev_pow2() {
        println!(
            "{i:>2}    {input:>w$}    {prev:>w$}",
            i = "i",
            input = "((1 << i) + 1)",
            prev = "prev_pow2",
            w = 10
        );
        println!("{}", "-".repeat(34));
        for i in 0..31 {
            if i > 0 {
                assert_eq!(BitTwiddlingHacks::prev_pow2((1 << i) - 1), 1 << (i - 1));
            }
            assert_eq!(BitTwiddlingHacks::prev_pow2(1 << i), 1 << i);

            if 0 < i && i < 30 {
                println!(
                    "{i:>2} .. {input:>w$} .. {prev:>w$}",
                    i = i,
                    input = ((1 << i) + 1),
                    prev = BitTwiddlingHacks::prev_pow2((1 << i) + 1),
                    w = 10
                );
                assert_eq!(BitTwiddlingHacks::prev_pow2((1 << i) + 1), 1 << i);
            }
        }
    }

    #[test]
    fn test_next_combination() {
        assert_eq!(BitTwiddlingHacks::next_combination(1), 2);
        assert_eq!(BitTwiddlingHacks::next_combination(0x300), 0x401);
    }

    #[test]
    fn test_interleave2() {
        for x in 0..100 {
            for y in 0..100 {
                let h = BitTwiddlingHacks::interleave2(x, y);
                assert_eq!(BitTwiddlingHacks::deinterleave2(h, 0), x);
                assert_eq!(BitTwiddlingHacks::deinterleave2(h, 1), y);
            }
        }
    }

    #[test]
    fn test_interleave3() {
        for x in 0..(25 + 1) {
            for y in 0..(25 + 1) {
                for z in 0..(25 + 1) {
                    let h = BitTwiddlingHacks::interleave3(x, y, z);
                    assert_eq!(BitTwiddlingHacks::deinterleave3(h, 0), x);
                    assert_eq!(BitTwiddlingHacks::deinterleave3(h, 1), y);
                    assert_eq!(BitTwiddlingHacks::deinterleave3(h, 2), z);
                }
            }
        }
    }
}
