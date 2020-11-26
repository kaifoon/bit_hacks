pub struct BitTwiddlingHacks {}

///not ~
#[inline]
pub fn not(v: i32) -> i32 {
    -(v + 1)
}

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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign() {
        assert_eq!(BitTwiddlingHacks::sign(-10), -1);
        assert_eq!(BitTwiddlingHacks::sign(101), 1);
        assert_eq!(BitTwiddlingHacks::sign(0), 0);
    }
    #[test]
    fn test_diff_sign() {
        assert_eq!(BitTwiddlingHacks::diff_sign(-10, 10), true);
        assert_eq!(BitTwiddlingHacks::diff_sign(-10, -10), false);
    }
}
