use group::ff::PrimeField;
// use group::ff::Field;
use std::{
    fmt::Debug,
    fmt::Display,
    marker::PhantomData,
    ops::{Add, Index, Mul, Sub},
};
use std::sync::Mutex;

// ------------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct FF {
    val: i64,
    modulo : i64,
}

// Define a static variable with Mutex
lazy_static::lazy_static! {
    static ref SHARED_VALUE: Mutex<i64> = Mutex::new(1);
}

impl FF {
    fn get_shared_mod_value() -> i64 {
        *SHARED_VALUE.lock().unwrap()
    }

    fn set_shared_mod_value(new_value: i64) {
        *SHARED_VALUE.lock().unwrap() = new_value;
    }

    pub fn new(val: i64) -> Self {
        FF { val : val , modulo : Self::get_shared_mod_value() }
    }
    pub fn as_raw(&self) -> i64 {
        self.val
    }

    // Returns x^y mod m.
    pub fn pow_mod(&self, mut y: i64) -> Self {
        let mut x = self.val;
        if y < 0 || self.modulo <= 0 {
            panic!("Invalid input parameters");
        }
        if !(0 <= x && x < self.modulo) {
            x = ((x % self.modulo) + self.modulo) % self.modulo;
        }
        let mut result = 1;
        while y != 0 {
            if y & 1 != 0 {
                result = (result * x) % self.modulo;
            }
            x = (x * x) % self.modulo;
            y >>= 1;
        }
        FF::new(result)
    }

    // Returns x^-1 mod m.
    pub fn reciprocal_mod(&self) -> Self {
        let mut x = self.val;
        if !(0 <= x && x < self.modulo) {
            panic!("Invalid input parameters");
        }
        // Based on a simplification of the extended Euclidean algorithm
        let mut y = x;
        x = self.modulo;
        let mut a = 0;
        let mut b = 1;
        while y != 0 {
            let temp = a - x / y * b;
            a = b;
            b = temp;
            let temp = x % y;
            x = y;
            y = temp;
        }
        if x == 1 {
            FF::new(((a % self.modulo) + self.modulo) % self.modulo)
        } else {
            panic!("Reciprocal does not exist");
        }
    }
}

impl Add<FF> for FF {
    type Output = FF;

    fn add(mut self, rhs: Self::Output) -> Self::Output {
        self.val = (self.val + rhs.val) % self.modulo;
        self
    }
}

impl Mul<FF> for FF {
    type Output = FF;

    fn mul(mut self, rhs: Self::Output) -> Self::Output {
        self.val = (self.val * rhs.val) % self.modulo;
        self
    }
}
impl Sub<FF> for FF {
    type Output = FF;

    fn sub(mut self, rhs: Self::Output) -> Self::Output {
        self.val = (self.val - rhs.val) % self.modulo;
        self
    }
}

// ------------------------------------------------------------------------------

struct NTT<const M : usize > {
}

impl<const M : usize> NTT<M> {
    pub fn new () -> Self {
        NTT{}
    }

    pub fn ntt_recursive(&self, invec: &Vec<FF>, root: &FF) -> Vec<FF> {
        let n = invec.len();
        if n == 1 {
            return vec![invec[0].clone()];
        }

        let half_n = n / 2;

        // Separate even and odd coefficients
        let (even, odd): (Vec<_>, Vec<_>) =
            invec.iter().enumerate().partition(|(i, _)| *i % 2 == 0);
        let even: Vec<_> = even.into_iter().map(|(_, v)| v.clone()).collect();
        let odd: Vec<_> = odd.into_iter().map(|(_, v)| v.clone()).collect();

        // Recursively compute NTT on even and odd parts
        let root_squared = root.pow_mod(2);
        let even_ntt = self.ntt_recursive(&even, &root_squared);
        let odd_ntt = self.ntt_recursive(&odd, &root_squared);

        // Combine the results
        let mut outvec = vec![FF::new(0); n];
        let mut current_root = FF::new(1);
        for i in 0..half_n {
            let t = current_root.clone() * odd_ntt[i].clone();
            outvec[i] = even_ntt[i].clone() + t.clone();
            outvec[i + half_n] = even_ntt[i].clone() - t.clone();
            current_root = current_root.clone() * root.clone();
        }

        outvec
    }

    pub fn intt_recursive(&self, invec: &Vec<FF>, root: &FF) -> Vec<FF> {
        let mut outvec = self.ntt_recursive(invec, &(root.reciprocal_mod()));
        let scaler = FF::new(invec.len() as i64).reciprocal_mod();

        for i in 0..outvec.len() {
            outvec[i] = outvec[i].clone() * scaler.clone();
        }
        outvec
    }
}

#[cfg(test)]
mod tests1 {
    use super::*;
    #[test]
    fn test_example_ntt_struct() {
        let inputs: Vec<i64> = vec![11, 42, 31, 43, -11, 12, 78, 37];
        let mut max_val = 0;
        for &x in inputs.iter() {
            if x > max_val {
                max_val = x;
            }
        }

        let min_mod = max_val * max_val * inputs.len() as i64 + 1;
        let modulus = find_modulus(inputs.len(), min_mod as i128);
        let root = find_primitive_root(inputs.len() as i128, modulus - 1, modulus);

        const M :usize = 1;


        FF::set_shared_mod_value(modulus as i64);

        let ins : Vec<FF> = inputs.iter().map(|x| FF::new(x.clone())).collect();

        let rr = FF::new(root as i64);

        let ntt = NTT::<M>::new();
        let res = ntt.ntt_recursive(&ins, &rr);

        println!("Post ntt:   {:?}", res);
        for (u,v) in res.iter().zip(inputs.iter()) {
            assert_ne!(u.val, v.clone());
        }

        let res1 = ntt.intt_recursive(&res, &rr);
        println!("Post intt:   {:?}", res);
        println!("---");
        println!("Actual Result:   {:?}", res1);
        println!("Expected Result: {:?}", inputs.clone());
        for (u,v) in res1.iter().zip(inputs.iter()) {
            assert_eq!(u.val, v.clone());
        }
        assert!(false);
    }

}



// ------------------------------------------------------------------------------

// Returns the forward number-theoretic transform of the given vector with
// respect to the given primitive nth root of unity under the given modulus.

// Returns the forward number-theoretic transform of the given vector with
// respect to the given primitive nth root of unity under the given modulus.
pub fn ntt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len();
    if n == 1 {
        return vec![invec[0]];
    }

    let half_n = n / 2;

    // Separate even and odd coefficients
    let (even, odd): (Vec<_>, Vec<_>) = invec.iter().enumerate().partition(|(i, _)| *i % 2 == 0);
    let even: Vec<_> = even.into_iter().map(|(_, v)| v.clone()).collect();
    let odd: Vec<_> = odd.into_iter().map(|(_, v)| v.clone()).collect();

    // Recursively compute NTT on even and odd parts
    let root_squared = pow_mod(root, 2, modulus);
    let even_ntt = ntt_recursive(&even, root_squared, modulus);
    let odd_ntt = ntt_recursive(&odd, root_squared, modulus);

    // Combine the results
    let mut outvec = vec![0; n];
    let mut current_root = 1;
    for i in 0..half_n {
        let t = modulo_using_subtraction(current_root * odd_ntt[i], modulus);
        outvec[i] = modulo_using_subtraction(even_ntt[i] + t, modulus);
        outvec[i + half_n] = modulo_using_subtraction(even_ntt[i] - t, modulus);
        current_root = modulo_using_subtraction(current_root * root, modulus);
    }

    outvec
}

pub fn intt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let mut outvec = ntt_recursive(invec, reciprocal_mod(root, modulus), modulus);
    let scaler = reciprocal_mod(invec.len() as i128, modulus);
    for i in 0..outvec.len() {
        outvec[i] = modulo_using_subtraction(outvec[i] * scaler, modulus);
    }
    outvec
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_example_ntt_rec() {
        let inputs: Vec<i128> = vec![11, 42, 31, 43, -11, 12, 78, 37];
        let mut max_val = 0;
        for &x in inputs.iter() {
            if x > max_val {
                max_val = x;
            }
        }

        let min_mod = max_val * max_val * inputs.len() as i128 + 1;
        let modulus = find_modulus(inputs.len(), min_mod);
        let root = find_primitive_root(inputs.len() as i128, modulus - 1, modulus);

        println!("modulus:   {:?}", modulus.clone());
        println!("root: {:?}", root.clone());

        let temp0 = ntt_recursive(&inputs, root, modulus);

        let res = intt_recursive(&temp0, root, modulus);

        println!("Actual Result:   {:?}", res);
        println!("Expected Result: {:?}", inputs.clone());
        assert_eq!(res, inputs.clone());
    }

    #[test]
    fn test_example_ntt() {
        let inputs: Vec<i128> = vec![4, 1, 4, 2, 1, 3, 5, 6];
        let mut max_val = 0;
        for &x in inputs.iter() {
            if x > max_val {
                max_val = x;
            }
        }

        let min_mod = max_val * max_val * inputs.len() as i128 + 1;
        let modulus = find_modulus(inputs.len(), min_mod);
        let root = find_primitive_root(inputs.len() as i128, modulus - 1, modulus);

        println!("modulus:   {:?}", modulus.clone());
        println!("root: {:?}", root.clone());

        let temp0 = ntt_recursive(&inputs, root, modulus);

        let res = intt_recursive(&temp0, root, modulus);

        println!("Actual Result:   {:?}", res);
        println!("Expected Result: {:?}", inputs.clone());
        assert_eq!(res, inputs.clone());

        //  test_circular_convole() {
        //
        // Example vectors for testing
        let vec0: Vec<i128> = vec![4, 1, 4, 2, 1, 3, 5, 6];
        let vec1: Vec<i128> = vec![6, 1, 8, 0, 3, 3, 9, 8];

        // Circular convolution
        let result = circular_convolve(&vec0, &vec1);

        // Debugging information
        println!("Actual Result:   {:?}", result);
        println!(
            "Expected Result: {:?}",
            vec![123, 120, 106, 92, 139, 144, 140, 124]
        );

        // Assert the result
        assert_eq!(result, vec![123, 120, 106, 92, 139, 144, 140, 124]);
    }
}

fn modulo_using_subtraction(dividend: i128, divisor: i128) -> i128 {
    dividend % divisor
}

pub fn ntt(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len();
    let mut outvec = Vec::with_capacity(n);
    for i in 0..n {
        let mut sum = 0;
        for j in 0..n {
            let k = modulo_using_subtraction(i as i128 * j as i128, n as i128);
            sum = modulo_using_subtraction(sum + invec[j] * pow_mod(root, k, modulus), modulus);
        }
        outvec.push(sum);
    }
    outvec
}

// Returns the inverse number-theoretic transform of the given vector with
// respect to the given primitive nth root of unity under the given modulus.
pub fn inverse_ntt(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let mut outvec = ntt(invec, reciprocal_mod(root, modulus), modulus);
    let scaler = reciprocal_mod(invec.len() as i128, modulus);
    for i in 0..outvec.len() {
        outvec[i] = modulo_using_subtraction(outvec[i] * scaler, modulus);
    }
    outvec
}

// Returns the circular convolution of the given vectors of integers.
// All values must be non-negative. Internally, a sufficiently large modulus
// is chosen so that the convolved result can be represented without overflow.
pub fn circular_convolve(vec0: &Vec<i128>, vec1: &Vec<i128>) -> Vec<i128> {
    if vec0.len() == 0 || vec0.len() != vec1.len() {
        panic!("Invalid vector lengths");
    }

    let mut max_val = 0;
    for &x in vec0.iter().chain(vec1.iter()) {
        if x > max_val {
            max_val = x;
        }
    }

    let min_mod = max_val * max_val * vec0.len() as i128 + 1;
    let modulus = find_modulus(vec0.len(), min_mod);
    let root = find_primitive_root(vec0.len() as i128, modulus - 1, modulus);
    let temp0 = ntt(vec0, root, modulus);
    let temp1 = ntt(vec1, root, modulus);
    let mut temp2 = Vec::with_capacity(temp0.len());
    for i in 0..temp0.len() {
        temp2.push(modulo_using_subtraction(temp0[i] * temp1[i], modulus));
    }
    inverse_ntt(&temp2, root, modulus)
}

// Returns the smallest modulus mod such that mod = i * veclen + 1
// for some integer i >= 1, mod > veclen, and mod is prime.
// Although the loop might run for a long time and create arbitrarily large numbers,
// Dirichlet's theorem guarantees that such a prime number must exist.
pub fn find_modulus(vec_len: usize, minimum: i128) -> i128 {
    if vec_len < 1 || minimum < 1 {
        panic!("invalid input parameters");
    }
    let vl = vec_len as i128;
    let mut start = (minimum - 1 + vl - 1) / vl;
    if start < 1 {
        start = 1;
    }
    let mut n = start * vl + 1;
    loop {
        if is_prime(n) {
            return n;
        }
        n += vl;
    }
}

// Returns an arbitrary primitive degree-th root of unity modulo mod.
// totient must be a multiple of degree. If mod is prime, an answer must exist.
pub fn find_primitive_root(degree: i128, totient: i128, modulus: i128) -> i128 {
    if !(0 <= degree && degree <= totient && totient < modulus)
        || modulo_using_subtraction(totient, degree) != 0
    {
        panic!("Invalid input parameters");
    }
    let gen = find_generator(totient, modulus);
    pow_mod(gen, (totient / degree) as i128, modulus)
}

// Returns an arbitrary generator of the multiplicative group of integers modulo mod.
// totient must equal the Euler phi function of mod. If mod is prime, an answer must exist.
pub fn find_generator(totient: i128, modulus: i128) -> i128 {
    if !(1 <= totient && totient < modulus) {
        panic!("Invalid input parameters");
    }
    for i in 1..modulus {
        if is_primitive_root(i, totient, modulus) {
            return i;
        }
    }
    panic!("No generator exists");
}

pub fn is_primitive_root(val: i128, degree: i128, modulus: i128) -> bool {
    if !(0 <= val && val < modulus) {
        panic!("Invalid input parameters");
    }
    if !(1 <= degree && degree < modulus) {
        panic!("Invalid input parameters");
    }
    pow_mod(val, degree, modulus) == 1
        && unique_prime_factors(degree)
            .into_iter()
            .all(|p| pow_mod(val, degree / p, modulus) != 1)
}

// Returns a list of unique prime factors of the given integer in
// ascending order. For example, unique_prime_factors(60) = vec![2, 3, 5].
pub fn unique_prime_factors(mut n: i128) -> Vec<i128> {
    if n < 1 {
        panic!("Invalid input parameter");
    }
    let mut result = Vec::new();
    for i in 2.. {
        if i * i > n {
            break;
        }
        if modulo_using_subtraction(n, i) == 0 {
            result.push(i);
            while modulo_using_subtraction(n, i) == 0 {
                n /= i;
            }
        }
    }
    if n > 1 {
        result.push(n);
    }
    result
}

// Tests whether the given integer n >= 2 is a prime number.
pub fn is_prime(n: i128) -> bool {
    if n <= 1 {
        panic!("Invalid input parameter");
    }
    if n == 2 {
        return true;
    }
    if modulo_using_subtraction(n, 2) == 0 {
        return false;
    }
    let mut i = 3;
    while i * i <= n {
        if modulo_using_subtraction(n, i) == 0 {
            return false;
        }
        i += 2;
    }
    true
}

// Returns floor(sqrt(x)) for the given integer x >= 0.
pub fn sqrt(x: i128) -> i128 {
    if x < 0 {
        panic!("Invalid input parameter");
    }
    let mut y = 0;
    let mut i = 1 << (64 - x.leading_zeros() - 2);
    while i > 0 {
        y |= i;
        if y * y > x {
            y ^= i;
        }
        i >>= 1;
    }
    y
}

// Returns x^y mod m.
pub fn pow_mod(mut x: i128, mut y: i128, modulus: i128) -> i128 {
    if y < 0 || modulus <= 0 {
        panic!("Invalid input parameters");
    }
    if !(0 <= x && x < modulus) {
        x = modulo_using_subtraction(modulo_using_subtraction(x, modulus) + modulus, modulus);
    }
    let mut result = 1;
    while y != 0 {
        if y & 1 != 0 {
            result = modulo_using_subtraction(result * x, modulus);
        }
        x = modulo_using_subtraction(x * x, modulus);
        y >>= 1;
    }
    result
}

// Returns x^-1 mod m.
pub fn reciprocal_mod(mut x: i128, modulus: i128) -> i128 {
    if !(0 <= x && x < modulus) {
        panic!("Invalid input parameters");
    }
    // Based on a simplification of the extended Euclidean algorithm
    let mut y = x;
    x = modulus;
    let mut a = 0;
    let mut b = 1;
    while y != 0 {
        let temp = a - x / y * b;
        a = b;
        b = temp;
        let temp = modulo_using_subtraction(x, y);
        x = y;
        y = temp;
    }
    if x == 1 {
        modulo_using_subtraction(modulo_using_subtraction(a, modulus) + modulus, modulus)
    } else {
        panic!("Reciprocal does not exist");
    }
}
