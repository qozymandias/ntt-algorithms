use std::{
    fmt::Debug,
    ops::{Add, Index, Mul, Sub},
};

// ------------------------------------------------------------------------------

#[derive(Clone)]
struct FF {
    val: i128,
    modulo: i128,
}
impl Debug for FF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.val)
    }
}

// todo:
//    - make trait
//

impl FF {
    pub fn new(val: i128, modulo: i128) -> Self {
        FF {
            val: val,
            modulo: modulo,
        }
    }
    pub fn as_raw(&self) -> i128 {
        self.val
    }

    // Returns x^y mod m.
    pub fn pow_mod(&self, y: i128) -> Self {
        FF::new(self.val.pow(y as u32) % self.modulo, self.modulo)
    }

    // Returns x^-1 mod m.
    pub fn reciprocal_mod(&self) -> Self {
        fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
            if a == 0 {
                (b, 0, 1)
            } else {
                let (g, x, y) = extended_gcd(b % a, a);
                (g, y - (b / a) * x, x)
            }
        }
        fn mod_inverse(x: i128, m: i128) -> Option<i128> {
            let (g, x_inv, _) = extended_gcd(x, m);
            if g == 1 {
                Some((x_inv % m + m) % m)
            } else {
                None
            }
        }

        fn reciprocal_mod(x: i128, m: i128) -> Option<i128> {
            mod_inverse(x, m)
        }

        FF::new(reciprocal_mod(self.val, self.modulo).unwrap(), self.modulo)
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
        // ensure the result is non negative
        self.val = (self.val - rhs.val + self.modulo) % self.modulo;
        self
    }
}

// ------------------------------------------------------------------------------

struct NTT {
    modulo: i128,
}

impl NTT {
    pub fn new(modulus: i128) -> Self {
        NTT { modulo: modulus }
    }

    pub fn ntt_recursive(&self, invec: &Vec<FF>, root: &FF, n_calls: usize) -> Vec<FF> {
        println!("Recursive call number {:?}", n_calls);
        println!("current in vec = {:?}", invec);

        println!("root_curr = {:?}", root);

        let n = invec.len();
        if n == 1 {
            println!("Base case");
            return vec![invec[0].clone()];
        }

        let half_n = n / 2;

        // Separate even and odd coefficients
        let (even, odd): (Vec<_>, Vec<_>) =
            invec.iter().enumerate().partition(|(i, _)| *i % 2 == 0);
        let even: Vec<_> = even.into_iter().map(|(_, v)| v.clone()).collect();
        let odd: Vec<_> = odd.into_iter().map(|(_, v)| v.clone()).collect();

        let root_squared = root.pow_mod(2);
        println!("root_next = {:?}", root_squared);
        // Recursively compute NTT on even and odd parts
        let even_ntt = self.ntt_recursive(&even, &root_squared, n_calls + 1);
        let odd_ntt = self.ntt_recursive(&odd, &root_squared, n_calls + 1);

        println!("Looping at recursion level = {:?}", n_calls);
        println!("w_in = {:?}", root);
        // Combine the results
        let mut outvec = vec![FF::new(0, self.modulo); n];
        let mut current_root = FF::new(1, self.modulo);
        for i in 0..half_n {
            println!("w = w_in^{:?} = {:?}", i, current_root);
            let t = current_root.clone() * odd_ntt[i].clone();
            outvec[i] = even_ntt[i].clone() + t.clone();
            outvec[i + half_n] = even_ntt[i].clone() - t.clone();
            println!("outvec[i] = {:?}", outvec[i]);
            println!("outvec[i+half_n] = {:?}", outvec[i + half_n]);
            current_root = current_root.clone() * root.clone();
        }
        println!("Done");

        outvec
    }

    pub fn intt_recursive(&self, invec: &Vec<FF>, root: &FF) -> Vec<FF> {
        let mut outvec = self.ntt_recursive(invec, &(root.reciprocal_mod()), 0);
        let scaler = FF::new(invec.len() as i128, self.modulo).reciprocal_mod();

        for i in 0..outvec.len() {
            outvec[i] = outvec[i].clone() * scaler.clone();
            println!("outvec[i] * s = {:?}", outvec[i]);
        }
        outvec
    }

    fn factorial(&self, n: usize) -> FF {
        FF::new((1..=n).product::<usize>() as i128, self.modulo)
    }

    pub fn taylor_shift(&self, p: &Vec<FF>, k: FF, root: FF) -> Vec<FF> {
        // p is the input polynomial, k is the shift amount
        //      u_i = (n - i)! * p_{n-i}
        //      v_i = k^i * (1/i!)
        //      g_i = conv(u_i, v_i)
        //      q_i = g_{n-i} / i!

        let n = p.len() - 1;

        let u: Vec<FF> = (0..=n)
            .map(|i| self.factorial(n - i) * p[n - i].clone())
            .collect();
        let v: Vec<FF> = (0..=n)
            .map(|i| k.pow_mod(i as i128) * self.factorial(i).reciprocal_mod())
            .collect();

        let untt = self.ntt_recursive(&u, &root, 0);
        let vntt = self.ntt_recursive(&v, &root, 0);

        let mut g = Vec::with_capacity(untt.len());
        for i in 0..=n {
            g.push(untt[i].clone() * vntt[i].clone());
        }

        g = self.intt_recursive(&g, &root);

        let q: Vec<FF> = (0..=n)
            .map(|i| self.factorial(n).reciprocal_mod() * g[n - i].clone())
            .collect();
        q

        //   let n = p.len() - 1;
        //   let u = p
        //       .iter()
        //       .rev()
        //       .enumerate()
        //       .map(|(i, x)| x.clone() * self.factorial(n - i))
        //       .collect();
        //   let v = (0..=n)
        //       .map(|x| k.pow_mod(x as i128) * self.factorial(x).reciprocal_mod())
        //       .collect();

        //   let untt = self.ntt_recursive(&u, &root, 0);
        //   let vntt = self.ntt_recursive(&v, &root, 0);
        //   let mut g = (0..=n).map(|i| untt[i].clone() * vntt[i].clone()).collect();
        //   g = self.intt_recursive(&g, &root);

        //   let q: Vec<FF> = g
        //       .into_iter()
        //       .rev()
        //       .take(n + 1)
        //       .map(|x| x * self.factorial(n).reciprocal_mod())
        //       .collect();
        //   q
    }
}
// use rustfft::{num_complex::Complex, Fft, FftPlanner, num_traits::Zero};
// use std::cmp::max;
//
// fn fft_convolve(c1: Vec<Complex<f64>>, c2: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
//     let mut planner = FftPlanner::<f64>::new();
//     let fft = planner.plan_fft(c1.len() + c2.len());
//
//     let mut input1: Vec<Complex<f64>> = c1.into_iter().map(|x| x / (c1.len() as f64).sqrt()).collect();
//     let mut input2: Vec<Complex<f64>> = c2.into_iter().map(|x| x / (c2.len() as f64).sqrt()).collect();
//
//     fft.process(&mut input1);
//     fft.process(&mut input2);
//
//     let result: Vec<Complex<f64>> = input1.into_iter().zip(input2).map(|(a, b)| a * b).collect();
//
//     let mut inverse_fft = planner.plan_fft(result.len()).make_inverse();
//     inverse_fft.process(&mut result);
//
//     result.into_iter().map(|x| x / (result.len() as f64).sqrt()).collect()
// }
//
// fn tay_shift_fourier(coeffs: Vec<f64>, a: f64) -> Vec<f64> {
//     let nco = coeffs.len() - 1;
//
//     let convolved = fft_convolve(
//         coeffs.iter().rev().enumerate().map(|(i, &x)| x * (1..=nco - i).product::<usize>() as f64).collect(),
//         (0..=nco).map(|x| a.powi(x as i32) / (1..=x).product::<usize>() as f64).collect(),
//     );
//
//     convolved.into_iter().rev().take(nco + 1).map(|x| x / (1..=nco).product::<usize>() as f64).collect()
// }

// fn main() { // Example usage
// let coeffs = vec![1.0, 2.0, 3.0];
// let a = 2.0;
//
// let result = tay_shift_fourier(coeffs, a);
//
// println!("{:?}", result);
// }

#[cfg(test)]
mod tests1 {
    use super::*;
    #[test]
    fn test_convolve_4d() {
        // Example vectors for testing
        let vec0: Vec<i128> = vec![4, 1, 4, 2, 1, 3, 5, 6];
        let vec1: Vec<i128> = vec![6, 1, 8, 0, 3, 3, 9, 8];
        let expected_out: Vec<i128> = vec![123, 120, 106, 92, 139, 144, 140, 124];

        let mut max_val = 0;
        for &x in vec0.iter().chain(vec1.iter()) {
            if x > max_val {
                max_val = x;
            }
        }

        let min_mod = max_val * max_val * vec0.len() as i128 + 1;
        let modulus = find_modulus(vec0.len(), min_mod as i128);
        let root = find_primitive_root(vec0.len() as i128, modulus - 1, modulus);

        let ntt = NTT::new(modulus);

        let ins0: Vec<FF> = vec0
            .iter()
            .map(|x| FF::new(x.clone(), ntt.modulo))
            .collect();
        let ins1: Vec<FF> = vec1
            .iter()
            .map(|x| FF::new(x.clone(), ntt.modulo))
            .collect();

        let rr = FF::new(root as i128, ntt.modulo);

        let res0 = ntt.ntt_recursive(&ins0, &rr, 0);
        let res1 = ntt.ntt_recursive(&ins1, &rr, 0);
        println!("Post ntt:   {:?}", res0);
        println!("Post ntt:   {:?}", res1);

        let mut res = Vec::with_capacity(vec0.len());
        for i in 0..vec0.len() {
            res.push(res0[i].clone() * res1[i].clone());
        }

        let res2 = ntt.intt_recursive(&res, &rr);
        println!("Post intt:  {:?}", res2);

        for i in 0..res2.len() {
            println!("{:?} vs {:?}", res2[i], expected_out[i]);
            assert_eq!(res2[i].val, FF::new(expected_out[i], ntt.modulo).val);
        }

        // @TODO: convolve only works for 2d vectors, so it broken somewhere in the recursion
        // Assert the result
        // assert_eq!(res, expected_out);
    }

    #[test]
    fn test_example_ntt_struct() {
        let mut inputs: Vec<i128> = vec![0; 8];
        for i in 0..inputs.len() {
            inputs[i] = i as i128 + 1;
        }

        let mut max_val = 0;
        for &x in inputs.iter() {
            if x > max_val {
                max_val = x;
            }
        }

        let min_mod = max_val * max_val * inputs.len() as i128 + 1;
        let modulus = find_modulus(inputs.len(), min_mod as i128);
        let root = find_primitive_root(inputs.len() as i128, modulus - 1, modulus);

        let ntt = NTT::new(modulus);

        let ins: Vec<FF> = inputs
            .iter()
            .map(|x| FF::new(x.clone(), ntt.modulo))
            .collect();
        let rr = FF::new(root as i128, ntt.modulo);

        let res = ntt.ntt_recursive(&ins, &rr, 0);

        println!("Post ntt:   {:?}", res);
        for (u, v) in res.iter().zip(inputs.iter()) {
            assert_ne!(u.val, v.clone());
        }

        let res1 = ntt.intt_recursive(&res, &rr);
        println!("---");
        println!("Post intt:   {:?}", res1);
        println!("---");
        println!("modulus = {:?}", modulus);
        println!("root = {:?}", root);
        println!("---");
        println!("Actual Result:   {:?}", res1);
        println!("Expected Result: {:?}", inputs.clone());
        println!("---");

        for (u, v) in res1.iter().zip(inputs.iter()) {
            assert_eq!(u.val, v.clone());
            assert_eq!(u.modulo, modulus.clone());
        }
    }

    #[test]
    fn test_example_ntt_taylor_shift() {
        let inputs: Vec<i128> = vec![1, 2, 3, 4];

        let mut max_val = 0;
        for &x in inputs.iter() {
            if x > max_val {
                max_val = x;
            }
        }

        let min_mod = max_val * max_val * inputs.len() as i128 + 1;
        let modulus = find_modulus(inputs.len(), min_mod as i128);
        let root = find_primitive_root(inputs.len() as i128, modulus - 1, modulus);

        let ntt = NTT::new(modulus);

        let ins: Vec<FF> = inputs
            .iter()
            .map(|x| FF::new(x.clone(), ntt.modulo))
            .collect();
        let rr = FF::new(root as i128, ntt.modulo);

        let shift = FF::new(2, ntt.modulo);

        let res = ntt.taylor_shift(&ins, shift, rr);

        println!("taylor_shift:   {:?}", res);
        for (u, v) in res.iter().zip(inputs.iter()) {
            assert_ne!(u.val, v.clone());
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
