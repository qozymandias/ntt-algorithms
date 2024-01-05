// 
// fn taylor_shift_horner_matrix(n: usize, a: &mut [f64], shift: f64) {
//     if shift == 0.0 {
//         return; // No shift, no change
//     }
// 
//     let m = n + 1;
//     let mut t: Vec<Vec<f64>> = vec![vec![0.0; m]; m];
// 
//     for i in 0..=n {
//         t[0][i] = a[i];
//         if i != 0 {
//             t[i][0] = a[0];
//         }
//     }
// 
//     for i in 1..=n {
//         a[i] = t[n - i + 1][i];
//     }
// 
//     for j in 1..=n {
//         for i in 1..=n - j + 1 {
//             t[j][i] = shift * t[j][i - 1] + t[j - 1][i];
//         }
//     }
// 
//     for i in 1..=n {
//         a[i] = t[n - i + 1][i];
//     }
// }
// 
// fn main() {
//     {
//     let mut a = vec![1.0, 2.0, 3.0, 4.0];
//     let shift = 2.0;
// 
//     taylor_shift_horner_matrix(a.len() - 1, &mut a, shift);
// 
//     println!("{:?}", a);
//     }
//     {
//     let mut a = vec![4.0, 3.0, 2.0, 1.0];
//     let shift = 2.0;
// 
//     taylor_shift_horner_matrix(a.len() - 1, &mut a, shift);
// 
//     println!("{:?}", a);
//     }
// }
// fn taylor_shift(n: usize, a: &mut [f64], shift: f64) {
//     if shift == 0.0 {
//         return; // No shift, no change
//     }
// 
//     let m = n + 1;
//     let mut b: Vec<f64> = vec![0.0; m];
//     b[0] = a[1] * shift.powi((n - 1) as i32);
// 
//     for i in 1..=n {
//         b[i] = a[0] * shift.powi(n as i32);
//     }
// 
//     for j in 1..=n {
//         for i in (1..=j).rev() {
//             b[i] += b[i - 1];
//         }
// 
//         if j == n {
//             b[0] = 0.0;
//         } else {
//             b[0] = a[j + 1] * shift.powi((n - j - 1) as i32);
//         }
//     }
// 
//     for i in 0..n {
//         a[n - i] = b[i + 1] / shift.powi(i as i32);
//     }
// }
// 
// fn main() {
//     let mut coefficients = vec![1.0, 2.0, 3.0, 4.0]; // Replace with your polynomial coefficients
//     let shift_value = 2.0;
// 
//     println!("Original Coefficients: {:?}", coefficients);
// 
//     taylor_shift(coefficients.len() - 1, &mut coefficients, shift_value);
// 
//     println!("Shifted Coefficients: {:?}", coefficients);
// }



// 
// 
// 
// fn taylor_shift(n: usize, a: &mut [f64], shift: f64) {
//     if shift == 0.0 {
//         return; // No shift, no change
//     }
// 
//     let m = n + 1;
//     let mut t: Vec<Vec<f64>> = vec![vec![0.0; m]; m];
// 
//     for i in 0..n {
//         t[i][0] = a[i + 1] * shift.powi((n - i - 1) as i32);
//         t[i][i + 1] = a[0] * shift.powi(n as i32);
//     }
// 
//     for j in 0..n {
//         for i in j + 1..=n {
//             t[i][j + 1] = t[i - 1][j] + t[i - 1][j + 1];
//         }
//     }
// 
//     for i in 0..n {
//         a[n - i] = t[n][i + 1] / shift.powi(i as i32);
//     }
// }
// 
// fn main() {
//     let mut coefficients = vec![1.0, 2.0, 3.0, 4.0]; // Replace with your polynomial coefficients
//     let shift_value = 2.0;
// 
//     println!("Original Coefficients: {:?}", coefficients);
// 
//     taylor_shift(coefficients.len() - 1, &mut coefficients, shift_value);
// 
//     println!("Shifted Coefficients: {:?}", coefficients);
// }







// 
// fn taylor_shift_horner_matrix(n: usize, a: &mut [f64], shift: f64) {
//     if shift == 0.0 {
//         return; // No shift, no change
//     }
// 
//     let m = n + 1;
//     let mut t: Vec<Vec<f64>> = vec![vec![0.0; m]; m];
// 
//     for i in 0..=n {
//         t[0][i] = a[i];
//         if i != 0 {
//             t[i][0] = a[0];
//         }
//     }
// 
//     for i in 1..=n {
//         a[i] = t[n - i + 1][i];
//     }
// 
//     for j in 1..=n {
//         for i in 1..=n - j + 1 {
//             t[j][i] = shift * t[j][i - 1] + t[j - 1][i];
//         }
//     }
// 
//     for i in 1..=n {
//         a[i] = t[n - i + 1][i];
//     }
// }
// 
// fn main() {
//     {
//     let mut a = vec![1.0, 2.0, 3.0, 4.0];
//     let shift = 2.0;
// 
//     taylor_shift_horner_matrix(a.len() - 1, &mut a, shift);
// 
//     println!("{:?}", a);
//     }
//     {
//     let mut a = vec![4.0, 3.0, 2.0, 1.0];
//     let shift = 2.0;
// 
//     taylor_shift_horner_matrix(a.len() - 1, &mut a, shift);
// 
//     println!("{:?}", a);
//     }
// }





//fn taylor_shift_horner(n: usize, a: &mut [f64], shift: f64) {
//    if shift == 0.0 {
//        return; // No shift, no change
//    }
//
//    for j in 1..=n {
//        for i in (1..=n - j + 1).rev() {
//            a[i] += shift * a[i - 1];
//        }
//    }
//}
//
//fn main() {
//    let mut a = vec![1.0, 2.0, 3.0, 4.0];
//    let shift = 2.0;
//
//    taylor_shift_horner(a.len() - 1, &mut a, shift);
//
//    println!("{:?}", a);
//}


// use std::vec::Vec;
// 
// fn factorial(n: usize) -> usize {
//     (1..=n).product()
// }
// 
// fn convolve(f: &Vec<i32>, h: &Vec<i32>) -> Vec<i32> {
//     let mut g = vec![0; f.len() + h.len() - 1];
//     for (hindex, &hval) in h.iter().enumerate() {
//         for (findex, &fval) in f.iter().enumerate() {
//             g[hindex + findex] += fval * hval;
//         }
//     }
//     g
// }
// 
// fn shift(f: &Vec<i32>, a: i32) -> Vec<i32> {
//     let n = f.len() - 1;
//     let u: Vec<i32> = (0..=n).map(|i| factorial(i) as i32 * f[n - i]).collect();
//     let v: Vec<i32> = (0..=n)
//         .map(|i| factorial(n) as i32 * a.pow(i as u32) / factorial(i) as i32)
//         .collect();
//     let g = convolve(&u, &v);
//     let result: Vec<i32> = (0..=n)
//         .map(|i| g[n - i] / (factorial(n) as i32 * factorial(i) as i32))
//         .collect();
//     result
// }
// 
// fn main() {
//     let f = vec![4,3,2,1];
//     let result = shift(&f, 2);
//     println!("{:?}", result);
// 
//     assert_eq!(result, [49, 62, 27, 4]);
// }

