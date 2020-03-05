#[allow(dead_code)]
mod utils {
    extern crate libc;
    use itertools::izip;
    pub use std::f64::consts::E;
    use std::fmt;
    use std::ops::Add;

    mod cmath {
        use libc::c_double;

        extern "C" {
            pub fn nextafter(x: c_double, y: c_double) -> c_double;
        }
    }

    /// Returns the next floating-point number after x in the direction of y
    /// uses an unsafe call to C::Math
    pub fn next_after(x: f64, y: f64) -> f64 {
        println!("");
        unsafe { cmath::nextafter(x, y) }
    }

    /// Find argmin of slice
    /// Returns index and value of first ocurring minimum
    /// ```
    /// let a = vec![5.0, 2.0, 3.0, 4.0];
    /// let c = utils::argmin(&a);
    /// assert_eq!(c.0, 1);
    /// assert_eq!(c.1, 2.0);
    /// ```
    pub fn argmin<T>(u: &[T]) -> (usize, T)
    where
        T: Copy + PartialOrd,
    {
        assert!(u.len() != 0);

        let mut min_index = 0;

        let mut min = u[min_index];

        for (i, v) in u.iter().enumerate().skip(1) {
            if min > *v {
                min_index = i;
                min = *v
            }
        }
        (min_index, min)
    }

    /// Custom version of np's linspace
    /// returns evenly spaced numbers over a specified interval
    /// returns num evenly spaced samples calculated over the interval [start, stop]
    pub fn linspace(start: i64, stop: i64, num: i64) -> std::vec::Vec<f64> {
        if num == 1 {
            return vec![stop as f64];
        }
        let h = (stop as f64 - start as f64) / num as f64;
        let mut values = vec![];
        for i in 0..(num + 1) {
            values.push(start as f64 + h as f64 * i as f64);
        }
        return values;
    }

    pub fn diff(a: Bin, b: Bin, weighted: bool) -> f64 {
        let mut diff = b.value - a.value;
        if weighted {
            diff *= (E + a.count.min(b.count) as f64).ln();
        }
        return diff;
    }

    pub fn bin_diff(data: Vec<Bin>, weighted: bool) -> Vec<f64> {
        let mut result = vec![];
        for window in data.windows(2) {
            result.push(diff(window[0], window[1], weighted));
        }
        return result;
    }

    /// PERF: Make some modifications so that no floating point comparisons are done
    pub fn bin_sums(data: Vec<Bin>, less: Option<i64>) -> Vec<f64> {
        let mut result = vec![];
        for a in data.windows(2) {
            match less {
                Some(x) if a[1].value <= x as f64 => {
                    println!("-- {:?}{:?} --1\n", result, a);
                    result.push((a[0].count + a[1].count) as f64 / 2.0)
                }
                None => {
                    println!("-- {:?}{:?} --2\n", result, a);
                    result.push((a[0].count + a[1].count) as f64 / 2.0)
                }
                _ => (),
            }
        }
        return result;
    }

    /// Straightforward quadratic solver
    pub fn roots(a: f64, b: f64, c: f64) -> [f64; 2] {
        let d = b.powf(2.0) - (4.0 * a * c);
        if d < 0.0 {
            panic!("This eq has no real solution!");
        } else if d == 0.0 {
            let x = (-b + d.sqrt()) / (2.0 * a);
            return [x, x];
        } else {
            let x1 = (-b + d.sqrt()) / (2.0 * a);
            let x2 = (-b - d.sqrt()) / (2.0 * a);
            return [x1, x2];
        }
    }

    pub fn find_z(a: f64, b: f64, c: f64) -> f64 {
        let candidate_roots = roots(a, b, c);
        let mut result_root: f64 = 0.0;
        for candidate_root in candidate_roots.iter() {
            println!("candidate_root:{:?}", candidate_root);
            if *candidate_root >= 0.0 && *candidate_root <= 1.0 {
                println!("candidate_root:{:?}", candidate_root);
                result_root = *candidate_root;
                break;
            }
        }
        return result_root;
    }

    pub fn compute_quantile(x: f64, bin_i: Bin, bin_il: Bin, prev_sum: f64) -> f64 {
        let d = x - prev_sum;
        let a = bin_il.count - bin_i.count;
        if a == 0 {
            let offset = d / (bin_i.count + bin_il.count) as f64 / 2.0;
            let u = bin_i.value + (offset * (bin_il.value - bin_i.value));
            return u;
        } else {
            let b = 2.0 * bin_i.count as f64;
            let c = -2.0 * d;
            let z = find_z(a as f64, b, c);
            let u = bin_i.value + (bin_il.value - bin_i.value) * z;
            return u;
        }
    }

    ///Finding the density starting from the sum.
    /// s = p + (1/2 + r - r^2/2)*i + r^2/2*i1
    /// r = (x - m) / (m1 - m)
    /// s_dx = i - (i1 - i) * (x - m) / (m1 - m)
    pub fn compute_density(p: f64, bin_i: Bin, bin_il: Bin) -> f64 {
        let b_diff = p - bin_i.value;
        let p_diff = bin_il.value - bin_i.value;
        let bp_ratio = b_diff / p_diff;

        let inner = (bin_il.count - bin_i.count) as f64 * bp_ratio;
        let result = (bin_i.count as f64 + inner) * (1.0 / p_diff);
        return result;
    }

    /// Compute sum -- WIP details
    pub fn compute_sum(x: f64, bin_i: Bin, bin_il: Bin, prev_sum: f64) -> f64 {
        let b_diff = x - bin_i.value;
        let p_diff = bin_il.value - bin_i.value;
        let bp_ratio = b_diff / p_diff;

        let ilTerm = 0.5 * bp_ratio.powf(2.0);
        let iTerm = bp_ratio - ilTerm;

        let first = prev_sum + bin_i.count as f64 * iTerm;
        let ss = first + bin_il.count as f64 * ilTerm;

        return ss;
    }

    // ------------------------------------------------------------------------- //
    //                     Streaming Histogram implementation                    //
    // ------------------------------------------------------------------------- //
    #[derive()]
    pub struct StreamHist {
        pub maxbins: i64,
        pub bins: Vec<f64>,
        pub total: i64,
        pub weighted: bool,
        pub min: f64,
        pub max: f64,
        pub freeze: bool,
        pub missing_count: i64,
    }

    impl StreamHist {}

    // ------------------------------------------------------------------------- //
    //                             Bin Implementation                            //
    // ------------------------------------------------------------------------- //

    /// A Bin with a given mean and count:
    /// value: The mean of the bin.
    /// count: The number of points in this bin. It is assumed that there are
    ///        `count` points surrounding `value`, of which `count/2` points are
    ///        to the left and `count/2` points are to the right.
    #[derive(Debug, Copy, Clone)]
    pub struct Bin {
        pub value: f64,
        pub count: i64,
    }

    impl Bin {
        pub fn from_dict() -> Bin {
            return Bin {
                value: 1.0,
                count: 0i64,
            };
        }

        /// Test wheter two Bins are equal,
        /// comparison is done through the mean  
        pub fn bin_equality(&self, y: &Bin) -> bool {
            return self.value == y.value;
        }

        /// Test whether this bin has a lower mean
        /// than another Bin
        pub fn lower_than(&self, y: &Bin) -> bool {
            return self.value < y.value;
        }

        /// Test whether this bin has a higher mean
        /// than another Bin
        pub fn greater_than(&self, y: &Bin) -> bool {
            return self.value > y.value;
        }
    }

    /// Merge this bin with another bin and return the result
    /// This method implements Step 7 from Algorithm 1 in ref:
    /// Ben-Haim & Tom-Tov's Streaming Parallel Decision Tree Algorithm
    impl Add for Bin {
        type Output = Bin;
        fn add(self, other: Bin) -> Self::Output {
            // Sum Heights
            let count: i64 = self.count + other.count;
            // Weighted average
            let mut value = self.value * self.count as f64 + other.value + other.count as f64;
            value /= count as f64;
            return Bin {
                value: value,
                count: count,
            };
        }
    }

    /// Simple representation of a histogram bin.
    /// where value and count are the bin's
    /// stored mean and count.
    impl fmt::Display for Bin {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "Bin(value={}, count={})", self.value, self.count)
        }
    }

    /// Iterate over the mean and count of this bin
    impl Iterator for Bin {
        type Item = Bin;
        fn next(&mut self) -> Option<Self::Item> {
            unimplemented!("");
        }
    }
}

fn main() {
    let a = vec![5.0, 2.0, 3.0, 4.0];
    let b = utils::Bin {
        value: 1.0,
        count: 10,
    };
    let c = utils::Bin {
        value: 4.0,
        count: 20,
    };
    println!("{:?}", utils::find_z(-1.5, -1.0, 1.0));
    println!("compute_sum: {:?}", utils::compute_sum(1.5, b, c, 0.0));
    println!(
        "compute_quantile: {:?}",
        utils::compute_quantile(3.0, b, c, 0.0)
    );
}

#[cfg(test)]
mod tests {
    #[test]
    fn arg_min() {
        let a = vec![5.0, 2.0, 3.0, 4.0];
        let c = super::utils::argmin(&a);
        assert_eq!(c.0, 1);
        assert_eq!(c.1, 2.0);
    }

    #[test]
    fn bin_eq() {
        let a = super::utils::Bin {
            value: 1.0,
            count: 0,
        };
        let b = super::utils::Bin {
            value: 1.0,
            count: 0,
        };
        let c = super::utils::Bin {
            value: 4.0,
            count: 0,
        };

        assert_eq!(a.bin_equality(&b), true);
        assert_eq!(b.bin_equality(&c), false);
    }

    #[test]
    fn bin_lt() {
        let a = super::utils::Bin {
            value: 1.0,
            count: 0,
        };
        let b = super::utils::Bin {
            value: 1.0,
            count: 0,
        };
        let c = super::utils::Bin {
            value: 4.0,
            count: 0,
        };

        assert_eq!(a.lower_than(&b), false);
        assert_eq!(b.lower_than(&c), true);
    }

    #[test]
    fn bin_gt() {
        let a = super::utils::Bin {
            value: 1.0,
            count: 0,
        };
        let b = super::utils::Bin {
            value: 1.0,
            count: 0,
        };
        let c = super::utils::Bin {
            value: 4.0,
            count: 0,
        };

        assert_eq!(a.greater_than(&b), false);
        assert_eq!(c.greater_than(&b), true);
    }

    #[test]
    fn linspace_test() {
        let a = super::utils::linspace(0, 10, 10);
        assert_eq!(
            a,
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        );
    }

    #[test]
    fn bin_sums_test() {
        let b = super::utils::Bin {
            value: 1.0,
            count: 10,
        };
        let c = super::utils::Bin {
            value: 4.0,
            count: 20,
        };

        let result_none = super::utils::bin_sums(vec![b, c, b, b, c], None::<i64>);
        let result_lt = super::utils::bin_sums(vec![b, c, b, b, c], Some(100));
        let result_gt = super::utils::bin_sums(vec![b, c, b, b, c], Some(0));
        assert_eq! {result_none, vec![15.0, 15.0, 10.0, 15.0]}
        assert_eq! {result_lt, vec![15.0, 15.0, 10.0, 15.0]}
        assert_eq! {result_gt, vec![]}
    }

    #[test]
    fn bin_diff_test() {
        let b = super::utils::Bin {
            value: 1.0,
            count: 10,
        };
        let c = super::utils::Bin {
            value: 4.0,
            count: 20,
        };

        let result_true = super::utils::bin_diff(vec![b, c, b, b, c], true);
        let result_false = super::utils::bin_diff(vec![b, c, b, b, c], false);
        assert_eq! {result_true, vec![7.629121417228003, -7.629121417228003, 0.0, 7.629121417228003]}
        assert_eq! {result_false, vec![3.0, -3.0, 0.0, 3.0]}
    }
}
