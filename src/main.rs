mod utils {
    extern crate libc;
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

    /// Straightforward quadratic solver
    pub fn roots(a: f64, b: f64, c: f64) -> (f64, f64) {
        let d = b.powf(2.0) - (4.0 * a * c);
        if d < 0.0 {
            panic!("This eq has no real solution!");
        } else if d == 0.0 {
            let x = (-b + d.sqrt()) / (2.0 * a);
            return (x, x);
        } else {
            let x1 = (-b + d.sqrt()) / (2.0 * a);
            let x2 = (-b - d.sqrt()) / (2.0 * a);
            return (x1, x2);
        }
    }

    // ------------------------------------------------------------------------- //
    //                             Bin Implementation                            //
    // ------------------------------------------------------------------------- //

    /// A Bin with a given mean and count:
    /// value: The mean of the bin.
    /// count: The number of points in this bin. It is assumed that there are
    ///        `count` points surrounding `value`, of which `count/2` points are
    ///        to the left and `count/2` points are to the right.
    #[derive(Copy, Clone)]
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
    println!("Hello, world!");
    println!("{}", utils::next_after(1.1, 0.1));
    println!("{}", utils::E);

    let a = vec![5.0, 2.0, 3.0, 4.0];
    println!("{:?}", utils::argmin(&a));
    println!("{}", utils::argmin(&a).0);
    let b = utils::Bin {
        value: 1.0,
        count: 0,
    };
    let c = utils::Bin {
        value: 4.0,
        count: 0,
    };

    println!("diff:{}", utils::diff(b, c, true));
    println!("{}", b);
    println!("{}", b.bin_equality(&c));
    println!("{:?}", utils::linspace(0, 10, 10));
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
}
