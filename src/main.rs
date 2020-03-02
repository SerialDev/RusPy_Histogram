mod utils {
    extern crate libc;
    pub use std::f64::consts::E;
    use std::fmt;
    use std::ops::Index;

    mod cmath {
        use libc::c_double;

        extern "C" {
            pub fn nextafter(x: c_double, y: c_double) -> c_double;
        }
    }

    /// Returns the next floating-point number after x in the direction of y
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

    // ------------------------------------------------------------------------- //
    //                             Bin Implementation                            //
    // ------------------------------------------------------------------------- //

    /// A Bin with a given mean and count:
    /// value: The mean of the bin.
    /// count: The number of points in this bin. It is assumed that there are
    ///        `count` points surrounding `value`, of which `count/2` points are
    ///        to the left and `count/2` points are to the right.
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

    // impl<T> Index<T> for Bin {
    //     type Output = &T;

    //     fn index<'a>(&self, bin: usize) -> Self::Output {
    //         // unimplemented!("");
    //         !match bin {
    //             0 => &self.value,
    //             1 => &self.count,
    //         }
    //     }
    // }
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

    println!("{}", b);
    println!("{}", b.bin_equality(&c));
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
}
