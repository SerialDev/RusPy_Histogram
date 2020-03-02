mod utils {
    extern crate libc;
    pub use std::f64::consts::E;
    use std::fmt;

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
    }

    impl fmt::Display for Bin {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "Bin(value={}, count={})", self.value, self.count)
        }
    }

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
    println!("{}", b);
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
}
