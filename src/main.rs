extern crate libc;

mod cmath {
    use libc::c_double;

    extern "C" {
        pub fn nextafter(x: c_double, y: c_double) -> c_double;
    }
}

/// Returns the next floating-point number after x in the direction of y
fn next_after(x: f64, y: f64) -> f64 {
    println!("");
    unsafe { cmath::nextafter(x, y) }
}

fn main() {
    println!("Hello, world!");
    println!("{}", next_after(1.1, 0.1))
}

// _E = 2.718281828459045
