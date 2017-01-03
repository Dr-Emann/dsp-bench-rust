extern crate rand;
extern crate time;


use std::f64::consts::PI;
use rand::Rng;
use time::precise_time_ns;
use std::env;
use std::fmt;

#[derive(Copy, Clone, PartialEq)]
struct Biquad {
    b0: f64,
    b1: f64,
    b2: f64,
    a1: f64,
    a2: f64,

    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
}

impl fmt::Debug for Biquad {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(writeln!(f, "Biquad coefficents:"));
        try!(writeln!(f, "b0={:+.13},", self.b0));
        try!(writeln!(f, "b1={:+.13},", self.b1));
        try!(writeln!(f, "b2={:+.13},", self.b2));
        try!(writeln!(f, "a1={:+.13},", self.a1));
        try!(write!  (f, "a2={:+.13},", self.a2));
        Ok(())
    }
}


fn peak_eq_coefs(fs: f64, f0: f64, q: f64, db_gain: f64) -> Biquad {
    let a = 10.0_f64.powf(db_gain / 40.0);
    let omega = 2.0 * PI * f0 / fs;

    let alpha = omega.sin() / (2.0 * q);

    let a0 = 1.0 + alpha / a;

    let b0 = (1.0 + alpha * a) / a0;
    let b1 = (-2.0 * omega.cos()) / a0;
    let b2 = (1.0 - alpha * a) / a0;
    let a2 = (1.0 - alpha / a) / a0;

    Biquad {
        b0: b0,
        b1: b1,
        b2: b2,
        a1: b1,
        a2: a2,
        x1: 0.0,
        x2: 0.0,
        y1: 0.0,
        y2: 0.0,
    }
}

fn white_noise(length: usize) -> Vec<f64> {
    let mut rng = rand::thread_rng();
    let mut vec: Vec<f64> = Vec::new();
    vec.resize(length, 0.0);

    for i in 0..vec.len() {
        vec[i] = rng.gen::<f64>()
    }

    return vec;
}


fn white_noise_array() -> [f64; 4096] {
    let mut rng = rand::thread_rng();

    let mut arr = [0.0; 4096];
    for i in 0..arr.len() {
        arr[i] = rng.gen::<f64>()
    }

    return arr;
}


fn iir_vec(input: &Vec<f64>, output: &mut Vec<f64>, bq: &mut Biquad) {
    for i in 0..input.len() {
        output[i] = (bq.b0 * input[i]) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) - (bq.a1 * bq.y1) -
                    (bq.a2 * bq.y2);

        bq.x2 = bq.x1;
        bq.x1 = input[i];

        bq.y2 = bq.y1;
        bq.y1 = output[i];
    }
}

fn iir_vec_zip(input: &Vec<f64>, output: &mut Vec<f64>, bq: &mut Biquad) {

    for (x, y) in input.iter().zip(output.iter_mut()) {
        *y = (bq.b0 * *x) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) - (bq.a1 * bq.y1) - (bq.a2 * bq.y2);

        bq.x2 = bq.x1;
        bq.x1 = *x;

        bq.y2 = bq.y1;
        bq.y1 = *y;
    }
}

fn iir_vec_zip_2(input: &Vec<f64>, output: &mut Vec<f64>, bq: &mut Biquad) {

    for (&x, y) in input.iter().zip(output.iter_mut()) {
        *y = (bq.b0 * x) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) - (bq.a1 * bq.y1) - (bq.a2 * bq.y2);

        bq.x2 = bq.x1;
        bq.x1 = x;

        bq.y2 = bq.y1;
        bq.y1 = *y;
    }
}

fn iir_vec_2(input: &Vec<f64>, output: &mut Vec<f64>, bq: &mut Biquad) {
    let len = input.len();

    for i in 0..len {
        let inval = input[i];

        let outval = (bq.b0 * inval) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) - (bq.a1 * bq.y1) -
                     (bq.a2 * bq.y2);

        output[i] = outval;

        bq.x2 = bq.x1;
        bq.x1 = inval;

        bq.y2 = bq.y1;
        bq.y1 = outval;
    }
}

fn iir_array(input: &[f64; 4096], output: &mut [f64; 4096], bq: &mut Biquad) {
    for i in 0..input.len() {
        output[i] = (bq.b0 * input[i]) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) - (bq.a1 * bq.y1) -
                    (bq.a2 * bq.y2);

        bq.x2 = bq.x1;
        bq.x1 = input[i];

        bq.y2 = bq.y1;
        bq.y1 = output[i];
    }
}

fn iir_slice(input: &[f64], output: &mut [f64], bq: &mut Biquad) {
    for i in 0..input.len() {
        output[i] = (bq.b0 * input[i]) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) - (bq.a1 * bq.y1) -
                    (bq.a2 * bq.y2);

        bq.x2 = bq.x1;
        bq.x1 = input[i];

        bq.y2 = bq.y1;
        bq.y1 = output[i];
    }
}

fn iir_slice_unsafe(input: &[f64], output: &mut [f64], bq: &mut Biquad) {
    unsafe {
        for i in 0..input.len() {
            output[i] = (bq.b0 * *input.get_unchecked(i)) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) -
                        (bq.a1 * bq.y1) - (bq.a2 * bq.y2);

            bq.x2 = bq.x1;
            bq.x1 = *input.get_unchecked(i);

            bq.y2 = bq.y1;
            bq.y1 = *output.get_unchecked(i);
        }
    }
}

fn main() {
    println!("DSP bench rust");

    let buffer_length = env::args().nth(1).and_then(|arg| arg.parse::<usize>().ok()).unwrap_or(4096);


    let bench_loops = 200000;

    let mut bq = peak_eq_coefs(48000.0, 200.0, 2.0, 6.0);
    println!("{:?}", bq);

    let input = white_noise(buffer_length);
    let mut output: Vec<f64> = vec![0.0; input.len()];

    {
        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_vec(&input, &mut output, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_vec:\t\t\t{} ns per loop", elapsed / bench_loops);
    }

    {
        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_vec_zip(&input, &mut output, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_vec_zip:\t\t\t{} ns per loop", elapsed / bench_loops);
    }

    {
        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_vec_zip_2(&input, &mut output, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_vec_zip_2:\t\t\t{} ns per loop", elapsed / bench_loops);
    }

    {
        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_vec_2(&input, &mut output, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_vec_2:\t\t\t{} ns per loop", elapsed / bench_loops);
    }

    {
        let input_array = white_noise_array();
        let mut output_array = [0.0; 4096];

        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_array(&input_array, &mut output_array, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_array (4096):\t\t{} ns per loop", elapsed / bench_loops);
    }

    {
        let input_array = white_noise_array();
        let mut output_array = [0.0; 4096];

        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_slice(&input_array[0..buffer_length],
                      &mut output_array[0..buffer_length],
                      &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_slice (array):\t\t{} ns per loop",
                 elapsed / bench_loops);
    }

    {
        let input_array = white_noise_array();
        let mut output_array = [0.0; 4096];

        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_slice(&input_array, &mut output_array, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_slice (4096 array):\t\t{} ns per loop",
                 elapsed / bench_loops);
    }

    {
        let input_array = white_noise_array();
        let mut output_array = [0.0; 4096];

        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_slice_unsafe(&input_array, &mut output_array, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_slice_unsafe (4096 array):\t{} ns per loop",
                 elapsed / bench_loops);
    }

    {
        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_slice(&input, &mut output, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_slice (vec):\t\t{} ns per loop", elapsed / bench_loops);
    }

    {
        let start = precise_time_ns();
        for _ in 0..bench_loops {
            iir_slice_unsafe(&input, &mut output, &mut bq);
        }
        let elapsed = precise_time_ns() - start;
        println!("iir_slice_unsafe (vec):\t\t{} ns per loop",
                 elapsed / bench_loops);
    }

}
