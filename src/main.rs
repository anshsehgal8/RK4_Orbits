use std::ops::Mul;
use std::ops::Add;
use hdf5::{File, Group, H5Type};
use kepler_two_body::{OrbitalElements, OrbitalState};

static GRAVITATIONAL_CONSTANT_G: f64 = 1.0; //6.67E-2;
static MASS1:f64 = 1.0;
static MASS2:f64 = 1.0;

// ============================================================================
#[derive(Copy, Clone)]

struct TwoBodyState
{
    x1:  f64,
    y1:  f64,
    vx1: f64,
    vy1: f64,
    x2:  f64,
    y2:  f64,
    vx2: f64,
    vy2: f64,
}


#[derive(Copy, Clone)]
struct DerivativeState{
    vx1: f64,
    vy1: f64,
    ax1: f64,
    ay1: f64,
    vx2: f64,
    vy2: f64,
    ax2: f64,
    ay2: f64,
}


impl DerivativeState{
    fn to_twobodystate(&self) -> TwoBodyState {
        TwoBodyState{
            x1:  self.vx1,
            y1:  self.vy1,
            vx1: self.ax1,
            vy1: self.ay1,
            x2:  self.vx2,
            y2:  self.vy2,
            vx2: self.ax2,
            vy2: self.ay2,
            
        }
    }
}



impl Mul<f64> for TwoBodyState
{
    type Output = TwoBodyState;
    fn mul(self, scale: f64) -> TwoBodyState
    {
        TwoBodyState
        {
            x1:  self.x1  * scale,
            x2:  self.x2  * scale,
            vx1: self.vx1 * scale,
            vx2: self.vx2 * scale,
            y1:  self.y1  * scale,
            y2:  self.y2  * scale,
            vy1: self.vy1 * scale,
            vy2: self.vy2 * scale,
        }
    }
}



impl Mul<f64> for DerivativeState
{
    type Output = DerivativeState;
    fn mul(self, scale: f64) -> DerivativeState
    {
        DerivativeState
        {
            vx1: self.vx1 * scale,
            vy1: self.vy1 * scale,
            ax1: self.ax1 * scale,
            ay1: self.ay1 * scale,
            vx2: self.vx2 * scale,
            vy2: self.vy2 * scale,
            ax2: self.ax2 * scale,
            ay2: self.ay2 * scale,
        }
    }
}



impl Add for TwoBodyState
{
    type Output = TwoBodyState;
    fn add(self, other: TwoBodyState) -> TwoBodyState
    {
        TwoBodyState
        {
            x1:  self.x1  + other.x1,
            y1:  self.y1  + other.y1,
            vx1: self.vx1 + other.vx1,
            vy1: self.vy1 + other.vy1,
            x2:  self.x2  + other.x2,
            y2:  self.y2  + other.y2,
            vx2: self.vx2 + other.vx2,
            vy2: self.vy2 + other.vy2,
        }
    }
}



fn derivative(state_vector: &TwoBodyState) -> DerivativeState{
    let x1  = state_vector.x1;
    let y1  = state_vector.y1;
    let vx1 = state_vector.vx1;
    let vy1 = state_vector.vy1;
    let x2  = state_vector.x2;
    let y2  = state_vector.y2;
    let vx2 = state_vector.vx2;
    let vy2 = state_vector.vy2;  
    let r   = ((x1-x2).powf(2.0) + (y1-y2).powf(2.0)).powf(0.5);
    let ax1 = -(GRAVITATIONAL_CONSTANT_G * MASS2 / (r.powi(3) )) * (x1-x2);
    let ay1 = -(GRAVITATIONAL_CONSTANT_G * MASS2 / (r.powi(3) )) * (y1-y2);
    let ax2 = -(GRAVITATIONAL_CONSTANT_G * MASS1 / (r.powi(3) )) * (x2-x1);
    let ay2 = -(GRAVITATIONAL_CONSTANT_G * MASS1 / (r.powi(3) )) * (y2-y1);

    DerivativeState{
        vx1: vx1,
        vy1: vy1,
        ax1: ax1,
        ay1: ay1,
        vx2: vx2,
        vy2: vy2,
        ax2: ax2,
        ay2: ay2,
    }
}




fn rk4(state_vector: TwoBodyState, dt: f64) -> TwoBodyState
{
    let k1 = derivative(&state_vector).to_twobodystate() * dt;
    let a2 = state_vector + k1 * 0.5;
    let k2 = derivative(&a2).to_twobodystate() * dt;
    let a3 = state_vector + k2 * 0.5;
    let k3 = derivative(&a3).to_twobodystate() * dt;
    let a4 = state_vector + k3 * 1.0;
    let k4 = derivative(&a4).to_twobodystate() * dt;
    let u = state_vector + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (1.0/6.0);
    u
}


fn main()
{
    let mut s1 = TwoBodyState{
        x1:  -0.5,
        y1:   0.0,
        vx1:  0.0,
        vy1: -0.5,
        x2:   0.5,
        y2:   0.0,
        vx2:  0.0,
        vy2:  0.5,
    };
    let mut t = 0.0;
    let dt = 0.01;
    let mut x1_results = Vec::new();
    let mut y1_results = Vec::new();
    let mut x2_results = Vec::new();
    let mut y2_results = Vec::new();
    while t < 20.0 {
        s1 = rk4(s1,dt);
        println!("{:+0.8} {:+0.8} {:+0.8} {:+0.8} {:+0.8}", t, s1.x1, s1.y1, s1.x2, s1.y2);
        x1_results.push(s1.x1);
        y1_results.push(s1.y1);
        x2_results.push(s1.x2);
        y2_results.push(s1.y2);
        t = t + dt;
    }
    
    let file = File::create("output.h5").unwrap();
    let group = file.create_group("data").unwrap();
    group.new_dataset::<f64>().create("x1",x1_results.len()).unwrap().write_raw(&x1_results);
    group.new_dataset::<f64>().create("y1",y1_results.len()).unwrap().write_raw(&y1_results);
    group.new_dataset::<f64>().create("x2",x2_results.len()).unwrap().write_raw(&x2_results);
    group.new_dataset::<f64>().create("y2",y2_results.len()).unwrap().write_raw(&y2_results);




}