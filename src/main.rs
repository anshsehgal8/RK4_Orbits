use std::ops::Mul;
use std::ops::Add;
//use std::num::Float;

static GRAVITATIONAL_CONSTANT_G: f64 = 6.67E-8;
static MASS_RATIO:f64 = 0.5;

// ============================================================================

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




impl Mul<f64> for TwoBodyState
{
    type Output = TwoBodyState;
    fn mul(self, scale: f64) -> TwoBodyState
    {
        TwoBodyState
        {
            x1: self.x1 * scale,
            x2: self.x2 * scale,
            vx1: self.vx1 * scale,
            vx2: self.vx2 * scale,
            y1: self.y1 * scale,
            y2: self.y2 * scale,
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
            ax2: self.ax1 * scale,
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
            x1: self.x1 + other.x1,
            y1: self.y1 + other.x2,
            vx1: self.vx1 + other.vx1,
            vy1: self.vy1 + other.vy2,
            x2: self.x2 + other.x2,
            y2: self.y2 + other.y2,
            vx2: self.vx2 + other.vx2,
            vy2: self.vy2 + other.vy2,
        }
    }
}



fn derivative(state_vector: TwoBodyState) -> DerivativeState
{

    let x1 = state_vector.x1;
    let y1 = state_vector.y1;
    let vx1 = state_vector.vx1;
    let vy1 = state_vector.vy1;
    let x2 = state_vector.x2;
    let y2 = state_vector.y2;
    let vx2 = state_vector.vx2;
    let vy2 = state_vector.vy2;
    
    let r_squared_1 = x1.powi(2) + y1.powi(2) ;
    let r_squared_2 = x2.powi(2) + y2.powi(2) ;

    let rhat_x1  = x1 / (x1.powi(2) + y1.powi(2)).sqrt();
    let rhat_y1  = y1 / (x1.powi(2) + y1.powi(2)).sqrt();
    let rhat_x2  = x2 / (x2.powi(2) + y2.powi(2)).sqrt();
    let rhat_y2  = y2 / (x2.powi(2) + y2.powi(2)).sqrt();

    let ax1 = -(GRAVITATIONAL_CONSTANT_G / (r_squared_1 * MASS_RATIO)) * rhat_x1;
    let ay1 = -(GRAVITATIONAL_CONSTANT_G / (r_squared_1 * MASS_RATIO)) * rhat_y1;
    let ax2 = -(GRAVITATIONAL_CONSTANT_G / (r_squared_2 * MASS_RATIO)) * rhat_x2;
    let ay2 = -(GRAVITATIONAL_CONSTANT_G / (r_squared_2 * MASS_RATIO)) * rhat_y2;

    DerivativeState{
        vx1: vx1,
        vy1: vx2,
        ax1: ax1,
        ay1: ay1,
        vx2: vx2,
        vy2: vy2,
        ax2: ax2,
        ay2: ay2,
    }
}





fn rk4(state_vector: TwoBodyState, dt: f64) -> [f64; 8]
{
    let k1 = dt * derivative(state_vector);
    let k2 = dt * derivative(state_vector + k1 * 0.5);
    let k3 = dt * derivative(state_vector + k2 * 0.5);
    let k4 = dt * derivative(state_vector + k3 * 1.0);


    let u = state_vector + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    u
}




fn main()
{
    let s1 = TwoBodyState{
        x1: 0.2,
        y1: 0.3,
        vx1: 0.1,
        vy1: 0.1,
        x2: 1.0,
        y2: 0.0,
        vx2: 0.2,
        vy2: 0.3,
    };


    let mut t = 0.0;
    let dt = 0.01;
    while t < 2.0 {
        let mut u = rk4(s1,dt);
        println!("{:+0.8} {:+0.8} {:+0.8}", t, u[0], u[1]);
        t = t + dt;
    }
}