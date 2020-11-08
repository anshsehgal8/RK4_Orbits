use std::ops::Mul;

static GRAVITATIONAL_CONSTANT_G: f64 = 6.67E-8;
static MASS_RATIO:f64 = 0.01;

// ============================================================================

struct TwoBodyState
{
    x1: f64,
    y1: f64,
    vx1: f64,
    vy1: f64,
    x2: f64,
    y2: f64,
    vx2: f64,
    vy2: f64,
}


// pub trait AddState {
//     fn add_pos_x(&self) -> f64;
//     fn add_pos_y(&self) -> f64;
//     fn add_vel_x(&self) -> f64;
//     fn add_vel_y(&self) -> f64;
// }

pub trait Multipy{
    fn mul(&self, scale: f64) -> f64;
}


// impl AddState for TwoBodyState{
//     fn add_pos_x(&self) -> f64 {
//         return self.x1 + self.x2 ;
//     }
//     fn add_pos_y(&self) -> f64 {
//         return self.y1 + self.y2 ;
//     }
//     fn add_vel_x(&self) -> f64 {
//         return self.vx1 + self.vx2 ; 
//     }
//     fn add_vel_y(&self) -> f64 {
//         return self.vy1 + self.vy2 ; 
//     }
// }


impl std::ops::Mul<f64> for TwoBodyState
{
    type Output = TwoBodyState;
    fn mul(&self, scale: f64) -> TwoBodyState
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




// fn main()
// {
//     let mut s1 = TwoBodyState
//     {
//         x1: 0.,
//         y1: 0.,
//         vx1: 0.,
//         vy1: 0.,
//         x2: 1.,
//         y2: 0.,
//         vx2: 0.2,
//         vy2: 0.3,

//     };

//     let mut s2 = TwoBodyState 
//     {
//         x1: 0.2,
//         y1: 0.3,
//         vx1: 0.1,
//         vy1: 0.1,
//         x2: 1.2,
//         y2: 0.2,
//         vx2: 0.1,
//         vy2: 0.4,
//     };

//     s1.mul(5.0);
//     println!("{}", s1.x1 );

// }




fn derivative_pos_x(ux: [f64;2]) -> f64
{
    let vx = ux[1];
    vx
}




fn derivative_pos_y(uy: [f64;2]) -> f64
{
    let vy = uy[1];
    vy
}




fn derivative_vel_x(ux: [f64;2]) -> f64
{
    let ax = - (GRAVITATIONAL_CONSTANT_G * (1. + MASS_RATIO)) / (ux[0].powi(2));
    ax
}




fn derivative_vel_y(uy: [f64;2]) -> f64
{
    let ay = - (GRAVITATIONAL_CONSTANT_G * (1. + MASS_RATIO)) / (uy[0].powi(2));
    ay
}




fn rk4(ux: [f64; 2], uy: [f64;2], dt: f64) -> [f64; 4]
{
    let k1x = dt * derivative_pos_x(ux);
    let k1y = dt * derivative_pos_y(uy);
    let l1x = dt * derivative_vel_x(ux);
    let l1y = dt * derivative_vel_y(uy);
    let k2x = dt * derivative_pos_x([ux[0] + k1x * 0.5, ux[1] + l1x * 0.5]);
    let k2y = dt * derivative_pos_y([uy[0] + k1y * 0.5, uy[1] + l1y * 0.5]);
    let l2x = dt * derivative_vel_x([ux[0] + k1x * 0.5, ux[1] + l1x * 0.5]);
    let l2y = dt * derivative_vel_y([uy[0] + k1y * 0.5, uy[1] + l1y * 0.5]);
    let k3x = dt * derivative_pos_x([ux[0] + k2x * 0.5, ux[1] + l2x * 0.5]);
    let k3y = dt * derivative_pos_y([uy[0] + k2y * 0.5, uy[1] + l2y * 0.5]);
    let l3x = dt * derivative_vel_x([ux[0] + k2x * 0.5, ux[1] + l2x * 0.5]);
    let l3y = dt * derivative_vel_y([uy[0] + k2y * 0.5, uy[1] + l2y * 0.5]);
    let k4x = dt * derivative_pos_x([ux[0] + k3x * 1.0, ux[1] + l3x * 1.0]);
    let k4y = dt * derivative_pos_y([uy[0] + k3y * 1.0, uy[1] + l3y * 1.0]);
    let l4x = dt * derivative_vel_x([ux[0] + k3x * 1.0, ux[1] + l3x * 1.0]);
    let l4y = dt * derivative_vel_y([uy[0] + k3y * 1.0, uy[1] + l3y * 1.0]);
    let x_new = ux[0] + (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
    let y_new = uy[0] + (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
    let vx_new = ux[1] + (l1x + 2.0 * l2x + 2.0 * l3x + l4x) / 6.0;
    let vy_new = uy[1] + (l1y + 2.0 * l2y + 2.0 * l3y + l4y) / 6.0;
    let u = [x_new, y_new, vx_new, vy_new];
    u
}




fn main()
{
    let s1 = TwoBodyState{
        x: 0.0,
        y: 0.0,
        vx: 0.0,
        vy: 0.0,
    };

    let s2 = TwoBodyState{
        x: 1.0,
        y: 0.0,
        vx: 0.2,
        vy: 0.3,
    };
    let mut u = [1.0,0.0,0.5,0.1]; // [x, y, xdot, ydot]
    let mut t = 0.0;
    let dt = 0.01;
    while t < 2.0 {
        u = rk4([u[0],u[2]], [u[1],u[3]], dt);
        println!("{:+0.8} {:+0.8} {:+0.8}", t, u[0], u[1]);
        t = t + dt;
    }
}