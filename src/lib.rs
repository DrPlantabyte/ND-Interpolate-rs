#![crate_name = "nd_interpolate"]
/// interpolation functions using 64-bit (aka double-precision) floating point numbers
pub mod f64_data {

	#![allow(non_snake_case)]
	/// two-point x,y linear interpolation
	/// # Arguments
	/// * `x` - x position to interpolate the value at
	/// * `x0` - reference point before x
	/// * `y0` - Y value at x0
	/// * `xp1` - reference point after x
	/// * `yp1` - Y value at xp1
	/// # Returns
	/// The interpolated Y value at x
	///
	pub fn linear_1D(x: f64, x0: f64, y0: f64, xp1: f64, yp1: f64) -> f64 {
		let w = (x-x0) / (xp1-x0);
		return yp1 * w + y0 * (1. - w);
	}

	pub fn linear_2D_grid(x: f64, y: f64, local_2x2: &[[f64;2];2]) -> f64 {
		let ix = f64::floor(x) as usize;
		let xx = x - ix as f64;
		let iy = f64::floor(y) as usize;
		let yy = y - iy as f64;
		let y0 = linear_1D(xx, 0f64, local_2x2[0][0], 1f64, local_2x2[1][0]);
		let yp1 = linear_1D(xx, 0f64, local_2x2[0][1], 1f64, local_2x2[1][1]);
		return linear_1D(yy, 0f64, y0, 1f64, yp1);
	}


	pub fn linear_3D_grid(x: f64, y: f64, z: f64, local_2x2x2: &[[[f64;2];2];2]) -> f64 {
		let ix = f64::floor(x) as usize;
		let xx = x - ix as f64;
		let iy = f64::floor(y) as usize;
		let yy = y - iy as f64;
		let iz = f64::floor(z) as usize;
		let zz = z - iz as f64;
		let zy0 = linear_2D_grid(xx, yy, &local_2x2x2[0]);
		let zy1 = linear_2D_grid(xx, yy, &local_2x2x2[1]);
		return linear_1D(zz, 0f64, zy0, 1f64, zy1);
	}

	pub fn linear_4D_grid(x: f64, y: f64, z: f64, w: f64, local_2x2x2x2: &[[[[f64;2];2];2];2]) ->
																							   f64 {
		let ix = f64::floor(x) as usize;
		let xx = x - ix as f64;
		let iy = f64::floor(y) as usize;
		let yy = y - iy as f64;
		let iz = f64::floor(z) as usize;
		let zz = z - iz as f64;
		let iw = f64::floor(w) as usize;
		let ww = w - iw as f64;
		let wy0 = linear_3D_grid(xx, yy, zz,&local_2x2x2x2[0]);
		let wy1 = linear_3D_grid(xx, yy, zz,&local_2x2x2x2[1]);
		return linear_1D(ww, 0f64, wy0, 1f64, wy1);
	}

	#[cfg(test)]
	mod tests {
		use crate::f64_data;
		#[test]
		fn it_works() {
			assert_eq!(f64_data::linear_1D(0.25, 0., 1., 1., 3.), 1.5);
		}
	}
}
