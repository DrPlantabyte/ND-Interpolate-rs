#![crate_name = "nd_interpolate"]
#![warn(missing_docs)]
//! This crate provides linear and cubic interpolators for up to 4 dimensions

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
	/// four-point X,Y bi-linear interpolation in a grid, using the assumption that x and y are
	/// between the provided reference grid points
	/// # Arguments
	/// * `x` - x position to interpolate the value at (will be normalized by the function
	/// `X' = X - floor(X)` )
	/// * `y` - y position to interpolate the value at(will be normalized by the function
	/// `Y' = Y - floor(Y)` )
	/// * `local_2x2` - the grid points immediately around point (X,Y)
	/// # Returns
	/// The interpolated value at (X,Y)
	///
	pub fn linear_2D_grid(x: f64, y: f64, local_2x2: &[[f64;2];2]) -> f64 {
		let ix = f64::floor(x) as usize;
		let xx = x - ix as f64;
		let iy = f64::floor(y) as usize;
		let yy = y - iy as f64;
		let y0 = linear_1D(xx, 0f64, local_2x2[0][0], 1f64, local_2x2[1][0]);
		let yp1 = linear_1D(xx, 0f64, local_2x2[0][1], 1f64, local_2x2[1][1]);
		return linear_1D(yy, 0f64, y0, 1f64, yp1);
	}


	/// eight-point X,Y,Z tri-linear interpolation in a grid, using the assumption that (X,Y,Z) is
	/// between the provided reference grid points
	/// # Arguments
	/// * `x` - x position to interpolate the value at (will be normalized by the function
	/// `X' = X - floor(X)` )
	/// * `y` - y position to interpolate the value at(will be normalized by the function
	/// `Y' = Y - floor(Y)` )
	/// * `z` - z position to interpolate the value at(will be normalized by the function
	/// `Z' = Z - floor(Z)` )
	/// * `local_2x2x2` - the grid points immediately around point (X,Y,Z)
	/// # Returns
	/// The interpolated value at (X,Y,Z)
	///
	pub fn linear_3D_grid(x: f64, y: f64, z: f64, local_2x2x2: &[[[f64;2];2];2]) -> f64 {
		let iz = f64::floor(z) as usize;
		let zz = z - iz as f64;
		let zy0 = linear_2D_grid(x, y, &local_2x2x2[0]);
		let zy1 = linear_2D_grid(x, y, &local_2x2x2[1]);
		return linear_1D(zz, 0f64, zy0, 1f64, zy1);
	}
	/// sixteen-point X,Y,Z,W tetra-linear interpolation in a grid, using the assumption that
	/// (X,Y,Z,W) is between the provided reference grid points
	/// # Arguments
	/// * `x` - x position to interpolate the value at (will be normalized by the function
	/// `X' = X - floor(X)` )
	/// * `y` - y position to interpolate the value at(will be normalized by the function
	/// `Y' = Y - floor(Y)` )
	/// * `z` - z position to interpolate the value at(will be normalized by the function
	/// `Z' = Z - floor(Z)` )
	/// * `w` - w position to interpolate the value at(will be normalized by the function
	/// `W' = W - floor(W)` )
	/// * `local_2x2x2x2` - the grid points immediately around point (X,Y,Z,W)
	/// # Returns
	/// The interpolated value at (X,Y,Z,W)
	///
	pub fn linear_4D_grid(
		x: f64, y: f64, z: f64, w: f64, local_2x2x2x2: &[[[[f64;2];2];2];2]
	) -> f64 {
		let iw = f64::floor(w) as usize;
		let ww = w - iw as f64;
		let wy0 = linear_3D_grid(x, y, z,&local_2x2x2x2[0]);
		let wy1 = linear_3D_grid(x, y, z,&local_2x2x2x2[1]);
		return linear_1D(ww, 0f64, wy0, 1f64, wy1);
	}


	/// Four-point cubic interpolation of point (X) in a grid of Y-values. This function assumes
	/// that X is in the middle of the provided grid array (that is, X lies between `local_4[1]`
	/// and `local_4[2]`). X will be normalized with `X - floor(X)` where 0 equals `local_4[1]`
	/// and 0.999... approaches `local_4[2]`, interpolating with a cubic spline.
	/// # Arguments
	/// * `x` - fractional grid coordinate to a interpolate at
	/// * `local_4` - Reference points for cubic interpolation surrounding the coordinate of
	/// interest
	/// # Returns
	/// Returns the cubic-spline interpolated value of Y at X location between `local_4[1]` and
	/// `local_4[2]`
	pub fn cubic_1D_grid(x:f64, local_4: &[f64;4]) -> f64{
		// from https://github.com/DrPlantabyte/Cyanos-Noise-Library/blob/master/cchall.noise/src/cchall/noise/math/CubicInterpolator.java
		let yn2 = local_4[0];
		let yn1 = local_4[1];
		let yp1 = local_4[2];
		let yp2 = local_4[3];
		let w = x - f64::floor(x);
		let O1 = -0.5 * yn2;
		let O2 = 0.5 * yp2;
		let O3 = w * w;
		let A = O1 + 1.5 * yn1 - 1.5 * yp1 + O2;
		let B = yn2 - 2.5 * yn1 + 2.0 * yp1 - O2;
		let C = O1 + 0.5 * yp1;
		let D = yn1;
		return A * O3 * w + B * O3 + C * w + D;
	}
	/// 16-point bi-cubic interpolation of a coordinate in a grid of values. This function assumes
	/// that `coord` is in the middle of the provided grid array, such that `floor(coord)` represents
	/// index\[1,1\] in the grid and `floor(coord)+1` represents index\[2,2\] in the grid
	/// # Arguments
	/// * `coord` - coordinate position within the grid (each dimension will be normalized with
	/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
	/// larger grid)
	/// * `local_4x4` - Reference points for cubic interpolation surrounding the coordinate of
	/// interest
	/// # Returns
	/// Returns the cubic-spline interpolated value at the provided coordinate
	pub fn cubic_2D_grid(coord: [f64;2], local_4x4: &[[f64;4];4]) -> f64{
		let mut local4 = [0f64;4];
		for n in 0..4{
			local4[n] = cubic_1D_grid(coord[1], &local_4x4[n]);
		}
		return cubic_1D_grid(coord[0], &local4);
	}
	/// 64-point tri-cubic interpolation of a coordinate in a grid of values. This function assumes
	/// that `coord` is in the middle of the provided grid array, such that `floor(coord)` represents
	/// index\[1,1,...\] in the grid and `floor(coord)+1` represents index\[2,2,...\] in the grid
	/// # Arguments
	/// * `coord` - coordinate position within the grid (each dimension will be normalized with
	/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
	/// larger grid)
	/// * `local_4x4x4` - Reference points for cubic interpolation surrounding the coordinate of
	/// interest
	/// # Returns
	/// Returns the cubic-spline interpolated value at the provided coordinate
	pub fn cubic_3D_grid(coord: [f64;3], local_4x4x4: &[[[f64;4];4];4]) -> f64{
		let mut local4 = [0f64;4];
		let subcoord = [coord[1], coord[2]];
		for n in 0..4{
			local4[n] = cubic_2D_grid(subcoord, &local_4x4x4[n]);
		}
		return cubic_1D_grid(coord[0], &local4);
	}
	/// 256-point tetra-cubic interpolation of a coordinate in a grid of values. This function
	/// assumes that `coord` is in the middle of the provided grid array, such that
	/// `floor(coord)` represents index\[1,1,...\] in the grid and `floor(coord)+1` represents
	/// index\[2,2,...\] in the grid
	/// # Arguments
	/// * `coord` - coordinate position within the grid (each dimension will be normalized with
	/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
	/// larger grid)
	/// * `local_4x4x4x4` - Reference points for cubic interpolation surrounding the coordinate of
	/// interest
	/// # Returns
	/// Returns the cubic-spline interpolated value at the provided coordinate
	pub fn cubic_4D_grid(coord: [f64;4], local_4x4x4x4: &[[[[f64;4];4];4];4]) -> f64{
		let mut local4 = [0f64;4];
		let subcoord = [coord[1], coord[2], coord[3]];
		for n in 0..4{
			local4[n] = cubic_3D_grid(subcoord, &local_4x4x4x4[n]);
		}
		return cubic_1D_grid(coord[1], &local4);
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
