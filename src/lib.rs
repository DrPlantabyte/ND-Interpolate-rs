#![crate_name = "nd_interpolate"]
#![warn(missing_docs)]
#![allow(non_snake_case)]
//! This crate provides linear and cubic interpolators for up to 4 dimensions

use num_traits::Float;
extern crate num_traits;

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
pub fn linear_1D<T: Float>(x: T, x0: T, y0: T, xp1: T, yp1: T) -> T {
	let w = (x-x0) / (xp1-x0);
	return yp1 * w + y0 * (T::one() - w);
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
pub fn linear_2D_grid<T: Float>(x: T, y: T, local_2x2: &[[T;2];2]) -> T {
	let ix = T::floor(x);
	let xx = x - ix;
	let iy = T::floor(y);
	let yy = y - iy;
	let xy0 = linear_1D(yy, T::zero(), local_2x2[0][0], T::one(), local_2x2[0][1]);
	let xy1 = linear_1D(yy, T::zero(), local_2x2[1][0], T::one(), local_2x2[1][1]);
	let yx0 = linear_1D(xx, T::zero(), local_2x2[0][0], T::one(), local_2x2[1][0]);
	let yx1 = linear_1D(xx, T::zero(), local_2x2[0][1], T::one(), local_2x2[1][1]);
	return T::from(0.5).unwrap() * (linear_1D(xx, T::zero(), xy0, T::one(), xy1) + linear_1D(yy, T::zero(), yx0, T::one(), yx1));
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
pub fn linear_3D_grid<T: Float>(x: T, y: T, z: T, local_2x2x2: &[[[T;2];2];2]) -> T {
	let ix = T::floor(x);
	let xx = x - ix;
	let xy0 = linear_2D_grid(y, z, &local_2x2x2[0]);
	let xy1 = linear_2D_grid(y, z, &local_2x2x2[1]);
	return linear_1D(xx, T::zero(), xy0, T::one(), xy1);
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
pub fn linear_4D_grid<T: Float>(
	x: T, y: T, z: T, w: T, local_2x2x2x2: &[[[[T;2];2];2];2]
) -> T {
	let ix = T::floor(x);
	let xx = x - ix;
	let xy0 = linear_3D_grid(y, z, w, &local_2x2x2x2[0]);
	let xy1 = linear_3D_grid(y, z, w, &local_2x2x2x2[1]);
	return linear_1D(xx, T::zero(), xy0, T::one(), xy1);
}

/// 32-point hyperdimensional-linear interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[0,0,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[1,1,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_2x_grid` - Reference points for linear interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the linear interpolated value at the provided coordinate
pub fn linear_5D_grid<T: Float>(coord: [T;5], local_2x_grid: &[[[[[T;2];2];2];2];2]) -> T {
	let x = coord[0];
	let ix = T::floor(x);
	let xx = x - ix;
	let mut local2 = [T::zero();2];
	for n in 0..2{
		local2[n] = linear_4D_grid(coord[1], coord[2], coord[3], coord[4], &local_2x_grid[n]);
	}
	return linear_1D(xx, T::zero(), local2[0], T::one(), local2[1]);
}
/// 64-point hyperdimensional-linear interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[0,0,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[1,1,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_2x_grid` - Reference points for linear interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the linear interpolated value at the provided coordinate
pub fn linear_6D_grid<T: Float>(coord: [T;6], local_2x_grid: &[[[[[[T;2];2];2];2];2];2]) -> T {
	let x = coord[0];
	let ix = T::floor(x);
	let xx = x - ix;
	let mut local2 = [T::zero();2];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5]];
	for n in 0..2{
		local2[n] = linear_5D_grid(subcoord, &local_2x_grid[n]);
	}
	return linear_1D(xx, T::zero(), local2[0], T::one(), local2[1]);
}
/// 128-point hyperdimensional-linear interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[0,0,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[1,1,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_2x_grid` - Reference points for linear interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the linear interpolated value at the provided coordinate
pub fn linear_7D_grid<T: Float>(coord: [T;7], local_2x_grid: &[[[[[[[T;2];2];2];2];2];2];2]) -> T {
	let x = coord[0];
	let ix = T::floor(x);
	let xx = x - ix;
	let mut local2 = [T::zero();2];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6]];
	for n in 0..2{
		local2[n] = linear_6D_grid(subcoord, &local_2x_grid[n]);
	}
	return linear_1D(xx, T::zero(), local2[0], T::one(), local2[1]);
}
/// 256-point hyperdimensional-linear interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[0,0,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[1,1,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_2x_grid` - Reference points for linear interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the linear interpolated value at the provided coordinate
pub fn linear_8D_grid<T: Float>(coord: [T;8], local_2x_grid: &[[[[[[[[T;2];2];2];2];2];2];2];2]) -> T {
	let x = coord[0];
	let ix = T::floor(x);
	let xx = x - ix;
	let mut local2 = [T::zero();2];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6], coord[7]];
	for n in 0..2{
		local2[n] = linear_7D_grid(subcoord, &local_2x_grid[n]);
	}
	return linear_1D(xx, T::zero(), local2[0], T::one(), local2[1]);
}
/// 512-point hyperdimensional-linear interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[0,0,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[1,1,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_2x_grid` - Reference points for linear interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the linear interpolated value at the provided coordinate
pub fn linear_9D_grid<T: Float>(coord: [T;9], local_2x_grid: &[[[[[[[[[T;2];2];2];2];2];2];2];2];2]) -> T {
	let x = coord[0];
	let ix = T::floor(x);
	let xx = x - ix;
	let mut local2 = [T::zero();2];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6], coord[7], coord[8]];
	for n in 0..2{
		local2[n] = linear_8D_grid(subcoord, &local_2x_grid[n]);
	}
	return linear_1D(xx, T::zero(), local2[0], T::one(), local2[1]);
}
/// 1024-point hyperdimensional-linear interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[0,0,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[1,1,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_2x_grid` - Reference points for linear interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the linear interpolated value at the provided coordinate
pub fn linear_10D_grid<T: Float>(coord: [T;10], local_2x_grid: &[[[[[[[[[[T;2];2];2];2];2];2];2];2];2];2]) -> T {
	let x = coord[0];
	let ix = T::floor(x);
	let xx = x - ix;
	let mut local2 = [T::zero();2];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6], coord[7], coord[8], coord[9]];
	for n in 0..2{
		local2[n] = linear_9D_grid(subcoord, &local_2x_grid[n]);
	}
	return linear_1D(xx, T::zero(), local2[0], T::one(), local2[1]);
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
pub fn cubic_1D_grid<T: Float>(x: T, local_4: &[T;4]) -> T {
	// from https://github.com/DrPlantabyte/Cyanos-Noise-Library/blob/master/cchall.noise/src/cchall/noise/math/CubicInterpolator.java
	let yn2 = local_4[0];
	let yn1 = local_4[1];
	let yp1 = local_4[2];
	let yp2 = local_4[3];
	let w = x - T::floor(x);
	let O1 = T::from(-0.5).unwrap() * yn2;
	let O2 =  T::from(0.5).unwrap() * yp2;
	let O3 = w * w;
	let A = O1 +  T::from(1.5).unwrap() * yn1 -  T::from(1.5).unwrap() * yp1 + O2;
	let B = yn2 -  T::from(2.5).unwrap() * yn1 +  T::from(2.0).unwrap() * yp1 - O2;
	let C = O1 +  T::from(0.5).unwrap() * yp1;
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
pub fn cubic_2D_grid<T: Float>(coord: [T;2], local_4x4: &[[T;4];4]) -> T {
	let mut local4 = [T::zero();4];
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
pub fn cubic_3D_grid<T: Float>(coord: [T;3], local_4x4x4: &[[[T;4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2]];
	for n in 0..4{
		local4[n] = cubic_2D_grid(subcoord, &local_4x4x4[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 256-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_4D_grid<T: Float>(coord: [T;4], local_4x_grid: &[[[[T;4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3]];
	for n in 0..4{
		local4[n] = cubic_3D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 1024-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_5D_grid<T: Float>(coord: [T;5], local_4x_grid: &[[[[[T;4];4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3], coord[4]];
	for n in 0..4{
		local4[n] = cubic_4D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 4096-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_6D_grid<T: Float>(coord: [T;6], local_4x_grid: &[[[[[[T;4];4];4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5]];
	for n in 0..4{
		local4[n] = cubic_5D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 16,384-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_7D_grid<T: Float>(coord: [T;7], local_4x_grid: &[[[[[[[T;4];4];4];4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6]];
	for n in 0..4{
		local4[n] = cubic_6D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 65,536-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_8D_grid<T: Float>(coord: [T;8], local_4x_grid: &[[[[[[[[T;4];4];4];4];4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6], coord[7]];
	for n in 0..4{
		local4[n] = cubic_7D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 262,144-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # WARNING: Large stack memory usage!
/// You will need to increase the stack size limit to more than 2 MB to use this function:
/// ```
/// use std::thread;
/// use nd_interpolate::cubic_9D_grid;
/// let builder = thread::Builder::new()
///   .name("big-stack-thread".into())
///   .stack_size(8 * 1024 * 1024); // 8MB of stack space
/// let handler = builder.spawn(|| {
///   let mut grid = [[[[[[[[[0f64;4];4];4];4];4];4];4];4];4];
///   let target_coord = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9];
///   let interpolated_value = cubic_9D_grid(target_coord, &grid);
/// }).unwrap();
/// handler.join().unwrap(); 
/// ```
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_9D_grid<T: Float>(coord: [T;9], local_4x_grid: &[[[[[[[[[T;4];4];4];4];4];4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6], coord[7], coord[8]];
	for n in 0..4{
		local4[n] = cubic_8D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}

/// 1,048,576-point hyperdimensional-cubic interpolation of a coordinate in a 
/// grid of values. This function assumes that `coord` is in the middle of the 
/// provided grid array, such that `floor(coord)` represents index\[1,1,...\] 
/// in the grid and `floor(coord)+1` represents
/// index\[2,2,...\] in the grid
/// # WARNING: Large stack memory usage!
/// You will need to increase the stack size limit to more than 8 MB to use this function:
/// ```
/// use std::thread;
/// use nd_interpolate::cubic_10D_grid;
/// let builder = thread::Builder::new()
///   .name("big-stack-thread".into())
///   .stack_size(16 * 1024 * 1024); // 16MB of stack space
/// let handler = builder.spawn(|| {
///   let mut grid = [[[[[[[[[[0f64;4];4];4];4];4];4];4];4];4];4];
///   let target_coord = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9];
///   let interpolated_value = cubic_10D_grid(target_coord, &grid);
/// }).unwrap();
/// handler.join().unwrap(); 
/// ```
/// # Arguments
/// * `coord` - coordinate position within the grid (each dimension will be normalized with
/// `X = X - floor(X)` so that you don't need to correct the position when subsampling from a
/// larger grid)
/// * `local_4x_grid` - Reference points for cubic interpolation surrounding the coordinate of
/// interest
/// # Returns
/// Returns the cubic-spline interpolated value at the provided coordinate
pub fn cubic_10D_grid<T: Float>(coord: [T;10], local_4x_grid: &[[[[[[[[[[T;4];4];4];4];4];4];4];4];4];4]) -> T {
	let mut local4 = [T::zero();4];
	let subcoord = [coord[1], coord[2], coord[3], coord[4], coord[5], coord[6], coord[7], coord[8], coord[9]];
	for n in 0..4{
		local4[n] = cubic_9D_grid(subcoord, &local_4x_grid[n]);
	}
	return cubic_1D_grid(coord[0], &local4);
}




#[cfg(test)]
mod tests {
	use super::*;
	#[test]
	fn linear_test_f64() {
		assert_eq!(linear_1D(0.25, 0., 1., 1., 3.), 1.5);
		// a=numpy.asarray([[0,1],[2,3]])
		// interpolator = scipy.interpolate.LinearNDInterpolator(
		//   numpy.asarray([[0,0],[1,0],[0,1],[1,1]]), 
		//   numpy.asarray([a[0][0], a[1][0], a[0][1], a[1][1]]))
		// interpolator(0.25, 0.75) => 1.25
		assert_eq!(linear_2D_grid(0.25, 0.75, &[[0.,1.],[2.,3.]]), 1.25);
	}
	#[test]
	fn cubic_test_f64() {
		let a = 0.019; let b = -0.15; let c = 1.0; let d = -13.0;
		let x = 1.64; let y = a*x*x*x + b*x*x + c*x + d;
		let local_4 = [d, a + b + c + d, a*8.0 + b*4.0 + c*2.0 + d, a*27.0 + b*9.0 + c*3.0 + d];
		let yi = cubic_1D_grid(x, &local_4);
		let percent_delta = 100. * ((yi-y)/y).abs();
		println!("{} ?= {} (% delta = {}%)", yi, y, percent_delta);
		assert!(percent_delta < 0.1);
	}
}



