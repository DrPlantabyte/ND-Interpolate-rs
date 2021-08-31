# ND-Interpolate-rs
1-10 dimensional interpolation implemented in Rust. It supports both f64 and f32 data types for interpolation, and has both linear and cubic spline interpolation for up to 10 dimensions.

This library is particularly useful for Perlin Noise applications.

# Example
Here's an example of using 3D cubic interpolation to make a marble-pattern texture for a sphere:
```rust
use nd_interpolate::f64_data::*;
use image::*; // image = "0.23.14"
use rand::prelude::*; // rand = "0.8.4"
use std::f64::consts::PI;

fn sphere_marble_texture_example() {
	println!("Making marbled sphere texture...");
	let imgsize = 512; let w = 2*imgsize; let h = imgsize;
	let gridsize = 16; let gridcenter = gridsize / 2;
	let radius = (gridsize / 2 - 3) as f64;
	let mut img: RgbImage = ImageBuffer::new(w, h);
	let mut rand_grid = vec![vec![vec![0f64;gridsize];gridsize];gridsize];
	for x in 0..rand_grid.len(){ for y in 0..rand_grid[x].len(){ for z in 0..rand_grid[x][y].len() {
		let r: f64 = random();
		rand_grid[x][y][z] = r;
	} } }
	for py in 0..h { for px in 0..w {
		let lon = (px as f64 / w as f64) * (2.0*PI) - PI;
		let lat = ((h-py) as f64 / h as f64) * PI - 0.5*PI;
		let coord = [
			radius * lon.sin() * lat.cos(), // X
			radius * lat.sin(), // Y
			radius * lon.cos() * lat.cos(), // Z
		];
		let grid_coord: [usize; 3] = [
			(coord[0].floor() as i32 + gridcenter as i32) as usize,
			(coord[1].floor() as i32 + gridcenter as i32) as usize,
			(coord[2].floor() as i32 + gridcenter as i32) as usize
		];
		let mut local_4x4x4 = [[[0f64;4];4];4];
		for ix in 0..4{for iy in 0..4{for iz in 0..4{
			let gx = grid_coord[0];
			let gy = grid_coord[1];
			let gz = grid_coord[2];
			local_4x4x4[ix][iy][iz] = rand_grid[ix+gx-1][iy+gy-1][iz+gz-1];
		}}}
		let v = ((cubic_3D_grid(coord, &local_4x4x4)) * 255f64) as u8;
		let pixel = Rgb([v,v,v]);
		img.put_pixel(px, py, pixel);
	} }
	img.save("marbled_sphere_texture.png").unwrap();
	println!("...Done!");
}
```
The above code produced the following image:
![3D_cubic_image_interpolation](https://user-images.githubusercontent.com/1922739/131440111-1e3bebe2-23d7-48e4-9eb7-6bdade7d04bd.png)



# Notes
This Rust library is based on [my Java noise library interpolator](https://github.com/DrPlantabyte/Cyanos-Noise-Library/blob/master/cchall.noise/src/cchall/noise/math/CubicInterpolator.java)
