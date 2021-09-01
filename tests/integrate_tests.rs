use nd_interpolate;
mod test_manager;

extern crate image;
use image::*;
use rand::prelude::*;
use nd_interpolate::*;
use std::f64::consts::PI;
use std::thread;

#[test]
fn linear_image_interpolation() {
	test_manager::setup();

	let size = 512;
	let mut img: RgbImage = ImageBuffer::new(size, size);
	let mut rand_grid = [[0f64; 3]; 3];
	for x in 0..rand_grid.len(){
		for y in 0..rand_grid[x].len(){
			let r: f64 = random();
			rand_grid[x][y] = r;
		}
	}
	for py in 0..size{
		for px in 0..size{
			let x = rand_grid.len() as f64 * px as f64 / size as f64;
			let y = rand_grid[0].len() as f64 * py as f64 / size as f64;
			let ix = f64::floor(x) as usize % rand_grid.len();
			let iy = f64::floor(y) as usize % rand_grid[ix].len();
			let ixp1 = (ix+1) % rand_grid.len();
			let iyp1 = (iy+1) % rand_grid[ixp1].len();
			let v = ((linear_2D_grid(x, y, &[[ rand_grid[ix][iy],rand_grid[ix][iyp1] ],[
				rand_grid[ixp1][iy],rand_grid[ixp1][iyp1] ]])) * 255f64) as u8;
			let v = v.max(0).min(255);
			let pixel = Rgb([v,v,v]);
			img.put_pixel(px, py, pixel);
		}
	}
	img.save("linear_image_interpolation.png").unwrap();
}

#[test]
fn cubic_image_interpolation() {
	test_manager::setup();
	println!("tri-cubic interpolation test");
	let imgsize = 512;
	let w = 2*imgsize;
	let h = imgsize;
	let gridsize = 16;
	let gridcenter = gridsize / 2;
	let radius = (gridsize / 2 - 3) as f64;
	let mut img: RgbImage = ImageBuffer::new(w, h);
	let mut rand_grid = vec![vec![vec![0f64;gridsize];gridsize];gridsize];
	for x in 0..rand_grid.len(){
		for y in 0..rand_grid[x].len(){
			for z in 0..rand_grid[x][y].len() {
				let r: f64 = random();
				rand_grid[x][y][z] = r;
			}
		}
	}
	for py in 0..h {
		for px in 0..w {
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
			let mut local = [[[0f64;4];4];4];
			for ix in 0..4{for iy in 0..4{for iz in 0..4{
				let gx = grid_coord[0];
				let gy = grid_coord[1];
				let gz = grid_coord[2];
				local[ix][iy][iz] = rand_grid[ix+gx-1][iy+gy-1][iz+gz-1];
			}}}
			let v = ((cubic_3D_grid(coord, &local)) * 255f64) as u8;
			let pixel = Rgb([v,v,v]);
			img.put_pixel(px, py, pixel);
		}
	}
	img.save("3D_cubic_image_interpolation.png").unwrap();
}


#[test]
#[allow(non_snake_case)]
fn test_10d_linear_interpolation(){
	let target_coord = [0.9, 0.01, 0.8, 0.1, 0.7, 0.2, 0.6, 0.3, 0.5, 0.4];
	let slopes  = [-6.4, -8.5, -1.5, 6.4, 4.1, -1.9, 6.3, 9.6, 4.7, 5.2];
	let offsets = [8.9, -0.7, 4.1, -6.0, -3.1, -2.1, -7.9, 1.4, 5.0, -4.8];
	let target_value = linfunc(&target_coord, &slopes, &offsets);
	let mut grid = [[[[[[[[[[0f64;2];2];2];2];2];2];2];2];2];2];
	for d0 in 0..2{ for d1 in 0..2{ for d2 in 0..2{ for d3 in 0..2{ for d4 in 0..2{ for d5 in 0..2{ for d6 in 0..2{ for d7 in 0..2{ for d8 in 0..2{ for d9 in 0..2{
		// OMFG! 10D, you crazy
		let v = linfunc(&[d0 as f64, d1 as f64, d2 as f64, d3 as f64, d4 as f64, d5 as f64, d6 as f64, d7 as f64, d8 as f64, d9 as f64], &slopes, &offsets);
		grid[d0][d1][d2][d3][d4][d5][d6][d7][d8][d9] = v;
	}}}}}}}}}}
	let ival = linear_10D_grid(target_coord, &grid);
	let percent_delta = 100.0 * f64::abs((target_value - ival)/target_value);
	println!("10D linear target value = {}, interpolated = {}, % delta = {}%", target_value, ival, percent_delta);
	assert!(percent_delta < 0.1);
}
fn linfunc(coord: &[f64;10], slope: &[f64;10], offset: &[f64;10]) -> f64{
	let mut sum = 0.;
	for n in 0..10{
		sum += coord[n] * slope[n] + offset[n];
	}
	return sum;
}


#[test]
#[allow(non_snake_case)]
fn test_10d_cubic_interpolation(){
	// 10D is insane, not even going to try to visualize it
	// instead, we'll use fixed values
	// if every point in the grid is a function of its distance from the origin,
	// then the interpolated value should be roughly the same as it's distance function

	// NOTE: standard stack size limit is not large enough to hold a size 4 10D array :-(
	let builder = thread::Builder::new()
	.name("big-stack-thread".into())
	.stack_size(16 * 1024 * 1024); // 16MB of stack space
	let handler = builder.spawn(|| {
		// stack-intensive operations
		let A = 0.011; let B = -0.13; let C = 1.5; let D = -17.0;
		let target_coord = [1.9, 1.01, 1.8, 1.1, 1.7, 1.2, 1.6, 1.3, 1.5, 1.4];
		let target_value = cubic(dist_10d(&target_coord), A, B, C, D);
		let mut grid = [[[[[[[[[[0f64;4];4];4];4];4];4];4];4];4];4];
		for d0 in 0..4{ for d1 in 0..4{ for d2 in 0..4{ for d3 in 0..4{ for d4 in 0..4{ for d5 in 0..4{ for d6 in 0..4{ for d7 in 0..4{ for d8 in 0..4{ for d9 in 0..4{
			// OMFG! 10D, you crazy
			let d = dist_10d(&[d0 as f64, d1 as f64, d2 as f64, d3 as f64, d4 as f64, d5 as f64, d6 as f64, d7 as f64, d8 as f64, d9 as f64]);
			grid[d0][d1][d2][d3][d4][d5][d6][d7][d8][d9] = cubic(d, A, B, C, D);
		}}}}}}}}}}
		let ival = cubic_10D_grid(target_coord, &grid);
		let percent_delta = 100.0 * f64::abs((target_value - ival)/target_value);
		println!("10D cubic target value = {}, interpolated = {}, % delta = {}%", target_value, ival, percent_delta);
		assert!(percent_delta < 0.1);
	}).unwrap();
	handler.join().unwrap(); 
}
fn dist_10d(coord: &[f64;10]) -> f64{
	let mut sqsum = 0f64;
	for n in 0..coord.len(){
		sqsum += coord[n] * coord[n];
	}
	return f64::sqrt(sqsum);
}
#[allow(non_snake_case)]
fn cubic(x: f64, A: f64, B: f64, C: f64, D: f64) -> f64{
	return A*x*x*x+B*x*x+C*x+D;
}
