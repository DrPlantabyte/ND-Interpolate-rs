use nd_interpolate;
mod test_manager;

extern crate image;
use image::*;
use rand::prelude::*;
use nd_interpolate::f64_data::*;
use std::f64::consts::PI;

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
