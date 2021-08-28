use nd_interpolate;
mod test_manager;

extern crate image;
use image::*;
use rand::prelude::*;
use nd_interpolate::f64_data::linear_2D_grid;


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
	for px in 0..size{
		for py in 0..size{
			let x = rand_grid.len() as f64 * px as f64 / size as f64;
			let y = rand_grid[0].len() as f64 * py as f64 / size as f64;
			let ix = f64::floor(x) as usize % rand_grid.len();
			let iy = f64::floor(y) as usize % rand_grid[ix].len();
			let ixp1 = (ix+1) % rand_grid.len();
			let iyp1 = (iy+1) % rand_grid[ixp1].len();
			let v = ((linear_2D_grid(x, y, &[[ rand_grid[ix][iy],rand_grid[ix][iyp1] ],[
				rand_grid[ixp1][iy],rand_grid[ixp1][iyp1] ]])) * 255f64) as u8;
			let pixel = Rgb([v,v,v]);
			img.put_pixel(px, py, pixel);
		}
	}
	img.save("linear_image_interpolation.png").unwrap();
	println!("TODO: integration test")
}
