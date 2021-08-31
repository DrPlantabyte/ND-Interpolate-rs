use std::fs;
use std::env;

pub fn setup() {
	let test_dir = env!("CARGO_TARGET_TMPDIR");
	fs::create_dir_all(test_dir).unwrap();
	env::set_current_dir(test_dir).unwrap();
}