#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::rc::Rc;

use log::error;
use material::{Lambertian, Metal};
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::WindowBuilder;
use winit_input_helper::WinitInputHelper;

// I'm not sure that I'm doing this correctly.
mod util;
use util::{Color, Hittable, HittableList, Point3, Ray, Sphere, Vec3};

mod camera;
use camera::Camera;

mod material;

use rand::Rng;

fn clamp(x: f32, min: f32, max: f32) -> f32 {
    if x < min {
        return min;
    } else if x > max {
        return max;
    }

    x
}

fn calculate_color(color: Vec3, samples_per_pixel: i32) -> Vec3 {
    let mut r = color.x;
    let mut g = color.y;
    let mut b = color.z;

    let scale = 1.0 / samples_per_pixel as f32;

    // sqrt is for gamma correction
    r = f32::sqrt(r * scale);
    g = f32::sqrt(g * scale);
    b = f32::sqrt(b * scale);

    Vec3 {
        x: clamp(r, 0.0, 0.999) * 256.0,
        y: clamp(g, 0.0, 0.999) * 256.0,
        z: clamp(b, 0.0, 0.999) * 256.0,
    }
}

/// Determine the color of a pixel for a given ray.
fn color_pixel(ray: &Ray, world: &HittableList, depth: i32) -> Vec3 {
    if depth <= 0 {
        return Vec3::new(0.0, 0.0, 0.0);
    }

    if let Some(rec) = world.hit(*ray, 0.01, 99999999999.0) {
        if let Some((attenuation, scattered)) = rec.mat.scatter(ray, &rec) {
            // Don't overload the * operator to do dot product...
            let res = color_pixel(&scattered, world, depth - 1);
            let vec = Vec3 {
                x: res.x * attenuation.x,
                y: res.y * attenuation.y,
                z: res.z * attenuation.z,
            };
            return vec;
        }
        
        return Color::new(0.0, 0.0, 0.0);
    }

    let unit_direction = ray.direction.unit_vector();
    let t = 0.5 * (unit_direction.y + 1.0);

    // gradient background
    Vec3::new(1.0, 1.0, 1.0) * (1.0 - t) + Vec3::new(0.5, 0.7, 1.0) * t
}

fn generate_world(world: &mut HittableList) {
    let material_ground = Lambertian {
        albedo: Color::new(0.8, 0.8, 0.0),
    };
    let material_center = Lambertian {
        albedo: Color::new(0.7, 0.3, 0.3),
    };
    let material_left = Metal {
        albedo: Color::new(0.8, 0.8, 0.8),
        fuzz: 0.0,
    };
    let material_right = Metal {
        albedo: Color::new(0.8, 0.6, 0.2),
        fuzz: 0.0,
    };
    world.add(Sphere::new(
        Point3::new(0.0, -100.5, -1.0),
        100.0,
        Rc::new(material_ground),
    ));
    world.add(Sphere::new(
        Point3::new(0.0, 0.0, -1.0),
        0.5,
        Rc::new(material_center),
    ));
    world.add(Sphere::new(
        Point3::new(-1.0, 0.0, -1.0),
        0.5,
        Rc::new(material_left),
    ));
    world.add(Sphere::new(
        Point3::new(1.0, 0.0, -1.0),
        0.5,
        Rc::new(material_right),
    ));
}

fn generate_small_world(world: &mut HittableList) {
    let material_center = Lambertian {
        albedo: Color::new(0.7, 0.3, 0.3),
    };
    let material_right = Metal {
        albedo: Color::new(0.8, 0.6, 0.2),
        fuzz: 0.0,
    };

    world.add(Sphere::new(
        Point3::new(0.0, 0.0, -1.0),
        1.0,
        Rc::new(material_right),
    ));

}

fn main() -> Result<(), Error> {
    let aspect_ratio = 3.0 / 2.0;
    let width: u32 = 640;
    let height: u32 = (width as f32 / aspect_ratio) as u32;

    env_logger::init();
    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let window = {
        let size = LogicalSize::new(width as f64, height as f64);
        WindowBuilder::new()
            .with_title("Hello Pixels")
            .with_inner_size(size)
            .with_min_inner_size(size)
            .build(&event_loop)
            .unwrap()
    };

    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(width, height, surface_texture)?
    };

    let samples_per_pixel = 50;

    let mut rng = rand::thread_rng();

    let lookfrom = Point3::new(13.0, 2.0, 3.0);
    let lookat = Point3::new(0.0, 0.0, 0.0);
    let vup = Vec3::new(0.0, 1.0, 0.0);

    let camera = Camera::new(
        lookfrom,
        lookat,
        vup,
        20.0,
        width as f32 / height as f32,
        0.1,
        10.0,
    );

    let mut world = HittableList::new();
    generate_world(&mut world);

    event_loop.run(move |event, _, control_flow| {
        // Draw the current frame
        if let Event::RedrawRequested(_) = event {
            let frame = pixels.get_frame();
            for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
                // x and y are pixel coordinates
                let x = (i % width as usize) as i16;
                let y = (i / width as usize) as i16;

                // image was displaying upside down for some reason.
                let y = height as i16 - y;

                let mut pixel_color = Vec3::new(0.0, 0.0, 0.0);
                for _s in 0..samples_per_pixel {
                    // u and v are the how far, as a percentage, x and y are from
                    // the vertical and horizontal of our viewport. This is used
                    // to map our pixel coords to the "camera" coords.
                    let u = (x as f32 + rng.gen::<f32>()) as f32 / (width - 1) as f32;
                    let v = (y as f32 + rng.gen::<f32>()) as f32 / (height - 1) as f32;

                    // origin is the camera (0, 0 ,0) and direction is the point in
                    // the viewport whose color value we are calculating.
                    let ray = camera.get_ray(u, v);
                    pixel_color += color_pixel(&ray, &world, 50);
                }

                let color = calculate_color(pixel_color, samples_per_pixel);
                let ir = (color.x) as u8;
                let ig = (color.y) as u8;
                let ib = (color.z) as u8;

                let rgba = [ir, ig, ib, 0xff];

                pixel.copy_from_slice(&rgba);
            }
            println!("FRAME");

            if pixels
                .render()
                .map_err(|e| error!("pixels.render() failed: {}", e))
                .is_err()
            {
                *control_flow = ControlFlow::Exit;
                return;
            }
        }


        // Handle input events
        if input.update(&event) {
            // Close events
            if input.key_pressed(VirtualKeyCode::Escape) || input.quit() {
                *control_flow = ControlFlow::Exit;
                return;
            }

            // Resize the window
            if let Some(size) = input.window_resized() {
                pixels.resize_surface(size.width, size.height);
            }

            // Update internal state and request a redraw
            // world.update();
            window.request_redraw();
        }
    });
}
