#![deny(clippy::all)]
#![forbid(unsafe_code)]

// This seems incorrect.
use crate::util;
use util::{Point3, Ray};

use nalgebra::Vector3;

#[derive(Copy, Clone, Debug, Default)]
pub struct Camera {
    pub origin: Point3,
    lower_left_corner: Point3,
    horizontal: Vector3<f32>,
    vertical: Vector3<f32>,
}

impl Camera {
    pub fn new(aspect_ratio: f32) -> Camera {
        let viewport_height = 2.0;
        let viewport_width = aspect_ratio * viewport_height;

        let origin = Point3::new(0.0, 0.0, 0.0);
        let horizontal = Vector3::new(viewport_width, 0.0, 0.0);
        let vertical = Vector3::new(0.0, viewport_height, 0.0);
        let lower_left_corner =
            origin - (horizontal / 2.0) - (vertical / 2.0) - Vector3::new(0.0, 0.0, 1.0);

        Camera {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
        }
    }

    pub fn get_ray(self, s: f32, t: f32) -> Ray {
        Ray {
            origin: self.origin,
            direction: self.lower_left_corner + (self.horizontal * s) + (self.vertical * t)
                - self.origin,
        }
    }
}
