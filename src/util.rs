#![deny(clippy::all)]
#![forbid(unsafe_code)]

use std::{
    rc::Rc,
};

use nalgebra::Vector3;
use rand::Rng;

use crate::material;
use material::Material;

pub fn length_squared(v: Vector3<f32>) -> f32 {
    v.magnitude() * v.magnitude()
}

pub fn unit_vector(v: Vector3<f32>) -> Vector3<f32> {
    let len = v.magnitude();
    Vector3::new(
        v.x / len,
        v.y / len,
        v.z / len,
    )
}

pub fn random_vector3(min: f32, max: f32) -> Vector3<f32> {
    let mut rng = rand::thread_rng();
    Vector3::new(
        rng.gen_range(min..max),
        rng.gen_range(min..max),
        rng.gen_range(min..max),
    )
}

pub fn random_in_unit_sphere() -> Vector3<f32> {
    loop {
        let p = random_vector3(-1.0, 1.0);
        if length_squared(p) >= 1.0 {
            continue;
        }
        return p;
    }
}

pub fn random_unit_vector() -> Vector3<f32> {
    // Not sure about this. Should maybe just be some util function outside of the impl.
    unit_vector(random_in_unit_sphere())
}

pub fn near_zero(v: Vector3<f32>) -> bool {
    // return true if vector is close to zero is all dims
    const S: f32 = 1e-8;
    return f32::abs(v.x) < S && f32::abs(v.y) < S && f32::abs(v.z) < S;
}

pub fn reflect(v: Vector3<f32>, n: Vector3<f32>) -> Vector3<f32> {
    v - n * ((v.dot(&n)) * 2.0)
}

pub fn refract(v: Vector3<f32>, normal: Vector3<f32>, etai_over_etat: f32) -> Vector3<f32> {
    let cos_theta = f32::min(-v.dot(&normal), 1.0);
    let r_out_perp = (v + (normal * cos_theta)) * etai_over_etat;
    let r_out_parallel = normal * -f32::sqrt(f32::abs(1.0 - length_squared(r_out_perp)));

    r_out_perp + r_out_parallel
}

pub type Point3 = Vector3<f32>;
pub type Color = Vector3<f32>;

#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub origin: Point3,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn at(self, t: f32) -> Point3 {
        self.origin + (self.direction * t)
    }
}

#[derive(Clone)]
pub struct HitRecord {
    pub p: Point3,
    pub normal: Vector3<f32>,
    pub mat: Rc<dyn Material>,
    pub t: f32,
    pub front_face: bool,
}

impl HitRecord {
    pub fn set_face_normal(&mut self, ray: &Ray, outward_normal: Vector3<f32>) {
        // dot product of ray and outward normal tells us if we are hitting
        // the inside or ouside of the surface.
        self.front_face = (ray.direction.dot(&outward_normal)) < 0.0;
        self.normal = if self.front_face {
            outward_normal
        } else {
            -outward_normal
        };
    }
}

pub trait Hittable {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord>;
}

pub struct Triangle {
    pub v1: Point3,
    pub v2: Point3,
    pub v3: Point3,
    pub mat: Rc<dyn Material>,
}

impl Triangle {
    pub fn new(v1: Point3, v2: Point3, v3: Point3, mat: Rc<dyn Material>) -> Triangle {
        Triangle { v1, v2, v3, mat }
    }
}

impl Hittable for Triangle {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let aP: f32 = self.v1.x - self.v2.x;
        let bP: f32 = self.v1.y - self.v2.y;
        let cP: f32 = self.v1.z - self.v2.z;
        let dP: f32 = self.v1.x - self.v3.x;
        let eP: f32 = self.v1.y - self.v3.y;
        let fP: f32 = self.v1.z - self.v3.z;

        let gP: f32 = ray.direction.x;
        let hP: f32 = ray.direction.y;
        let iP: f32 = ray.direction.z;

        let jP: f32 = self.v1.x - ray.origin.x;
        let kP: f32 = self.v1.y - ray.origin.y;
        let lP: f32 = self.v1.z - ray.origin.z;

        let m = (aP * ((eP * iP) - (hP * fP)))
            + (bP * ((gP * fP) - (dP * iP)))
            + (cP * ((dP * hP) - (eP * gP)));

        let mut beta = (jP * ((eP * iP) - (hP * fP)))
            + (kP * ((gP * fP) - (dP * iP)))
            + (lP * ((dP * hP) - (eP * gP)));
        beta = beta / m;

        let mut gamma = (iP * ((aP * kP) - (jP * bP)))
            + (hP * ((jP * cP) - (aP * lP)))
            + (gP * ((bP * lP) - (kP * cP)));
        gamma = gamma / m;

        let mut t = (fP * ((aP * kP) - (jP * bP)))
            + (eP * ((jP * cP) - (aP * lP)))
            + (dP * ((bP * lP) - (kP * cP)));
        t = -t / m;

        if t < 0.0 || t < t_min || t > t_max { return None; }
        if gamma < 0.0 || gamma > 1.0 { return None; }
        if beta < 0.0 || beta > 1.0 - gamma { return None; }

        let at_ray = ray.at(t);
        let a = self.v3 - self.v1;
        let b = self.v2 - self.v1;
        let normal = a.cross(&b);
        let mut rec = HitRecord {
            t: t,
            p: at_ray,
            mat: self.mat.clone(),
            normal: unit_vector(normal),
            front_face: false,
        };
        rec.set_face_normal(&ray, normal);
        Some(rec)
    }
}

pub struct Sphere {
    pub center: Point3,
    pub radius: f32,
    pub mat: Rc<dyn Material>,
}

impl Sphere {
    pub fn new(center: Point3, radius: f32, material: Rc<dyn Material>) -> Sphere {
        Sphere {
            center,
            radius,
            mat: material,
        }
    }
}

impl Hittable for Sphere {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let oc: Vector3<f32> = ray.origin - self.center;
        let a = length_squared(ray.direction);
        let half_b = oc.dot(&ray.direction);
        let c = (length_squared(oc)) - (self.radius * self.radius);

        let discriminant: f32 = (half_b * half_b) - (a * c);

        // no roots (negative discriminant) = no intersction
        if discriminant < 0.0 {
            return None;
        }
        let sqrtd = f32::sqrt(discriminant);

        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || root > t_max {
            root = (-half_b + sqrtd) / a;
            if root < t_min || root > t_max {
                return None;
            }
        }

        let at_ray = ray.at(root);
        let mut rec = HitRecord {
            t: root,
            p: at_ray,
            mat: self.mat.clone(),
            normal: (at_ray - self.center) / self.radius,
            front_face: false,
        };

        let outward_normal = (rec.p - self.center) / self.radius;
        rec.set_face_normal(&ray, outward_normal);
        // hit_record.material = self.material;

        Some(rec)
    }
}

pub struct HittableList {
    // A vector of objects that implement the Hittable trait.
    // Box is needed because Hittable objects can be of different
    // sizes.
    objects: Vec<Box<dyn Hittable>>,
}

impl HittableList {
    pub fn new() -> HittableList {
        HittableList {
            objects: Vec::new(),
        }
    }

    pub fn add(&mut self, object: impl Hittable + 'static) {
        // I'm not 100% clear on if this is the correct way to do
        // this.
        self.objects.push(Box::new(object) as Box<dyn Hittable>);
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut temp_rec = None;
        let mut closest_so_far = t_max;

        for object in &self.objects {
            if let Some(rec) = object.hit(ray, t_min, closest_so_far) {
                closest_so_far = rec.t;
                temp_rec = Some(rec);
            }
        }

        temp_rec
    }
}
