use kiss3d::nalgebra::{Vector3, Point3};
use std::f32::consts::PI;
use rayon::prelude::*;
use crate::particle::Particle;
use crate::grid::Grid;

// Symulator SPH do symulacji płynów
pub struct Simulator {
    pub gravity: Vector3<f32>, // Siła grawitacji działająca na cząstki
    pub box_size: f32, // Rozmiar pudełka symulacji
    pub wall_damping: f32, // Współczynnik tłumienia przy kolizji ze ścianami
    pub h: f32, // Promień oddziaływania cząstek
    pub mass: f32, // Masa pojedynczej cząstki
    pub rest_density: f32, // Gęstość spoczynkowa płynu
    pub gas_const: f32, // Stała gazowa używana do obliczania ciśnienia
    pub viscosity: f32, // Współczynnik lepkości płynu
    grid: Grid, // Siatka przestrzenna do przyspieszenia wyszukiwania sąsiadów
}

// Implementacja symulatora SPH
impl Simulator {
    // Konstruktor symulatora
    pub fn new(box_size: f32, viscosity: f32) -> Self {
        let h = 1.2;
        Self {
            gravity: Vector3::new(0.0, -9.81, 0.0),
            box_size,
            wall_damping: 0.1,
            h,
            mass: 1.5,
            rest_density: 10.0,
            gas_const: 2000.0,
            viscosity,
            grid: Grid::new(h),
        }
    }

    // Główna funkcja aktualizująca symulację
    pub fn update(&mut self, dt: f32, particles: &mut [Particle]) {
        self.update_grid(particles);
        let density_pressure_data = self.calculate_density_and_pressure(particles);
        self.apply_density_and_pressure(particles, &density_pressure_data);
        let forces = self.calculate_forces(particles);
        self.integrate_and_collide(dt, particles, &forces);
    }

    // Aktualizacja siatki przestrzennej z cząstkami
    fn update_grid(&mut self, particles: &[Particle]) {
        self.grid.clear();
        self.grid.insert_particles(particles);
    }

    // Obliczanie gęstości i ciśnienia dla każdej cząstki
    fn calculate_density_and_pressure(&self, particles: &[Particle]) -> Vec<(f32, f32)> {
        particles.par_iter().map(|particle| {
            let mut neighbors = Vec::with_capacity(32);
            self.grid.get_neighbors(particle.position, &mut neighbors);
            let mut density_sum = 0.0;
            density_sum += self.mass * self.poly6(0.0);
            for &i in &neighbors {
                let particle_i = &particles[i];
                let distance_squared = (particle.position - particle_i.position).norm_squared();
                if distance_squared > 0.000001 && distance_squared < self.h * self.h {
                    density_sum += self.mass * self.poly6(distance_squared);
                }
            }
            let density = density_sum.max(0.001);
            let pressure = self.gas_const * (density - self.rest_density).max(0.0);
            (density, pressure)
        }).collect()
    }

    // Zastosowanie obliczonej gęstości i ciśnienia do cząstek
    fn apply_density_and_pressure(&self, particles: &mut [Particle], data: &[(f32, f32)]) {
        particles.par_iter_mut().zip(data.par_iter()).for_each
        (|(particle, (density, pressure))| {
            particle.density = *density;
            particle.pressure = *pressure;
        });
    }

    // Obliczanie sił działających na każdą cząstkę
    fn calculate_forces(&self, particles: &[Particle]) -> Vec<Vector3<f32>> {
        particles.par_iter().map(|p| {
            let mut total_force = self.gravity * self.mass;
            let mut neighbors = Vec::with_capacity(32);
            self.grid.get_neighbors(p.position, &mut neighbors);

            for &i in &neighbors {
                let particle_i = &particles[i];
                let difference = p.position - particle_i.position;
                let distance_squared = difference.norm_squared();

                // Unikanie dzielenia przez zero i sprawdzanie zasięgu
                if distance_squared > 0.000001 && distance_squared < self.h * self.h {
                    let distance = distance_squared.sqrt();

                    // Siła wynikająca z ciśnienia
                    let scalar = self.mass * (p.pressure + particle_i.pressure) / (2.0 * particle_i.density);
                    let grad_w = difference.normalize() * self.pressure_gradient(distance);
                    total_force -= grad_w * scalar;

                    // Siła wynikająca z lepkości
                    let velocity_difference = particle_i.velocity - p.velocity;
                    let visc_laplacian = self.viscosity_laplacian(distance);
                    total_force += velocity_difference * (self.mass / particle_i.density * self.viscosity * visc_laplacian);
                }
            }
            total_force
        }).collect()
    }

    // Integracja ruchu cząstek i obsługa kolizji ze ścianami
    fn integrate_and_collide(&self, dt: f32, particles: &mut [Particle], forces: &[Vector3<f32>]) {
        particles.par_iter_mut().zip(forces.par_iter()).for_each(|(particle, force)| {
            particle.force = *force;

            if particle.position[0].is_nan() {
                particle.position = Point3::origin();
                particle.velocity = Vector3::zeros();
            }

            let acceleration = force / self.mass;
            particle.velocity += acceleration * dt;
            particle.position += particle.velocity * dt;

            self.check_collisions(particle);
        });
    }

    // Sprawdzanie kolizji cząstek ze ścianami pudełka
    fn check_collisions(&self, particle: &mut Particle) {
        for i in 0..3 {
            if particle.position[i] < 0.0 {
                particle.position[i] = 0.0;
                particle.velocity[i] *= -self.wall_damping;
            } else if particle.position[i] > self.box_size {
                particle.position[i] = self.box_size;
                particle.velocity[i] *= -self.wall_damping;
            }
        }
    }

    // Funkcje jądra SPH
    fn poly6(&self, r_squared: f32) -> f32 {
        let h2 = self.h * self.h;
        let h9 = self.h.powi(9);
        let coefficient = 315.0 / (64.0 * PI * h9);
        if r_squared < h2 {
            let term = h2 - r_squared;
            return coefficient * term * term * term;
        }
        0.0
    }

    // Gradient do obliczania ciśnienia
    fn pressure_gradient(&self, r: f32) -> f32 {
        let h6 = self.h.powi(6);
        let coefficient = -45.0 / (PI * h6);
        if r < self.h {
            let term = self.h - r;
            return coefficient * term * term;
        }
        0.0
    }

    // Laplasjan do obliczania lepkości
    fn viscosity_laplacian(&self, r: f32) -> f32 {
        let h6 = self.h.powi(6);
        let coefficient = 45.0 / (PI * h6);
        if r < self.h { return coefficient * (self.h - r); }
        0.0
    }
}