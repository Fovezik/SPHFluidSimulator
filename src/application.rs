use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use kiss3d::text::Font;
use kiss3d::nalgebra::{Translation3, Point3, Point2};
use crate::particle::Particle;
use crate::simulator::Simulator;
use crate::config::{BOX_SIZE, PARTICLES_PER_AXIS, START_POSITION, PARTICLE_COLOR, BOX_COLOR, VISCOSITY};

// Główna struktura aplikacji
pub struct Application {
    window: Window, // Główne okno aplikacji
    camera: ArcBall, // Kamera do nawigacji w scenie 3D
    solver: Simulator, // Symulator SPH
    particles: Vec<Particle>, // Wektor cząstek w symulacji
    visual_particles: Vec<SceneNode>, // Wektor wizualnych reprezentacji cząstek
    is_paused: bool, // Flaga pauzy symulacji
}

// Implementacja aplikacji
impl Application {
    // Konstruktor aplikacji
    pub fn new() -> Self {
        let mut window = Window::new("SPH Fluid Simulation");
        let eye = Point3::new(BOX_SIZE, BOX_SIZE * 2.0, BOX_SIZE * 3.0);
        let at = Point3::new(BOX_SIZE / 2.0, 0.0, 0.0);
        let camera = ArcBall::new(eye, at);
        window.set_light(Light::Absolute(Point3::new(BOX_SIZE * 2.0, BOX_SIZE * 2.0, BOX_SIZE * 2.0)));
        let particles = Self::initialize_particles();
        let solver = Simulator::new(BOX_SIZE, VISCOSITY);
        let visual_particles = Self::initialize_visual_particles(&mut window, &particles);
        Self::create_bounding_box(&mut window);
        Self {
            window,
            camera,
            solver,
            particles,
            visual_particles,
            is_paused: false,
        }
    }

    // Główna pętla uruchamiająca aplikację
    pub fn run(&mut self) {
        let dt = 0.008;
        while self.window.render_with_camera(&mut self.camera) {
            if !self.is_paused {
                self.solver.update(dt, &mut self.particles);
                for (particle, visual_particle) in
                    self.particles.iter().zip(self.visual_particles.iter_mut()) {
                    let translation = Translation3::new(
                        particle.position[0],
                        particle.position[1],
                        particle.position[2]
                    );
                    visual_particle.set_local_translation(translation);
                }
            }
            self.draw_text();
            if self.window.get_key(kiss3d::event::Key::Space) == kiss3d::event::Action::Press {
                self.is_paused = !self.is_paused;
            }
        }
    }

    // Inicjalizacja cząstek w scenie
    fn initialize_particles() -> Vec<Particle> {
        let mut particles = Vec::new();
        for x in 0..PARTICLES_PER_AXIS {
            for y in 0..PARTICLES_PER_AXIS {
                for z in 0..PARTICLES_PER_AXIS {
                    let px = START_POSITION.0 + (x as f32) * 0.6;
                    let py = START_POSITION.1 + (y as f32) * 0.6;
                    let pz = START_POSITION.2 + (z as f32) * 0.6;
                    particles.push(Particle::new(px, py, pz));
                }
            }
        }
        particles
    }

    // Inicjalizacja wizualnych reprezentacji cząstek
    fn initialize_visual_particles(window: &mut Window, particles: &[Particle]) -> Vec<SceneNode> {
        let mut visual_particles = Vec::with_capacity(particles.len());
        for _ in particles {
            let mut sphere = window.add_sphere(0.4);
            sphere.set_color(PARTICLE_COLOR.0, PARTICLE_COLOR.1, PARTICLE_COLOR.2);
            visual_particles.push(sphere);
        }
        visual_particles
    }

    // Utworzenie obramowania sceny
    fn create_bounding_box(window: &mut Window) {
        let half_size = BOX_SIZE / 2.0;
        let mut bounding_box = window.add_cube(BOX_SIZE, BOX_SIZE, BOX_SIZE);
        bounding_box.set_local_translation(Translation3::new(half_size, half_size, half_size));
        bounding_box.set_color(BOX_COLOR.0, BOX_COLOR.1, BOX_COLOR.2);
        bounding_box.set_lines_width(1.0);
        bounding_box.set_surface_rendering_activation(false);
    }

    // Rysowanie interfejsu użytkownika
    fn draw_text(&mut self) {
        let text = format!("Particles: {}", self.particles.len());
        let text2 = format!("Viscosity: {:.2}", self.solver.viscosity);
        self.window.draw_text(
            &text,
            &Point2::new(10.0, 10.0),
            50.0,
            &Font::default(),
            &Point3::new(1.0, 1.0, 1.0)
        );
        self.window.draw_text(
            &text2,
            &Point2::new(10.0, 60.0),
            50.0,
            &Font::default(),
            &Point3::new(1.0, 1.0, 1.0)
        );
    }
}