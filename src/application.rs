use kiss3d::camera::ArcBall;
use kiss3d::event::Action;
use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use kiss3d::text::Font;
use kiss3d::nalgebra::{Translation3, Point3, Point2};
use crate::particle::Particle;
use crate::simulator::Simulator;
use crate::config::{BOX_SIZE, PARTICLES_PER_AXIS, START_POSITION, PARTICLE_COLOR, BOX_COLOR,
                    VISCOSITY, RADIUS, REST_DENSITY, GAS_CONST, TIME_STEP, WINDOW_SIZE};

// Główna struktura aplikacji
pub struct Application {
    window: Window, // Główne okno aplikacji
    camera: ArcBall, // Kamera do nawigacji w scenie 3D
    simulator: Simulator, // Symulator SPH
    particles: Vec<Particle>, // Wektor cząstek w symulacji
    visual_particles: Vec<SceneNode>, // Wektor wizualnych reprezentacji cząstek
    is_paused: bool, // Flaga pauzy symulacji
}

// Implementacja aplikacji
impl Application {
    // Konstruktor aplikacji
    pub fn new() -> Self {
        let mut window = Window::new_with_size("SPH Fluid Simulation", WINDOW_SIZE.0, WINDOW_SIZE.1);
        let eye = Point3::new(BOX_SIZE, BOX_SIZE, BOX_SIZE * 3.0);
        let at = Point3::new(BOX_SIZE / 2.0, 0.0, 0.0);
        let camera = ArcBall::new(eye, at);
        window.set_light(Light::Absolute(Point3::new(BOX_SIZE * 2.0, BOX_SIZE * 2.0, BOX_SIZE * 2.0)));
        let particles = Self::initialize_particles();
        let simulator = Simulator::new(BOX_SIZE, VISCOSITY, REST_DENSITY, GAS_CONST, RADIUS);
        let visual_particles = Self::initialize_visual_particles(&mut window, &particles);
        Self::create_bounding_box(&mut window);
        Self {
            window,
            camera,
            simulator,
            particles,
            visual_particles,
            is_paused: false,
        }
    }

    // Główna pętla uruchamiająca aplikację
    pub fn run(&mut self) {
        while self.window.render_with_camera(&mut self.camera) {
            if !self.is_paused {
                self.simulator.update(TIME_STEP, &mut self.particles);
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
            for event in self.window.events().iter() {
                if let kiss3d::event::WindowEvent::Key(kiss3d::event::Key::Space, Action::Press, _) = event.value {
                    self.is_paused = !self.is_paused;
                }
                if let kiss3d::event::WindowEvent::Key(kiss3d::event::Key::Escape, Action::Press, _) = event.value {
                    self.window.close();
                }
            }
        }
    }

    // Inicjalizacja cząstek w scenie
    fn initialize_particles() -> Vec<Particle> {
        let spacing = 1.5 * RADIUS;
        let mut particles = Vec::new();
        for x in 0..PARTICLES_PER_AXIS {
            for y in 0..PARTICLES_PER_AXIS {
                for z in 0..PARTICLES_PER_AXIS {
                    let px = START_POSITION.0 + (x as f32) * spacing;
                    let py = START_POSITION.1 + (y as f32) * spacing;
                    let pz = START_POSITION.2 + (z as f32) * spacing;
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
            let mut sphere = window.add_sphere(RADIUS);
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
        if self.is_paused{
            self.window.draw_text("Status: Paused", &Point2::new(10.0, 160.0),
                                  50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
            );
            self.window.draw_text("Press SPACE to Resume", &Point2::new(10.0, 210.0),
                                  50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
            );
        } else {
            self.window.draw_text("Status: Running...", &Point2::new(10.0, 160.0),
                                  50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
            );
            self.window.draw_text("Press SPACE to Pause", &Point2::new(10.0, 210.0),
                                  50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
            );
        }
        self.window.draw_text(&format!("Particles: {}", self.particles.len()), &Point2::new(10.0, 10.0),
                              50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
        );
        self.window.draw_text(&format!("Viscosity: {:.2}", self.simulator.viscosity), &Point2::new(10.0, 60.0),
                              50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
        );
        self.window.draw_text(&format!("Time Step: {:.4}", TIME_STEP), &Point2::new(10.0, 110.0),
                              50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
        );
        self.window.draw_text("Press ESC to exit", &Point2::new(10.0, 260.0),
                              50.0, &Font::default(), &Point3::new(1.0, 1.0, 1.0)
        );
    }
}