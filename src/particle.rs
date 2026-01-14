use kiss3d::nalgebra::{Point3, Vector3};

// Struktura reprezentująca pojedynczą cząstkę w symulacji SPH
#[derive(Clone, Debug)]
pub struct Particle {
    pub position: Point3<f32>, // Pozycja cząstki w przestrzeni 3D
    pub velocity: Vector3<f32>, // Prędkość cząstki
    pub force: Vector3<f32>, // Siła działająca na cząstkę
    pub density: f32, // Gęstość cząstki
    pub pressure: f32, // Ciśnienie cząstki
}

// Implementacja metod dla struktury Particle
impl Particle {
    // Konstruktor tworzący nową cząstkę z podanymi współrzędnymi
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            position: Point3::new(x, y, z),
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            density: 0.0,
            pressure: 0.0,
        }
    }
}