use std::collections::HashMap;
use kiss3d::nalgebra::Point3;
use crate::particle::Particle;

// Struktura siatki przestrzennej do przyspieszenia wyszukiwania sąsiadów
pub struct Grid {
    cell_size: f32, // Rozmiar pojedynczej komórki siatki
    cells: HashMap<(i32, i32, i32), Vec<usize>>, // Mapa komórek siatki do indeksów cząstek
}

// Implementacja siatki przestrzennej
impl Grid {
    // Konstruktor siatki
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
        }
    }

    // Czyszczenie siatki przed ponownym wstawianiem cząstek
    pub fn clear(&mut self) {
        for (_,cell) in self.cells.iter_mut() {
            cell.clear();
        }
    }

    // Wstawianie cząstek do siatki
    pub fn insert_particles(&mut self, particles: &[Particle]) {
        for (i, particle) in particles.iter().enumerate() {
            let key = self.get_cell_coordinates(particle.position);
            self.cells.entry(key).or_insert_with(|| Vec::with_capacity(20)).push(i);
        }
    }

    // Pobieranie sąsiadów dla danej pozycji cząstki
    pub fn get_neighbors(&self, position: Point3<f32>, neighbors: &mut Vec<usize>) {
        neighbors.clear();
        let (x1, y1, z1) = self.get_cell_coordinates(position);
        for x2 in -1..=1 {
            for y2 in -1..=1 {
                for z2 in -1..=1 {
                    let key = (x1 + x2, y1 + y2, z1 + z2);
                    if let Some(indices) = self.cells.get(&key) {
                        neighbors.extend_from_slice(indices);
                    }
                }
            }
        }
    }

    // Obliczanie współrzędnych komórki siatki dla danej pozycji
    fn get_cell_coordinates(&self, position: Point3<f32>) -> (i32, i32, i32) {
        (
            (position[0] / self.cell_size).floor() as i32,
            (position[1] / self.cell_size).floor() as i32,
            (position[2] / self.cell_size).floor() as i32,
        )
    }
}