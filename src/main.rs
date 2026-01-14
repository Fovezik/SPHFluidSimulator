mod particle;
mod simulator;
mod grid;
mod config;
mod application;
use application::Application;

fn main() {
    let mut simulation_app = Application::new();
    simulation_app.run();
}