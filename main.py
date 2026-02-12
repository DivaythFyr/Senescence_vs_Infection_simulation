import numpy as np
from simulation_core import MatrixSimulation
from parameters import *
from visualization import *
import time
import numpy.typing as npt

def main():
    print("Starting Matrix-based Cellular Infection Simulation")
    
    simulation = MatrixSimulation(GRID_SIZE)
    timepoint_cell_states_matrices: list[np.ndarray] = []
    timepoint_virus_surface_matrices: list[np.ndarray] = []
    timepoint_interferon_surface_matrices: list[np.ndarray] = []
    
    start_time = time.time()
    
    for timepoint in range(MAX_SIMULATION_TIMEPOINTS):
        # Store current state for visualization
        timepoint_cell_states_matrices.append(simulation.grid.cell_states.copy())
        timepoint_virus_surface_matrices.append(simulation.grid.virus_surface.copy())
        timepoint_interferon_surface_matrices.append(simulation.grid.interferon_surface.copy())
        
        # Run simulation step
        should_continue = simulation.run_time_step()
        
        counts = simulation.grid.get_cell_counts()
        print(f"Timepoint {timepoint}: {counts}")
        
        if not should_continue:
            break
    
    timepoint_cell_states_matrices_array: npt.NDArray = np.array(timepoint_cell_states_matrices)
    timepoint_virus_surface_matrices_array: npt.NDArray = np.array(timepoint_virus_surface_matrices)
    timepoint_interferon_surface_matrices_array: npt.NDArray = np.array(timepoint_interferon_surface_matrices)
    
    # Create visualization
    print("Creating cell type animation...")
    visualize_infection_spread(timepoint_cell_states_matrices_array)
    print("Animation saved to ./output/infection_spread_animation.gif")
    
    print("Creating virus surface animation...")
    visualize_virus_surface(timepoint_virus_surface_matrices_array)
    print("Animation saved to ./output/virus_level_surface_animation.gif")
    
    print("Creating interferon surface animation...")
    visualize_interferon_surface(timepoint_interferon_surface_matrices_array)
    print("Animation saved to ./output/interferon_level_surface_animation.gif")
    
    end_time = time.time()
    print(f"Simulation completed in {end_time - start_time:.2f} seconds")
    print(f"Total timepoints: {simulation.time}")

if __name__ == "__main__":
    main()
