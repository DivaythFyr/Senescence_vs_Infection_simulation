# monte_carlo_simulations.py

import csv
import os
import time
from typing import Dict, Any, List
from contextlib import contextmanager

import numpy as np

# Import the modules we need to patch
import parameters as params
import simulation_core
import state_grid
from simulation_core import MatrixSimulation


@contextmanager
def override_parameters(overrides: Dict[str, Any]):
    """
    Temporarily override parameters in ALL modules that have imported them.
    """
    original_values: Dict[str, List[tuple]] = {}
    
    targets = [params, simulation_core, state_grid]
    
    for name, new_val in overrides.items():
        original_values[name] = []
        for module in targets:
            if hasattr(module, name):
                old_val = getattr(module, name)
                original_values[name].append((module, old_val))
                setattr(module, name, new_val)
                
    try:
        yield
    finally:
        for name, saved_list in original_values.items():
            for module, old_val in saved_list:
                setattr(module, name, old_val)


def sample_parameters(rng: np.random.Generator) -> Dict[str, Any]:
    """
    Sample parameters with values rounded to appropriate decimal places.
    
    Original 6 parameters (rounded to 3 decimals):
    - Decay rates:          0.000 → 0.500
    - Diffusion rates:      0.000 → 1.000
    - Edge filling values:  0.000 → 5.000
    
    Additional parameters:
    - BASE_VIRUS_RELEASE: 100 → 400 (integer)
    - BASE_INTERFERON_RELEASE: 100 → 300 (integer)
    - INFECTIVITY: 0.1 → 1.0 (rounded to 3 decimals)
    - INTERFERON_SENSITIVITY_SUSCEPTIBLE: 0.1 → 1.0 (rounded to 3 decimals)
    - INFECTED_LIFE_LIMIT: 5 → 20 (integer)
    """
    return {
        # Original 6 parameters (3 decimals)
        "VIRUS_DECAY_PER_TIME": round(rng.uniform(0.000, 0.500), 3),
        "INTERFERON_DECAY_PER_TIME": round(rng.uniform(0.000, 0.500), 3),
        
        "VIRUS_DIFFUSION_RATE": round(rng.uniform(0.000, 1.000), 3),
        "INTERFERON_DIFFUSION_RATE": round(rng.uniform(0.000, 1.000), 3),
        
        "VIRUS_EDGE_DIFFUSION_FILLING": round(rng.uniform(0.000, 5.000), 3),
        "INTERFERON_EDGE_DIFFUSION_FILLING": round(rng.uniform(0.000, 5.000), 3),
        
        # Additional parameters for virus production
        "BASE_VIRUS_RELEASE": int(rng.integers(100, 401)),  # 100-400
        "BASE_INTERFERON_RELEASE": int(rng.integers(100, 301)),  # 100-300
        
        # Infection and interferon sensitivity parameters (3 decimals)
        "INFECTIVITY": round(rng.uniform(0.1, 1.0), 3),
        "INTERFERON_SENSITIVITY_SUSCEPTIBLE": round(rng.uniform(0.1, 1.0), 3),
        "INTERFERON_SENSITIVITY_RESISTANT_TO_SENESCENT": round(rng.uniform(0.001, 0.05), 3),
        "INTERFERON_SENSITIVITY_RESISTANT_TO_SUSCEPTIBLE": round(rng.uniform(0.05, 0.3), 3),
        "INTERFERON_SENSITIVITY_SENESCENT_TO_MATURITY": round(rng.uniform(0.01, 0.1), 3),
        
        # Cell lifespan parameters (integers)
        "INFECTED_LIFE_LIMIT": int(rng.integers(5, 21)),  # 5-20
        "AGE_OF_SENESCENT_TRANSITION_FROM_RESISTANT": int(rng.integers(2, 11)),  # 2-10
        
        # Virus production timing (integers and floats)
        "PEAK_VIRUS_PRODUCTION_PART_OF_LIFESPAN": int(rng.integers(2, 5)),  # 2-4
        "VIRUS_PRODUCTION_START_RATIO": round(rng.uniform(0.3, 0.9), 2),  # 0.3-0.9
        
        # Interferon production decay (floats)
        "INTERFERON_PRODUCTION_DEGENERATION_IN_INFECTED": round(rng.uniform(1.0, 3.0), 1),  # 1.0-3.0
        
        # Diffusion sigma parameters (floats)
        "VIRUS_DIFFUSION_SIGMA": round(rng.uniform(0.5, 2.5), 1),  # 0.5-2.5
        "INTERFERON_DIFFUSION_SIGMA": round(rng.uniform(0.5, 3.0), 1),  # 0.5-3.0
    }


def run_single_simulation(sampled_params: Dict[str, Any]) -> Dict[str, Any]:
    with override_parameters(sampled_params):
        sim = MatrixSimulation(grid_size=params.GRID_SIZE)

        # Safety limit to prevent infinite loops
        MAX_STEPS = 2000
        for step in range(MAX_STEPS):
            if not sim.run_time_step():
                break

        counts = sim.grid.get_cell_counts()

    # Build result row
    row: Dict[str, Any] = dict(sampled_params)  # All varied parameters

    # Add final outcomes
    row["final_time"] = sim.time
    row["total_cells"] = params.GRID_SIZE ** 2

    for cell_type, count in counts.items():
        row[f"cell_count_{cell_type}"] = count

    # Useful derived metrics
    row["infected_fraction"] = counts["infected"] / row["total_cells"]
    row["cleared"] = (counts["infected"] == 0)

    return row


def run_monte_carlo(
    n_simulations: int = 200,
    output_csv_path: str = "./output/monte_carlo_results.csv",
    random_seed: int = 42,
) -> None:
    rng = np.random.default_rng(random_seed)
    rows: List[Dict[str, Any]] = []

    t0 = time.time()
    for i in range(n_simulations):
        sampled = sample_parameters(rng)
        row = run_single_simulation(sampled)
        row["simulation_index"] = i
        rows.append(row)

        inf_count = row["cell_count_infected"]
        total = row["total_cells"]
        cleared = "Yes" if row["cleared"] else "No"
        print(f"Run {i+1:3d}/{n_simulations} | Infected: {inf_count:4d}/{total} | Cleared: {cleared} | Time: {row['final_time']:4d}")

    # Consistent column order
    fieldnames = [
        "simulation_index",
        
        # Original 6 parameters
        "VIRUS_DECAY_PER_TIME",
        "INTERFERON_DECAY_PER_TIME",
        "VIRUS_DIFFUSION_RATE",
        "INTERFERON_DIFFUSION_RATE",
        "VIRUS_EDGE_DIFFUSION_FILLING",
        "INTERFERON_EDGE_DIFFUSION_FILLING",
        
        # Virus/interferon production parameters
        "BASE_VIRUS_RELEASE",
        "BASE_INTERFERON_RELEASE",
        
        # Infection parameters
        "INFECTIVITY",
        "INFECTED_LIFE_LIMIT",
        
        # Interferon sensitivity parameters
        "INTERFERON_SENSITIVITY_SUSCEPTIBLE",
        "INTERFERON_SENSITIVITY_RESISTANT_TO_SENESCENT",
        "INTERFERON_SENSITIVITY_RESISTANT_TO_SUSCEPTIBLE",
        "INTERFERON_SENSITIVITY_SENESCENT_TO_MATURITY",
        
        # Virus production timing
        "PEAK_VIRUS_PRODUCTION_PART_OF_LIFESPAN",
        "VIRUS_PRODUCTION_START_RATIO",
        
        # Interferon production decay
        "INTERFERON_PRODUCTION_DEGENERATION_IN_INFECTED",
        
        # Diffusion sigma parameters
        "VIRUS_DIFFUSION_SIGMA",
        "INTERFERON_DIFFUSION_SIGMA",
        
        # Senescence transition
        "AGE_OF_SENESCENT_TRANSITION_FROM_RESISTANT",
        
        # Outcomes
        "final_time",
        "total_cells",
        "cell_count_empty",
        "cell_count_susceptible",
        "cell_count_infected",
        "cell_count_resistant",
        "cell_count_senescent",
        "infected_fraction",
        "cleared",
    ]

    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    with open(output_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    elapsed = time.time() - t0
    cleared_count = sum(1 for r in rows if r["cleared"])
    print(f"\nMonte Carlo completed: {n_simulations} simulations in {elapsed:.2f} seconds")
    print(f"Cleared infection in {cleared_count}/{n_simulations} runs ({cleared_count/n_simulations*100:.1f}%)")
    print(f"Results saved to: {output_csv_path}")


if __name__ == "__main__":
    # You can adjust the number of simulations as needed
    run_monte_carlo(n_simulations=50)