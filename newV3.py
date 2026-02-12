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
    Sample ONLY the 6 specified parameters for Monte Carlo exploration.
    All other parameters remain exactly as defined in parameters.py.
    
    Randomized parameters and ranges:
    - VIRUS_DECAY_PER_TIME:           0.0 → 0.5
    - INTERFERON_DECAY_PER_TIME:      0.0 → 0.5
    - VIRUS_DIFFUSION_RATE:           0.0 → 1.0
    - INTERFERON_DIFFUSION_RATE:      0.0 → 1.0
    - VIRUS_EDGE_DIFFUSION_FILLING:   0.0 → 5.0  (integer or float values allowed)
    - INTERFERON_EDGE_DIFFUSION_FILLING: 0.0 → 5.0
    """
    return {
        "VIRUS_DECAY_PER_TIME": rng.uniform(0.0, 0.5),
        "INTERFERON_DECAY_PER_TIME": rng.uniform(0.0, 0.5),
        
        "VIRUS_DIFFUSION_RATE": rng.uniform(0.0, 1.0),
        "INTERFERON_DIFFUSION_RATE": rng.uniform(0.0, 1.0),
        
        "VIRUS_EDGE_DIFFUSION_FILLING": rng.uniform(0.0, 5.0),
        "INTERFERON_EDGE_DIFFUSION_FILLING": rng.uniform(0.0, 5.0),
    }


def run_single_simulation(sampled_params: Dict[str, Any]) -> Dict[str, Any]:
    with override_parameters(sampled_params):
        sim = MatrixSimulation(grid_size=params.GRID_SIZE)

        # Run until natural termination or safety limit
        MAX_STEPS = 2000
        for step in range(MAX_STEPS):
            if not sim.run_time_step():
                break

        counts = sim.grid.get_cell_counts()

    # Build result row
    row: Dict[str, Any] = dict(sampled_params)  # The 6 varied parameters

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

    # Define consistent column order
    fieldnames = [
        "simulation_index",
        # Varied parameters
        "VIRUS_DECAY_PER_TIME",
        "INTERFERON_DECAY_PER_TIME",
        "VIRUS_DIFFUSION_RATE",
        "INTERFERON_DIFFUSION_RATE",
        "VIRUS_EDGE_DIFFUSION_FILLING",
        "INTERFERON_EDGE_DIFFUSION_FILLING",
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
    run_monte_carlo(n_simulations=30)  # Adjust number as needed