import numpy as np
from scipy import ndimage
from state_grid import StateGrid, EMPTY, SUSCEPTIBLE, INFECTED, RESISTANT, SENESCENT
from parameters import *
import numpy.typing as npt # if you try to import directly from numpy it doesn't work 
import random

class MatrixSimulation:
    """Matrix-based simulation for high performance."""
    
    def __init__(self, grid_size: int = 55):
        self.grid = StateGrid(grid_size)
        self.time = 0
        # Track infection history for trend analysis
        self.infection_history: list[int] = []  # Store infected counts over time
        
    def run_time_step(self) -> bool:
        """Run one complete time step. Returns True if simulation should continue."""
        self.time += 1
        self.grid.cell_ages += 1
        
        # Decay the virus and interferon after all other processes have been done
        self.grid.decay_virus()
        self.grid.decay_interferon()
        
        self.grid.release_virus()
        self.grid.diffuse_virus()
        
        self.grid.release_interferon()
        self.grid.diffuse_interferon()
        
        # Process in logical order to avoid race conditions
        self._process_infected_cells()
        
        self._process_susceptible_cells() 
        self._process_resistant_cells()
        self._process_senescent_cells()
         
        return not self._check_termination()
    
    def _process_infected_cells(self) -> None:
        """Update all infected cells: age, infectivity, interferon, death."""
        infected_mask: npt.NDArray[np.bool_] = self.grid.cell_states == INFECTED
        
        if not np.any(infected_mask):
            return
        
        # Update virus released
        self._update_virus_production(infected_mask)
        
        # Update interferon levels
        self._update_interferon_from_infected_released(infected_mask)
        
        # Check for cell death
        death_mask = infected_mask & (self.grid.cell_ages >= INFECTED_LIFE_LIMIT)
        if np.any(death_mask):
            self.grid.cell_states[death_mask] = EMPTY
            self.grid.cell_ages[death_mask] = 0
            self.grid.actual_virus_production[death_mask] = 0
            self.grid.actual_interferon_production[death_mask] = 0
            
    def _update_virus_production(self, infected_mask) -> None:
        """How many viruses does each infected cell release based on its living time."""
        ages: np.ndarray = self.grid.cell_ages[infected_mask]
        peak_age: int = INFECTED_LIFE_LIMIT // PEAK_VIRUS_PRODUCTION_PART_OF_LIFESPAN
        
        # Initialize output array
        virus_productions: np.ndarray = np.zeros_like(ages, dtype=np.int16)
        
         # Case 1: Young cells (1 to peak_age-1) - parabolic increase
        young_mask: np.ndarray = ages < peak_age
        young_ages: np.ndarray = ages[young_mask].astype(np.int16)
        
        if young_ages.size > 0:
            m: float = float(peak_age)
            k: float = (1.0 - VIRUS_PRODUCTION_START_RATIO) / (m * m) if m > 0 else 0.0
            p: np.ndarray = 1.0 - k * (young_ages - m) ** 2
            p = np.clip(p, 0.0, 1.0)
            virus_productions[young_mask] = np.round(p * BASE_VIRUS_RELEASE).astype(np.int16)
            
        # Case 2: Plateau period (peak_age to death-1) - constant max infectivity
        plateau_mask: np.ndarray = (ages >= peak_age) & (ages <= INFECTED_LIFE_LIMIT)
        virus_productions[plateau_mask] = INFECTED_LIFE_LIMIT
        
        self.grid.actual_virus_production[infected_mask] = virus_productions
        
    def _update_interferon_from_infected_released(self, infected_mask: np.ndarray) -> None:
        """Vectorized interferon update for infected cells.
        They are parabolically decreased during living time.
        This method is not applied to senescent cells, because senescent cells have constant BASE_INTERFERON_RELEASE."""
        ages = self.grid.cell_ages[infected_mask].astype(np.float32)
        T = float(INFECTED_LIFE_LIMIT)
        # For parabolic: a coefficient that controls the curvature
        # If INTERFERON_PRODUCTION_DEGENERATION_IN_INFECTED is available, use it as the exponent
        # For parabolic decay, we'll use 2 (quadratic) or adjust based on parameter
        exponent: float = INTERFERON_PRODUCTION_DEGENERATION_IN_INFECTED  # 2.0 = Quadratic parabolic decay
        
        # Parabolic formula: scaled = 1 - (ages/T)^exponent
        # This gives 1 at age 0 and 0 at age T
        scaled = 1.0 - np.power(ages / T, exponent)
        scaled = np.maximum(0.0, scaled)
        
        self.grid.actual_interferon_production[infected_mask] = np.round(BASE_INTERFERON_RELEASE * scaled).astype(int)
        
    def _process_susceptible_cells(self) -> None:
        """Process all susceptible cells: infection, antiviral transition, reproduction."""
        susceptible_positions = np.where(self.grid.cell_states == SUSCEPTIBLE)
        
        for i, j in zip(susceptible_positions[0], susceptible_positions[1]):
            if self._attempt_infection(i, j):
                continue
            if self._attempt_antiviral_transition(i, j):
                continue
            self._attempt_reproduction(i, j)
        
    def _attempt_infection(self, i: int, j: int) -> bool:
        """Attempt to infect a susceptible cell at position (i,j).
        INFECTIVITY is the probability of one virus particle to induce infection.
        Therefore, we multiply this on quantity of viruses from virus_surface variable and the product could be higher than 1.0."""
        
        if INFECTIVITY * self.grid.virus_surface[i,j] > np.random.rand():
            self.grid.cell_states[i,j] = INFECTED
            self.grid.cell_ages[i,j] = 0
            return True
        else:
            return False
        
    def _attempt_antiviral_transition(self, i: int, j: int) -> bool:
        """Attempt to transit susceptible cell to resistant based on interferon.
        We calculate the probability based on susceptible transition sensitivity to one interferon molecule.
        Therefore we multiply the sensitivity to the quantity of interferon molecules."""
        
        if INTERFERON_SENSITIVITY_SUSCEPTIBLE * self.grid.interferon_surface[i,j] > np.random.rand():
            self.grid.cell_states[i, j] = RESISTANT
            self.grid.cell_ages[i,j] = 0
            return True
        else:
            return False
        
    def _attempt_reproduction(self, i: int, j: int) -> None:
        """Attempt reproduction of SUSCEPTIBLE cell into nearby empty spots."""
        empty_neighbors: list[tuple[int,int]] = []
        
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                ni, nj = i + di, j + dj
                if 0 <= ni < self.grid.size and 0 <= nj < self.grid.size:
                    if self.grid.cell_states[ni, nj] == EMPTY:
                        empty_neighbors.append((ni, nj))
        
        if not empty_neighbors:
            return
        
        if len(empty_neighbors) == 1:
            ni, nj = empty_neighbors[0]
            self.grid.cell_states[ni, nj] = SUSCEPTIBLE
            self.grid.cell_ages[ni, nj] = 0
            return
        
        # Weight by empty neighbors of the target spot
        weights = []
        for ni, nj in empty_neighbors:
            empty_count = 0
            for ddi in (-1, 0, 1):
                for ddj in (-1, 0, 1):
                    if ddi == 0 and ddj == 0:
                        continue
                    xi, xj = ni + ddi, nj + ddj
                    if 0 <= xi < self.grid.size and 0 <= xj < self.grid.size:
                        if self.grid.cell_states[xi, xj] == EMPTY:
                            empty_count += 1
            weights.append(1 + empty_count)
        
        weights = np.array(weights, dtype=float)
        probs = weights / weights.sum()
        chosen_idx = np.random.choice(len(empty_neighbors), p=probs)
        ni, nj = empty_neighbors[chosen_idx]
        self.grid.cell_states[ni, nj] = SUSCEPTIBLE
        self.grid.cell_ages[ni, nj] = 0
    
    def _process_resistant_cells(self) -> None:
        """Process all resistant cells: transition to susceptible or senescent."""
        resistant_positions: tuple = np.where(self.grid.cell_states == RESISTANT)
        
        for i, j in zip(resistant_positions[0], resistant_positions[1]):
            if self._attempt_senescent_transition(i, j):
                continue
            self._attempt_susceptible_transition(i, j)
            
    def _attempt_senescent_transition(self, i: int, j: int) -> bool:
        """Attempt to transit resistant cell to senescent.
        Based on INTERFERON_SENSITIVITY_RESISTANT_TO_SENESCENT."""
        
        if INTERFERON_SENSITIVITY_RESISTANT_TO_SENESCENT * self.grid.interferon_surface[i,j] > np.random.rand():
            self.grid.cell_states[i,j] = SENESCENT
            self.grid.cell_ages[i,j] = 0
            return True
        else:
            return False
        
    def _attempt_susceptible_transition(self, i: int, j: int) -> bool:
        """Attempt to transit resistant cell to senescent.
        Based on INTERFERON_SENSITIVITY_RESISTANT_TO_SUSCEPTIBLE."""
        # Get the strongest interferon influence (most relevant for decision)
        if INTERFERON_SENSITIVITY_RESISTANT_TO_SENESCENT * self.grid.interferon_surface[i,j] <= np.random.rand():
            self.grid.cell_states[i,j] = SUSCEPTIBLE
            self.grid.cell_ages[i,j] = 0
            return True
        else:
            return False
        
    def _process_senescent_cells(self) -> None:
        """Process all senescent cells: update actual interferon production.
        For young senescent cells (age < AGE_OF_SENESCENT_TRANSITION_FROM_RESISTANT),
        check if they should revert to resistant state based on interferon sensitivity."""
        
        # Find all senescent cells
        senescent_mask: np.ndarray = self.grid.cell_states == SENESCENT
        senescent_indices = np.where(senescent_mask)
        
        if len(senescent_indices[0]) == 0:
            return
        
        # Get ages of senescent cells
        senescent_ages = self.grid.cell_ages[senescent_mask]
        
        # Identify young senescent cells (age < threshold)
        young_senescent_mask = senescent_ages < AGE_OF_SENESCENT_TRANSITION_FROM_RESISTANT
        
        if np.any(young_senescent_mask):
            # Extract indices of young senescent cells
            young_senescent_indices = tuple(arr[young_senescent_mask] for arr in senescent_indices)
            
            # Get interferon levels at young senescent cell positions
            interferon_levels = self.grid.interferon_surface[young_senescent_indices]
            
            # Check probability for each young senescent cell
            random_values = np.random.rand(len(young_senescent_indices[0]))
            revert_mask = INTERFERON_SENSITIVITY_SENESCENT_TO_MATURITY * interferon_levels <= random_values
            
            # Find indices of cells that should revert to resistant
            revert_indices = tuple(arr[revert_mask] for arr in young_senescent_indices)
            
            if len(revert_indices[0]) > 0:
                # Revert these cells to resistant state
                self.grid.cell_states[revert_indices] = RESISTANT
                
                # Set their interferon production to 0
                self.grid.actual_interferon_production[revert_indices] = 0
                
                self.grid.cell_ages[revert_indices] = 0
            
            # Update the senescent mask to exclude reverted cells
            # Keep only young senescent cells that didn't revert
            kept_young_mask = ~revert_mask
            senescent_mask[young_senescent_indices] = kept_young_mask
            senescent_indices = np.where(senescent_mask)
        
        # Process remaining (mature and young that passed check) senescent cells
        if len(senescent_indices[0]) > 0:
            self._update_senescent_interferon_production(senescent_positions=senescent_indices)
        
    def _update_senescent_interferon_production(self, senescent_positions: tuple) -> None:
        """Update interferon production for senescent cells, which basically means assign BASE_INTERFERON_RELEASE to them."""
        self.grid.actual_interferon_production[senescent_positions] = BASE_INTERFERON_RELEASE
            
    def _check_termination(self) -> bool:
        """Check if simulation should terminate."""
        infected_count = np.sum(self.grid.cell_states == INFECTED)
        empty_count = np.sum(self.grid.cell_states == EMPTY)
        senescent_count = np.sum(self.grid.cell_states == SENESCENT)
        
        if infected_count + empty_count + senescent_count == self.grid.size * self.grid.size:
            print("Simulation finished: only empty spots, senescent cells or infected cells left.")
            return True
        elif infected_count == 0:
            print("Simulation finished: no infected cells left.")
            return True
        
        return False
        
            
            