import numpy as np
import numpy.typing as npt
from scipy import ndimage
from parameters import *

# Cell state constants
EMPTY = 0
SUSCEPTIBLE = 1
INFECTED = 2
RESISTANT = 3
SENESCENT = 4

class StateGrid:
    """Represents the grid of cells using numerical matrices for performance."""
    
    def __init__(self, size: int = 55):
        self.size = size
        self.time = 0
        
        # Main state matrix: stores cell types as integers
        self.cell_states = np.full((size, size), SUSCEPTIBLE, dtype=np.uint8)
        
        # Additional state matrices
        self.cell_ages: npt.NDArray[np.int16] = np.zeros((size, size), dtype=np.int16)
        self.infection_times: npt.NDArray[np.int16]  = np.zeros((size, size), dtype=np.int16)
        
        # specify exact new produced entities on the spot of cells that have done it, used each timepoint 
        self.actual_virus_production: npt.NDArray[np.int16]  = np.zeros((size, size), dtype=np.int16)
        self.actual_interferon_production: npt.NDArray[np.int16] = np.zeros((size, size), dtype=np.int16)
        
        # how many viruses and interferons distribute per timepoint by diffusion
        self.virus_surface: npt.NDArray[np.int16] = np.zeros((size, size), dtype=np.int16)
        self.interferon_surface: npt.NDArray[np.int16] = np.zeros((size, size), dtype=np.int16)
        
        # Initialize center as infected
        center = size // 2
        self.cell_states[center, center] = INFECTED
        self.infection_times[center, center] = 0
        
    def count_cell_type(self, cell_type: int) -> int:
        """Count the number of cells of a specific type."""
        return np.sum(self.cell_states == cell_type)
    
    def get_neighbors_mask(self, radius: int = 1) -> np.ndarray:
        """Create a mask for neighbors at given radius (Moore neighborhood)."""
        mask = np.ones((2*radius + 1, 2*radius + 1), dtype=bool)
        mask[radius, radius] = False  # Exclude center
        return mask
    
    def release_virus(self) -> None:
        """Add values from self.actual_virus_production to self.virus_surface on the positions of infected cells.
        Then there will be diffusion by another method on the same array"""
        self.virus_surface += self.actual_virus_production
        
    def release_interferon(self) -> None:
        """Add values from self.actual_interferon_production to self.interferon_surface on the positions of infected cells.
        Then there will be diffusion by another method on the same array"""
        self.interferon_surface += self.actual_interferon_production
        
    def decay_virus(self) -> None:
        """Substract values from self.virus_surface. Before diffusion.
        Affects only those points on the grid which have at least on virus particle."""
        decay_positions: tuple = np.where(self.virus_surface > 0)
        temporal_virus_surface = self.virus_surface.astype(np.float32).copy()
        temporal_virus_surface[decay_positions] -= temporal_virus_surface[decay_positions] * VIRUS_DECAY_PER_TIME
        self.virus_surface = np.floor(temporal_virus_surface).astype(np.int16)
        self.virus_surface = self.virus_surface.clip(min=0)
        
    def decay_interferon(self) -> None:
        """Substract values from self.interferon_surface. Before diffusion.
        Affects only those points on the grid which have at least on interferon molecule."""
        decay_positions: tuple = np.where(self.interferon_surface > 0)
        temporal_interferon_surface = self.interferon_surface.astype(np.float32).copy()
        temporal_interferon_surface[decay_positions] -= temporal_interferon_surface[decay_positions] * INTERFERON_DECAY_PER_TIME
        self.interferon_surface = np.floor(temporal_interferon_surface).astype(np.int16)
        self.interferon_surface = self.interferon_surface.clip(min=0)
        
    def diffuse_virus(self, diffusion_rate: float = VIRUS_DIFFUSION_RATE, sigma: float = VIRUS_DIFFUSION_SIGMA, cval: float = VIRUS_EDGE_DIFFUSION_FILLING) -> None:
        """Diffuse viruses using Gaussian filter (very fast and smooth).
        Args:
            diffusion_rate: Fraction of viruses that diffuse away from each position (0.0-1.0)
            sigma: standard deviation that defies how different is center with other cells. Could be higher than 1.0.
            cval: what would be added at the end of diffusion circle. Better be 0.0.
            """
        if diffusion_rate <= 0:
            return
        
        # Calculate viruses that stay vs diffuse
        viruses_staying = self.virus_surface * (1.0 - diffusion_rate)
        viruses_diffusing = self.virus_surface * diffusion_rate
        
        # Use Gaussian filter for smooth diffusion
        
        diffused = ndimage.gaussian_filter(viruses_diffusing.astype(np.float32), 
                                        sigma=sigma,  # Controls spread distance
                                        mode='constant', # 'constant': Pixels outside are filled with a constant value defined by cval
                                        cval=cval, # which value would be added on the edges. Better to keep value without residue, as 0.0 1.0 2.0 and etc
                                        )
        
        # Normalize to conserve mass (Gaussian filter doesn't automatically conserve)
        current_total = np.sum(viruses_diffusing)
        diffused_total = np.sum(diffused)
        if diffused_total > 0:
            diffused = diffused * (current_total / diffused_total)
        
        self.virus_surface = np.floor(viruses_staying).astype(np.int16) + np.floor(diffused).astype(np.int16)
        
    def diffuse_interferon(self, diffusion_rate: float = INTERFERON_DIFFUSION_RATE, sigma: float = INTERFERON_DIFFUSION_SIGMA, cval: float = INTERFERON_EDGE_DIFFUSION_FILLING) -> None:
        """Diffuse interferons using Gaussian filter.
        Args:
            diffusion_rate: Fraction of interferons that diffuse away from each position (0.0-1.0)
            sigma: standard deviation that defies how different is center with other cells. Could be higher than 1.0.
            cval: what would be added at the end of diffusion circle. Better be 0.0.
        """
        if diffusion_rate <= 0:
            return
        
        interferons_staying = self.interferon_surface * (1.0 - diffusion_rate)
        interferons_diffusing = self.interferon_surface * diffusion_rate
        
        from scipy import ndimage
        diffused = ndimage.gaussian_filter(interferons_diffusing.astype(np.float32), 
                                        sigma=sigma,
                                        mode='constant',
                                        cval=cval
                                        )
        
        # Normalize to conserve mass
        current_total = np.sum(interferons_diffusing)
        diffused_total = np.sum(diffused)
        if diffused_total > 0:
            diffused = diffused * (current_total / diffused_total)
        
        self.interferon_surface = np.floor(interferons_staying).astype(np.int16) + np.floor(diffused).astype(np.int16)
        
    def get_cell_counts(self) -> dict[str,int]:
        """Get counts of each cell type."""
        return {
            'empty': int(np.sum(self.cell_states == EMPTY)),
            'susceptible': int(np.sum(self.cell_states == SUSCEPTIBLE)),
            'infected': int(np.sum(self.cell_states == INFECTED)), 
            'resistant': int(np.sum(self.cell_states == RESISTANT)),
            'senescent': int(np.sum(self.cell_states == SENESCENT))
        }
        
        
    