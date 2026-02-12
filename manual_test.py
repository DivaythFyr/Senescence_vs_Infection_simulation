import numpy as np
import numpy.typing as npt
from state_grid import StateGrid

# check how matrix summation does work
# actual_virus_production: npt.NDArray[np.int16] = np.zeros((10, 10), dtype=np.int16)
# actual_virus_production[1, 2] = 150
# actual_virus_production[5, 9] = 150

# virus_surface: npt.NDArray[np.int16]  = np.zeros((10, 10), dtype=np.int16)
# print(virus_surface)

# virus_surface += actual_virus_production

# print(virus_surface)


# Check how multiplication of matrix on one number does work
# virus_surface: npt.NDArray[np.int16]  = np.full((10, 10),10, dtype=np.int16)

# print(virus_surface * 0.2)


# Check how clipping does happen
grid = StateGrid(size=5)
grid.virus_surface[0, 4] = 150
grid.virus_surface[3, 4] = 150
grid.decay_virus()

print(grid.virus_surface)
