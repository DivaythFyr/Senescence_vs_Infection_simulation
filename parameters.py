# parameters_4776.py
# simulation related parameters:
MAX_SIMULATION_TIMEPOINTS: int = 365000 # maximum timepoints after which the code would stop even if finish state haven't been achieved
GRID_SIZE: int = 55 # always set the number uneven, otherwise there woudn't be center

BASE_VIRUS_RELEASE: int = 102 # how many particles at maximum are released by one infected cell one timepoint
BASE_INTERFERON_RELEASE: int = 197 # how many particles at maximum are released by one infected or senescent cell one timepoint

# parameters related to virus production, to _update_virus_production
PEAK_VIRUS_PRODUCTION_PART_OF_LIFESPAN: int = 50 # at which part of life infected cell starts to produce
# BASE_VIRUS_RELEASE quantity of virus particles. If it is 2, means at the middle point of timelife
VIRUS_PRODUCTION_START_RATIO: float = 0.45 # Initial production of viruses as fraction of max (age 0)

# parameters related to interferon production, to _update_interferon_from_infected_released
INTERFERON_PRODUCTION_DEGENERATION_IN_INFECTED: float = 1.7 # Controls how quickly interferon production decays as infected cells age.
# Higher values mean faster decay (more exponential falloff). could be higher than 1.0.
# 2.0 = parabolic decrease

# Infection related parameters
INFECTIVITY: float = 0.17 # the probability of one viral particle to induce infection
INFECTED_LIFE_LIMIT: int = 30 # after which timepoint an infected cell dies

INTERFERON_SENSITIVITY_SUSCEPTIBLE: float = 0.978 # the probability of one interferon molecule to induce transition in healthy cell to resistant

INTERFERON_SENSITIVITY_RESISTANT_TO_SENESCENT: float = 0.045 # the probability of one interferon molecule to induce transition in resistant cell to senescent
INTERFERON_SENSITIVITY_RESISTANT_TO_SUSCEPTIBLE: float = 0.13 # the probability of one interferon molecule to induce transition in resistant cell to susceptible
INTERFERON_SENSITIVITY_SENESCENT_TO_MATURITY: float = 0.055 # the probability of one interferon molecule to prevent senescent cell to go back to resistant

VIRUS_DECAY_PER_TIME: float = 0.222 # how much percents does one grid point looses viral particles per timepoint. from 0.0 to 1.0
INTERFERON_DECAY_PER_TIME: float = 0.267 # how much percents does one grid point looses interferon molecules per timepoint. from 0.0 to 1.0

VIRUS_DIFFUSION_RATE: float = 0.468  # Fraction of viruses that diffuse to neighbors (0.0-1.0)
VIRUS_DIFFUSION_SIGMA: float = 1.0     # Moderate spread distance (how many cells are affected by the diffusion source)
VIRUS_EDGE_DIFFUSION_FILLING : float = 0.584 # What would be add on the edges of diffusion. Keep this value float but without residues, like 0.0 1.0

INTERFERON_DIFFUSION_RATE: float = 0.768  # Fraction of interferons that diffuse to neighbors (0.0-1.0)
INTERFERON_DIFFUSION_SIGMA: float = 1.6 # Interferons spread wider (how many cells are affected by the diffusion source)
INTERFERON_EDGE_DIFFUSION_FILLING : float = 4.06 # What would be add on the edges of diffusion. Keep this value float but without residues, like 0.0 1.0

AGE_OF_SENESCENT_TRANSITION_FROM_RESISTANT: int = 4 # Number of days after which senescent cell produces interferon.
# Before that time SENESCENT CELL is checked each day is interferon enough to stay in senescence, otherwise it comes back to resistant