import numpy as np
from matplotlib import pyplot as plt 
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
import os
# from state_grid import EMPTY, SUSCEPTIBLE, INFECTED, RESISTANT, SENESCENT
# from simulation_core import MatrixSimulation
from parameters import *
import numpy.typing as npt

def visualize_infection_spread(timepoint_cell_states_matrices_array: npt.NDArray[np.uint8]):
    """Visualizes the infection spread over time using matrix data."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Define colors in the correct order (0-4)
    colors = ['grey', 'pink', 'red', 'brown', 'gold']  # EMPTY=0, SUSCEPTIBLE=1, INFECTED=2, RESISTANT=3, SENESCENT=4
    custom_cmap = ListedColormap(colors)
    
    # Create legend in the same order
    legend_elements = [
        Patch(facecolor='grey', label='EmptySpot'),
        Patch(facecolor='pink', label='Susceptible'),
        Patch(facecolor='red', label='Infected'),
        Patch(facecolor='brown', label='Resistant'),
        Patch(facecolor='gold', label='Senescent'),
    ]

    # Set up the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title("Infection Spread Over Time")
    
    # Add legend
    ax.legend(handles=legend_elements, 
              loc='upper left', 
              bbox_to_anchor=(1.05, 1),
              title="Cell Types")
    
    # Set fixed aspect ratio and limits
    ax.set_aspect('equal')
    grid_size = timepoint_cell_states_matrices_array[0].shape[0]
    ax.set_xlim(-0.5, grid_size - 0.5)
    ax.set_ylim(-0.5, grid_size - 0.5)
    
    def init():
        im = ax.imshow(timepoint_cell_states_matrices_array[0], cmap=custom_cmap, 
                      vmin=0, vmax=4, origin='lower',
                      extent=[-0.5, grid_size-0.5, -0.5, grid_size-0.5])
        ax.set_title(f"Timepoint 0")
        return [im]
        
    def animate(i: int):
        ax.clear()
        # Reapply settings after clear
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(f"Timepoint {i}")
        ax.set_aspect('equal')
        ax.set_xlim(-0.5, grid_size - 0.5)
        ax.set_ylim(-0.5, grid_size - 0.5)
        
        # Add legend to each frame
        ax.legend(handles=legend_elements, 
                  loc='upper left', 
                  bbox_to_anchor=(1.05, 1),
                  title="Cell Types")
        
        im = ax.imshow(timepoint_cell_states_matrices_array[i], cmap=custom_cmap, 
                      vmin=0, vmax=4, origin='lower',
                      extent=[-0.5, grid_size-0.5, -0.5, grid_size-0.5])
        return [im]

    plt.tight_layout()
    os.makedirs("./output", exist_ok=True)
    
    anim = FuncAnimation(fig, func=animate, init_func=init, 
                        frames=len(timepoint_cell_states_matrices_array), interval=10000, 
                        blit=False, repeat=False)  # blit=False for better compatibility
    
    # Save with higher quality
    # anim.save("./output/infection_spread_animation.gif", 
    #          writer='pillow', fps=1, dpi=150)
    
    try:
        writer = FFMpegWriter(fps=1, metadata=dict(artist='Simulation'), bitrate=1800)
        filename = "./output/infection_spread_animation.mp4"
        anim.save(filename, writer=writer)
    except Exception as e:
        print(f"Failed to save MP4: {e}. Falling back to GIF.")
        writer = PillowWriter(fps=fps)
        filename = "./output/infection_spread_animation.gif"
        anim.save(filename, writer=writer)
            
    plt.close(fig)
    
def visualize_virus_surface(timepoint_virus_level_array: np.ndarray):
    """Visualizes the infection level surface over time using matrix data.
    On axes there is number of viral particles."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Define a colormap from white to red
    virus_cmap = LinearSegmentedColormap.from_list('virus_cmap', ['white', 'orange', 'red'])
    
    # Set up the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title("Virus Level Surface Over Time")
    
    global_max: int = np.max(timepoint_virus_level_array)
    
    # Create initial image and colorbar ONCE
    im = ax.imshow(timepoint_virus_level_array[0], cmap=virus_cmap, 
                  vmin=0.0, vmax=global_max, origin='lower')
    cbar = plt.colorbar(im, ax=ax, label='Quantity of viral particles')
    cbar.set_label('Quantity of viral particles', rotation=270, labelpad=15)
    
    def init():
        return [im]
        
    def animate(i: int):
        # Update only the image data, not the entire plot
        im.set_array(timepoint_virus_level_array[i])
        ax.set_title(f"Timepoint {i}")
        
        return [im]

    plt.tight_layout()
    os.makedirs("./output", exist_ok=True)
    
    anim = FuncAnimation(fig, func=animate, init_func=init, 
                        frames=len(timepoint_virus_level_array), interval=10000, 
                        blit=True, repeat=False)  # blit=True for better performance
    
    # anim.save("./output/virus_level_surface_animation.gif", 
    #          writer='pillow', fps=1, dpi=150)
    try:
        writer = FFMpegWriter(fps=1, metadata=dict(artist='Simulation'), bitrate=1800)
        filename = "./output/virus_level_surface_animation.mp4"
        anim.save(filename, writer=writer)
    except Exception as e:
        print(f"Failed to save MP4: {e}. Falling back to GIF.")
        writer = PillowWriter(fps=fps)
        filename = "./output/virus_level_surface_animation.gif"
        anim.save(filename, writer=writer)
    
    plt.close(fig)
    print("Virus level surface animation saved to ./output/virus_level_surface_animation.mp4")
    
    
def visualize_interferon_surface(timepoint_interferon_level_array: np.ndarray):
    """Visualizes the interferon level surface over time using matrix data.
    On axes there is number of interferon molecules per point."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Define a colormap from white to red
    interferon_cmap = LinearSegmentedColormap.from_list('interferon_cmap', ['white', 'orange', 'brown'])
    
    # Set up the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title("Interferon Level Surface Over Time")
    
    global_max: int = np.max(timepoint_interferon_level_array)
    
    # Create initial image and colorbar ONCE
    im = ax.imshow(timepoint_interferon_level_array[0], cmap=interferon_cmap, 
                  vmin=0.0, vmax=global_max, origin='lower')
    cbar = plt.colorbar(im, ax=ax, label='Interferon molecule number')
    cbar.set_label('Interferon molecule number', rotation=270, labelpad=15)
    
    def init():
        return [im]
        
    def animate(i: int):
        # Update only the image data, not the entire plot
        im.set_array(timepoint_interferon_level_array[i])
        
        # Update title only
        current_max = np.max(timepoint_interferon_level_array[i])
        ax.set_title(f"Timepoint {i}")
        
        return [im]

    plt.tight_layout()
    os.makedirs("./output", exist_ok=True)
    
    anim = FuncAnimation(fig, func=animate, init_func=init, 
                        frames=len(timepoint_interferon_level_array), interval=10000, 
                        blit=True, repeat=False)  # blit=True for better performance
    
    # anim.save("./output/interferon_level_surface_animation.gif", 
    #          writer='pillow', fps=1, dpi=150)
    
    try:
        writer = FFMpegWriter(fps=1, metadata=dict(artist='Simulation'), bitrate=1800)
        filename = "./output/interferon_level_surface_animation.mp4"
        anim.save(filename, writer=writer)
    except Exception as e:
        print(f"Failed to save MP4: {e}. Falling back to GIF.")
        writer = PillowWriter(fps=fps)
        filename = "./output/interferon_level_surface_animation.gif"
        anim.save(filename, writer=writer)
        
    plt.close(fig)
    print("Interferon level surface animation saved to ./output/interferon_level_surface_animation.mp4")
