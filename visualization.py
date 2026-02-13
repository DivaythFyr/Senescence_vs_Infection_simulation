import numpy as np
from matplotlib import pyplot as plt 
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
import os
import numpy.typing as npt
# from state_grid import EMPTY, SUSCEPTIBLE, INFECTED, RESISTANT, SENESCENT
# from simulation_core import MatrixSimulation

def visualize_infection_spread(timepoint_cell_states_matrices_array: npt.NDArray[np.uint8]):
    fig, ax = plt.subplots(figsize=(13, 12))  # taller figure

    colors = ['grey', 'pink', 'red', 'brown', 'gold']
    custom_cmap = ListedColormap(colors)

    legend_elements = [
        Patch(facecolor='grey',    label='EmptySpot'),
        Patch(facecolor='pink',    label='Susceptible'),
        Patch(facecolor='red',     label='Infected'),
        Patch(facecolor='brown',   label='Resistant'),
        Patch(facecolor='gold',    label='Senescent'),
    ]

    grid_size = timepoint_cell_states_matrices_array[0].shape[0]

    # IMPORTANT: reserve space at bottom for the count text
    plt.subplots_adjust(left=0.07, right=0.82, top=0.95, bottom=0.18)

    def init():
        im = ax.imshow(
            timepoint_cell_states_matrices_array[0],
            cmap=custom_cmap,
            vmin=0, vmax=4,
            origin='lower',
            extent=[-0.5, grid_size-0.5, -0.5, grid_size-0.5]
        )

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')
        ax.set_xlim(-0.5, grid_size - 0.5)
        ax.set_ylim(-0.5, grid_size - 0.5)

        ax.legend(handles=legend_elements,
                  loc='upper left', bbox_to_anchor=(1.05, 1),
                  title="Cell Types")

        ax.set_title("Timepoint 0", fontsize=14, pad=10)

        counts = np.bincount(timepoint_cell_states_matrices_array[0].ravel(), minlength=5)
        count_str = f"EmptySpot: {counts[0]:,}   Susceptible: {counts[1]:,}   Infected: {counts[2]:,}   Resistant: {counts[3]:,}   Senescent: {counts[4]:,}"
        fig.text(
            0.5, 0.04,                     
            count_str,
            ha='center', va='bottom',
            fontsize=13, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.95, pad=8, boxstyle='round,pad=0.5')
        )

        return [im]

    def animate(i: int):
        ax.clear()

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')
        ax.set_xlim(-0.5, grid_size - 0.5)
        ax.set_ylim(-0.5, grid_size - 0.5)

        ax.legend(handles=legend_elements,
                  loc='upper left', bbox_to_anchor=(1.05, 1),
                  title="Cell Types")

        im = ax.imshow(
            timepoint_cell_states_matrices_array[i],
            cmap=custom_cmap,
            vmin=0, vmax=4,
            origin='lower',
            extent=[-0.5, grid_size-0.5, -0.5, grid_size-0.5]
        )

        ax.set_title(f"Timepoint {i}", fontsize=14, pad=10)

        counts = np.bincount(timepoint_cell_states_matrices_array[i].ravel(), minlength=5)
        count_str = f"EmptySpot: {counts[0]:,}   Susceptible: {counts[1]:,}   Infected: {counts[2]:,}   Resistant: {counts[3]:,}   Senescent: {counts[4]:,}"
        fig.text(
            0.5, 0.04,
            count_str,
            ha='center', va='bottom',
            fontsize=13, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.95, pad=8, boxstyle='round,pad=0.5')
        )

        return [im]

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=len(timepoint_cell_states_matrices_array),
                         interval=2000, blit=False, repeat=False)

    os.makedirs("./output", exist_ok=True)

    try:
        writer = FFMpegWriter(fps=8, metadata=dict(artist='Simulation'), bitrate=2800)
        filename = "./output/infection_spread_with_COUNTS.mp4"
        anim.save(filename, writer=writer)
        print(f"Saved MP4: {filename}")
    except Exception as e:
        print(f"MP4 save failed: {e}")
        writer = PillowWriter(fps=8)
        filename = "./output/infection_spread_with_COUNTS.gif"
        anim.save(filename, writer=writer)
        print(f"Saved GIF: {filename}")

    plt.close(fig)
    
def visualize_virus_surface(timepoint_virus_level_array: np.ndarray):
    """Visualizes the infection level surface over time using matrix data.
    On axes there is number of viral particles."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    virus_cmap = LinearSegmentedColormap.from_list('virus_cmap', ['white', 'orange', 'red'])
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title("Virus Level Surface Over Time")
    
    global_max: int = np.max(timepoint_virus_level_array)
    
    im = ax.imshow(timepoint_virus_level_array[0], cmap=virus_cmap, 
                  vmin=0.0, vmax=global_max, origin='lower')
    cbar = plt.colorbar(im, ax=ax, label='Quantity of viral particles')
    cbar.set_label('Quantity of viral particles', rotation=270, labelpad=15)
    
    def init():
        return [im]
        
    def animate(i: int):
  
        im.set_array(timepoint_virus_level_array[i])
        ax.set_title(f"Timepoint {i}")
        
        return [im]

    plt.tight_layout()
    os.makedirs("./output", exist_ok=True)
    
    anim = FuncAnimation(fig, func=animate, init_func=init, 
                        frames=len(timepoint_virus_level_array), interval=10000, 
                        blit=True, repeat=False)  # blit=True for better performance
    
 
    try:
        writer = FFMpegWriter(fps=8, metadata=dict(artist='Simulation'), bitrate=1800)
        filename = "./output/virus_level_surface_animation.mp4"
        anim.save(filename, writer=writer)
    except Exception as e:
        print(f"Failed to save MP4: {e}. Falling back to GIF.")
        writer = PillowWriter(fps=8)
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
        writer = FFMpegWriter(fps=8, metadata=dict(artist='Simulation'), bitrate=1800)
        filename = "./output/interferon_level_surface_animation.mp4"
        anim.save(filename, writer=writer)
    except Exception as e:
        print(f"Failed to save MP4: {e}. Falling back to GIF.")
        writer = PillowWriter(fps=8)
        filename = "./output/interferon_level_surface_animation.gif"
        anim.save(filename, writer=writer)
        
    plt.close(fig)
    print("Interferon level surface animation saved to ./output/interferon_level_surface_animation.mp4")
