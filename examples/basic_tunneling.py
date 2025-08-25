import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fed_model.physics.tunneling import TunnelingCalculator
from fed_model.materials.database import materials_db
from fed_model.utils.plotting import plot_potential_profile, plot_transmission_coefficient

def main():
    # Create triangular barrier
    x = np.linspace(0, 5e-9, 1000)
    barrier_height = 2.0
    barrier_width = 3e-9
    
    V = np.zeros_like(x)
    mask = (x >= 1e-9) & (x <= 4e-9)
    x_barrier = x[mask] - 1e-9
    V[mask] = barrier_height * (1 - x_barrier / barrier_width)
    
    fig1, ax1 = plot_potential_profile(x, V, "Triangular Potential Barrier")
    plt.show()
    
    # Calculate transmission
    calc = TunnelingCalculator()
    energies = np.linspace(0.1, 2.5, 50)
    transmission = [calc.transmission_wkb(x, V, E) for E in energies]
    
    fig2, ax2 = plot_transmission_coefficient(energies, transmission, barrier_height)
    plt.show()
    
    # Show some materials
    print("Available electrode materials:")
    for mat in ['Au', 'Pt', 'Ti', 'Al', 'Cu']:
        wf = materials_db.get_work_function(mat)
        print(f"{mat}: {wf:.2f} eV")

if __name__ == "__main__":
    main()
