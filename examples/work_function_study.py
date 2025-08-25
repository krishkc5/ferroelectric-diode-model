import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fed_model.device.structure import FerroeelectricDiode
from fed_model.materials.database import materials_db
from fed_model.utils.plotting import plot_work_function_comparison

def main():
    # Compare electrode materials
    electrode_materials = ['Au', 'Pt', 'Ti', 'Al', 'Cu']
    work_functions = [materials_db.get_work_function(mat) for mat in electrode_materials]
    
    fig1, ax1 = plot_work_function_comparison(electrode_materials, work_functions)
    plt.show()
    
    # Test electrode combinations
    combinations = [('Ti', 'Au'), ('Au', 'Ti'), ('Al', 'Pt'), ('Cu', 'Cu')]
    
    built_in_voltages = []
    for bottom, top in combinations:
        device = FerroeelectricDiode()
        device.build_standard_structure(
            bottom_electrode=bottom,
            top_electrode=top,
            ferroelectric='HfO2',
            ferroelectric_thickness=5,
            oxide_thickness=2
        )
        
        vbi = device.calculate_built_in_voltage()
        built_in_voltages.append(vbi)
        
        print(f"{bottom}/{top}: Built-in voltage = {vbi:+.3f} V")
    
    # Plot built-in voltages
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    combo_labels = [f"{b}/{t}" for b, t in combinations]
    colors = ['red' if v > 0 else 'blue' for v in built_in_voltages]
    
    bars = ax2.bar(combo_labels, built_in_voltages, color=colors, alpha=0.7)
    
    for bar, vbi in zip(bars, built_in_voltages):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.02*np.sign(height),
                f'{vbi:+.3f} V', ha='center', va='bottom' if height > 0 else 'top')
    
    ax2.set_ylabel('Built-in Voltage (V)')
    ax2.set_title('Built-in Voltage for Different Electrode Combinations')
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    ax2.grid(True, alpha=0.3)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
