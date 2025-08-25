import matplotlib.pyplot as plt
import numpy as np

plt.style.use('default')

def plot_potential_profile(x, V, title="Potential Profile"):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x * 1e9, V, 'b-', linewidth=2)
    ax.fill_between(x * 1e9, V, alpha=0.3)
    ax.set_xlabel('Position (nm)')
    ax.set_ylabel('Potential Energy (eV)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    return fig, ax

def plot_iv_curve(voltages, currents, title="I-V Characteristic"):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(voltages, np.abs(currents), 'r-', linewidth=2, marker='o')
    ax.set_xlabel('Voltage (V)')
    ax.set_ylabel('|Current| (A)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    return fig, ax

def plot_transmission_coefficient(energies, transmission, barrier_height=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(energies, transmission, 'g-', linewidth=2, marker='o')
    
    if barrier_height:
        ax.axvline(barrier_height, color='red', linestyle='--', 
                  label=f'Barrier height: {barrier_height:.2f} eV')
    
    ax.set_xlabel('Electron Energy (eV)')
    ax.set_ylabel('Transmission Coefficient')
    ax.set_title('Tunneling Probability vs Energy')
    ax.grid(True, alpha=0.3)
    if barrier_height:
        ax.legend()
    return fig, ax

def plot_work_function_comparison(materials, work_functions):
    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(materials, work_functions, color='skyblue', edgecolor='navy')
    
    for bar, wf in zip(bars, work_functions):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{wf:.2f}', ha='center', va='bottom', fontweight='bold')
    
    ax.set_ylabel('Work Function (eV)')
    ax.set_title('Work Function Comparison')
    ax.set_ylim(0, max(work_functions) * 1.2)
    plt.xticks(rotation=45)
    plt.tight_layout()
    return fig, ax
