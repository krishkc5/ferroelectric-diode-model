import numpy as np
from ..materials.database import materials_db

class Layer:
    def __init__(self, material, thickness, name=None):
        self.material = material
        self.thickness = thickness * 1e-9  # nm to m
        self.name = name or material
        
        self.work_function = materials_db.get_work_function(material)
        self.electron_affinity = materials_db.get_electron_affinity(material)
        self.dielectric_constant = materials_db.get_dielectric_constant(material)
    
    def __repr__(self):
        return f"Layer({self.name}, {self.thickness*1e9:.1f}nm)"

class FerroeelectricDiode:
    def __init__(self):
        self.layers = []
        self.total_thickness = 0
        
    def add_layer(self, material, thickness, name=None):
        layer = Layer(material, thickness, name)
        self.layers.append(layer)
        if layer.thickness > 0:
            self.total_thickness += layer.thickness
        return layer
    
    def build_standard_structure(self, 
                                bottom_electrode='Ti',
                                ferroelectric='HfO2', 
                                ferroelectric_thickness=10,
                                dead_layer_thickness=1,
                                oxide='Al2O3',
                                oxide_thickness=2,
                                top_electrode='Au'):
        
        self.layers = []
        self.total_thickness = 0
        
        self.add_layer(bottom_electrode, 0, "Bottom Electrode")
        
        if dead_layer_thickness > 0:
            self.add_layer(ferroelectric, dead_layer_thickness, "Dead Layer")
            
        self.add_layer(ferroelectric, ferroelectric_thickness, "Ferroelectric")
        
        if oxide_thickness > 0:
            self.add_layer(oxide, oxide_thickness, "Insulating Layer")
            
        self.add_layer(top_electrode, 0, "Top Electrode")
        
        return self
    
    def get_position_grid(self, num_points=1000):
        active_layers = [l for l in self.layers if l.thickness > 0]
        
        positions = []
        current_pos = 0
        
        for layer in active_layers:
            layer_points = int(num_points * layer.thickness / self.total_thickness)
            layer_positions = np.linspace(current_pos, current_pos + layer.thickness, layer_points)
            positions.extend(layer_positions[:-1])
            current_pos += layer.thickness
        
        positions.append(current_pos)
        return np.array(positions)
    
    def get_layer_at_position(self, x):
        current_pos = 0
        for layer in self.layers:
            if layer.thickness == 0:
                continue
            if current_pos <= x <= current_pos + layer.thickness:
                return layer
            current_pos += layer.thickness
        return None
    
    def calculate_built_in_voltage(self):
        electrodes = [l for l in self.layers if 'Electrode' in l.name]
        
        if len(electrodes) != 2 or not all(e.work_function for e in electrodes):
            return 0.0
            
        return electrodes[0].work_function - electrodes[1].work_function
    
    def summary(self):
        print("Ferroelectric Diode Structure")
        print(f"Total active thickness: {self.total_thickness*1e9:.2f} nm")
        print(f"Built-in voltage: {self.calculate_built_in_voltage():.3f} V")
        print("\nLayers:")
        
        for i, layer in enumerate(self.layers):
            thickness_str = f"{layer.thickness*1e9:.1f} nm" if layer.thickness > 0 else "electrode"
            wf_str = f"{layer.work_function:.2f} eV" if layer.work_function else "N/A"
            print(f"  {layer.name} ({layer.material}): {thickness_str}, WF: {wf_str}")
