import pandas as pd
import numpy as np

# Work functions in eV - modify these as needed
WORK_FUNCTIONS = {
    'Au': 5.10, 'Pt': 5.65, 'Pd': 5.12, 'Ag': 4.26, 'Cu': 4.65,
    'Al': 4.28, 'Ti': 4.33, 'Cr': 4.50, 'Ni': 5.15, 'W': 4.55,
    'Mo': 4.36, 'Ta': 4.25, 'Si': 4.85, 'Ge': 5.00,
}

# Electron affinities in eV
ELECTRON_AFFINITIES = {
    'BaTiO3': 4.0, 'PbTiO3': 3.9, 'PZT': 4.1, 'HfO2': 2.4, 'BFO': 4.2,
    'SiO2': 0.9, 'Al2O3': 2.8, 'TiO2': 4.2,
}

# Dielectric constants
DIELECTRIC_CONSTANTS = {
    'BaTiO3': 1700, 'PbTiO3': 200, 'PZT': 1000, 'HfO2': 25, 'BFO': 100,
    'SiO2': 3.9, 'Al2O3': 9.0, 'TiO2': 80,
    'Au': 1.0, 'Pt': 1.0, 'Al': 1.0, 'Ti': 1.0, 'Cu': 1.0,
}

class MaterialDatabase:
    def __init__(self):
        self.work_functions = WORK_FUNCTIONS
        self.electron_affinities = ELECTRON_AFFINITIES  
        self.dielectric_constants = DIELECTRIC_CONSTANTS
        
    def get_work_function(self, material):
        return self.work_functions.get(material)
    
    def get_electron_affinity(self, material):
        return self.electron_affinities.get(material)
        
    def get_dielectric_constant(self, material):
        return self.dielectric_constants.get(material)
    
    def list_materials(self):
        all_materials = set()
        all_materials.update(self.work_functions.keys())
        all_materials.update(self.electron_affinities.keys())
        all_materials.update(self.dielectric_constants.keys())
        return sorted(list(all_materials))
    
    def to_dataframe(self):
        materials = self.list_materials()
        data = []
        
        for material in materials:
            row = {'Material': material}
            row['Work_Function_eV'] = self.work_functions.get(material, np.nan)
            row['Electron_Affinity_eV'] = self.electron_affinities.get(material, np.nan)
            row['Dielectric_Constant'] = self.dielectric_constants.get(material, np.nan)
            data.append(row)
        
        return pd.DataFrame(data)

materials_db = MaterialDatabase()
