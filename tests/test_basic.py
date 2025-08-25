import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fed_model.materials.database import MaterialDatabase
from fed_model.device.structure import FerroeelectricDiode

def test_material_database():
    db = MaterialDatabase()
    assert db.get_work_function('Au') == 5.10
    assert len(db.list_materials()) > 0

def test_device_structure():
    device = FerroeelectricDiode()
    device.build_standard_structure()
    assert len(device.layers) == 5
    assert device.total_thickness > 0

def test_built_in_voltage():
    device = FerroeelectricDiode()
    device.build_standard_structure(bottom_electrode='Ti', top_electrode='Au')
    vbi = device.calculate_built_in_voltage()
    expected = 4.33 - 5.10
    assert abs(vbi - expected) < 0.01

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
