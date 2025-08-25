# Ferroelectric Diode Model

Python package for modeling ferroelectric diodes with work function effects.

## Installation
```bash
pip install -e .
```

## Usage
```python
from fed_model.device.structure import FerroeelectricDiode

device = FerroeelectricDiode()
device.build_standard_structure(
    bottom_electrode='Ti',
    top_electrode='Au'
)

built_in_voltage = device.calculate_built_in_voltage()
device.summary()
```

## Examples
```bash
python examples/basic_tunneling.py
python examples/work_function_study.py
```
