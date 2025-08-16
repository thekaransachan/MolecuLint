# Molecule Classifier

A comprehensive molecular property calculator and drug-likeness evaluator for chemical compounds using SMILES notation.

## Overview

The Molecule Classifier is a Python-based tool that analyzes chemical compounds and evaluates their drug-likeness according to multiple established pharmaceutical rules. It processes SMILES (Simplified Molecular Input Line Entry System) files and calculates various molecular properties including physicochemical descriptors, structural characteristics, and drug-likeness parameters.

## Features

### Molecular Property Calculation
- **Topological Polar Surface Area (TPSA)**: Measures molecular polarity
- **Water-octanol partition coefficient (WlogP)**: Indicates lipophilicity
- **Molecular Weight (MW)**: Total atomic mass
- **Atom Count**: Total number of atoms and heavy atoms
- **Formal Charge**: Net molecular charge
- **Heteroatom Count**: Number of non-carbon atoms
- **Carbon Count**: Number of carbon atoms
- **Molecular Formula**: Chemical formula representation
- **CSP3 Fraction**: Carbon saturation measure
- **Ring Count**: Number of cyclic structures
- **Hydrogen Bond Donors/Acceptors (HBD/HBA)**: Pharmacophore features
- **Rotatable Bonds**: Molecular flexibility indicator
- **Molar Refractivity (MR)**: Molecular volume measure
- **NH/OH Count**: Specific functional group counting

### Drug-Likeness Rule Evaluation
The tool evaluates compounds against five major drug-likeness rules:

1. **Lipinski's Rule of Five**: Classical drug-likeness criteria
2. **Ghose Filter**: Extended molecular property ranges
3. **Veber Rule**: Oral bioavailability indicators
4. **Egan Rule**: Enhanced bioavailability criteria
5. **Muegge Rule**: Comprehensive drug-likeness assessment

## Files

- **`FastPredictor.py`**: Main molecular property calculator and rule evaluator
- **`PropertiesParser.py`**: Utility to convert property output to CSV format
- **`FDA_DB_SMILES.smi`**: FDA database compounds in SMILES format
- **`MoleculesLessThan200Chars.smi`**: Filtered compounds with SMILES < 200 characters
- **`new_properties.txt`**: Output file containing calculated properties

## Requirements

### Dependencies
```bash
pip install rdkit-pypi
```

### System Requirements
- Python 3.6+
- RDKit cheminformatics library
- Sufficient memory for large molecule databases

## Usage

### Basic Usage
```bash
python FastPredictor.py input_file.smi
```

### Advanced Usage
```bash
python FastPredictor.py input_file.smi -o output_file.txt
```

### Input Format
The input file should contain SMILES strings, optionally with compound names:
```
CC(C)CO ethanol
CCO methanol
```

### Output Format
The tool generates a detailed report for each compound:
```
NAME: ethanol
TPSA: 20.23
WlogP: -0.31
Atoms: 9
FormalCharge: 0
Heteroatoms: 1
Carbon: 2
Formula: C2H6O
MW: 46.07
HeavyAtoms: 3
CSP3: 1.0
Rings: 0
HBD: 1
HBA: 1
RotBonds: 1
MR: 12.71
NHOH: 1
NO: 1

Results for ethanol:
Lipinski Rules:
        No violations
Ghose Rules:
        No violations
Veber Rules:
        No violations
Egan Rules:
        No violations
Muegge Rules:
        No violations
```

### Converting Output to CSV
Use the PropertiesParser utility to convert the output to CSV format:
```python
from PropertiesParser import parse_records
parse_records('new_properties.txt', 'output.csv')
```

## Drug-Likeness Rules

### Lipinski's Rule of Five
- Molecular weight < 500 Da
- WlogP < 5
- NH/OH groups < 5
- N/O atoms < 10

### Ghose Filter
- Molecular weight: 160-480 Da
- WlogP: -0.4 to 5.6
- Molar refractivity: 40-130
- Atom count: 20-70

### Veber Rule
- Rotatable bonds < 10
- TPSA < 140 Å²

### Egan Rule
- WlogP < 5.88
- TPSA < 131.6 Å²

### Muegge Rule
- Molecular weight: 200-600 Da
- WlogP: -2 to 5
- TPSA < 150 Å²
- Rings < 7
- Carbon atoms ≥ 4
- Heteroatoms ≥ 1
- Rotatable bonds < 15
- HBA < 10

## Performance Features

- **Efficient Processing**: Single-pass file handling for large datasets
- **Memory Optimization**: Minimal memory footprint during processing
- **Error Handling**: Graceful handling of invalid SMILES strings
- **Batch Processing**: Processes entire SMILES files in one execution

## Applications

- **Drug Discovery**: Screening compound libraries for drug-likeness
- **Chemical Informatics**: Molecular property analysis and comparison
- **Pharmaceutical Research**: Lead compound optimization
- **Academic Research**: Chemical structure-property relationship studies
- **Quality Control**: Compound library validation

## Contributing

This project is designed for molecular property analysis and drug-likeness evaluation. Contributions for additional molecular descriptors, rule sets, or performance optimizations are welcome.

## License

This project is provided as-is for research and educational purposes.

## Citation

If you use this tool in your research, please cite the RDKit library and relevant drug-likeness rule publications.

## Support

For issues or questions regarding molecular property calculations or drug-likeness evaluation, please refer to the RDKit documentation or relevant cheminformatics literature.
