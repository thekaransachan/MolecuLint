import os
import argparse
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors, rdqueries

def calculate_properties(mol, molH=None):
    """Calculate all molecular properties at once using the appropriate molecule object"""
    if molH is None:
        molH = Chem.AddHs(mol)
        
    return {
        'TPSA': round(Descriptors.TPSA(mol, includeSandP=True), 2),
        'WlogP': round(Descriptors.MolLogP(molH), 2),
        'Atoms': molH.GetNumAtoms(),
        'FormalCharge': Chem.GetFormalCharge(mol),
        'Heteroatoms': rdMolDescriptors.CalcNumHeteroatoms(mol),
        'Carbon': len(mol.GetAtomsMatchingQuery(rdqueries.AtomNumEqualsQueryAtom(6))),
        'Formula': rdMolDescriptors.CalcMolFormula(mol),
        'MW': round(Descriptors.ExactMolWt(mol), 2),
        'HeavyAtoms': Lipinski.HeavyAtomCount(mol),
        'CSP3': round(Lipinski.FractionCSP3(mol), 2),
        'Rings': Lipinski.RingCount(mol),
        'HBD': Lipinski.NumHDonors(mol),
        'HBA': Lipinski.NumHAcceptors(mol),
        'RotBonds': Lipinski.NumRotatableBonds(mol),
        'MR': round(Crippen.MolMR(mol), 2),
        'NHOH': Lipinski.NHOHCount(mol),
        'NO': Lipinski.NOCount(mol)
    }

def evaluate_rules(mol_props):
    """Evaluate drug-likeness rules without unused parameters"""
    results = {}

    # Lipinski
    lipinski = []
    if mol_props['MW'] >= 500: lipinski.append("Molecular weight violation: MW > 500")
    if mol_props['WlogP'] >= 5: lipinski.append("WlogP violation: WlogP > 5")
    if mol_props['NHOH'] >= 5: lipinski.append("NH or OH violation: NH or OH > 5")
    if mol_props['NO'] >= 10: lipinski.append("N or O violation: N or O > 10")
    results['Lipinski'] = lipinski

    # Ghose
    ghose = []
    if mol_props['MW'] < 160: ghose.append("Molecular weight violation: MW < 160")
    if mol_props['MW'] > 480: ghose.append("Molecular weight violation: MW > 480")
    if mol_props['WlogP'] < -0.4: ghose.append("WlogP violation: WlogP < -0.4")
    if mol_props['WlogP'] > 5.6: ghose.append("WlogP violation: WlogP > 5.6")
    if mol_props['MR'] < 40: ghose.append("Molar Refractivity violation: MR < 40")
    if mol_props['MR'] > 130: ghose.append("Molar Refractivity violation: MR > 130")
    if mol_props['Atoms'] < 20: ghose.append("Atom No. violation: Atoms < 20")
    if mol_props['Atoms'] > 70: ghose.append("Atom No. violation: Atoms > 70")
    results['Ghose'] = ghose

    # Veber
    veber = []
    if mol_props['RotBonds'] >= 10: veber.append("Rotatable bonds violation: R.Bonds > 10")
    if mol_props['TPSA'] >= 140: veber.append("TPSA violation: TPSA > 140")
    results['Veber'] = veber

    # Egan
    egan = []
    if mol_props['WlogP'] >= 5.88: egan.append("WlogP violation: WlogP > 5.88")
    if mol_props['TPSA'] >= 131.6: egan.append("TPSA violation: TPSA > 131.6")
    results['Egan'] = egan

    # Muegge
    muegge = []
    if mol_props['MW'] < 200: muegge.append("Molecular weight violation: MW < 200")
    if mol_props['MW'] > 600: muegge.append("Molecular weight violation: MW > 600")
    if mol_props['WlogP'] < -2: muegge.append("WlogP violation: WlogP < -2")
    if mol_props['WlogP'] > 5: muegge.append("WlogP violation: WlogP > 5")
    if mol_props['TPSA'] > 150: muegge.append("TPSA violation: TPSA > 150")
    if mol_props['Rings'] >= 7: muegge.append("No. of Rings violation: Rings > 7")
    if mol_props['Carbon'] < 4: muegge.append("No. of Carbon violation: C < 4")
    if mol_props['Heteroatoms'] < 1: muegge.append("No. of Heteroatoms violation: Het. Atoms < 1")
    if mol_props['RotBonds'] >= 15: muegge.append("Rotatable Bonds violation: > 15")
    if mol_props['HBA'] >= 10: muegge.append("HBA violation: > 10")
    results['Muegge'] = muegge

    return results

def main():
    parser = argparse.ArgumentParser(description="Process a SMILES file and calculate molecular properties.")
    parser.add_argument("smi_file", help="Path to the SMILES input file")
    parser.add_argument("-o", "--output", default="new_properties.txt", help="Output file name")
    args = parser.parse_args()

    if not os.path.isfile(args.smi_file):
        print(f"Error: file not found: {args.smi_file}")
        sys.exit(1)

    # Clear the output file once at the beginning
    with open(args.output, "w") as out:
        pass

    # Process all molecules with a single file handle
    with open(args.output, "a") as out, open(args.smi_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line: 
                continue

            try:
                parts = line.split()
                smiles = parts[0]
                name = parts[1] if len(parts) > 1 else f"Compound_{line_num}"

                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"Skipping invalid SMILES: {smiles}")
                    continue

                # Add hydrogens once
                molH = Chem.AddHs(mol)
                
                # Calculate properties using the appropriate molecule object
                props = calculate_properties(mol, molH)
                props['Name'] = name

                # Write properties
                out.write(f"\nNAME: {props['Name']}\n")
                for k, v in props.items():
                    if k != "Name":
                        out.write(f"{k}: {v}\n")

                # Check and display rules
                results = evaluate_rules(props)
                print(f"\nResults for {name}:")
                for rule, violations in results.items():
                    print(f"{rule} Rules:")
                    if not violations:
                        print("\tNo violations")
                    else:
                        print("\t" + "\n\t".join(violations))
                        
            except Exception as e:
                print(f"Error processing line {line_num}: {e}")

if __name__ == "__main__":
    main()
