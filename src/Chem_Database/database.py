from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sqlite3

df = PandasTools.LoadSDF("data/Molecules11.sdf", smilesName="SMILES")

print(df)
print(f"\nDataFrame shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")

# Convert ROMol objects to SMILES strings (which are serializable)
if 'ROMol' in df.columns:
    df = df.drop('ROMol', axis=1)  # Drop the non-serializable ROMol column

# Connect to SQLite database (or create it if it doesn't exist)
con = sqlite3.connect('data/chem_database.db')
cur = con.cursor()
for row in df.itertuples():
    id = row.CdId
    smiles = row.SMILES
    mol = Chem.MolFromSmiles(smiles)
    bond = rdMolDescriptors.CalcNumRotatableBonds(mol)
    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    bonddonor = rdMolDescriptors.CalcNumHBD(mol)
    bondacceptor = rdMolDescriptors.CalcNumHBA(mol)
    ringCount = rdMolDescriptors.CalcNumRings(mol)
    molWeight = rdMolDescriptors.CalcExactMolWt(mol)

#df.to_sql(name='molecules', con=con, if_exists='replace', index=False)
con.commit()
con.close()
print(f"\nDatabase created and {len(df)} rows inserted successfully.")