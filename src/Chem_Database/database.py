import pandas as pd
from rdkit.Chem import PandasTools
from rdkit import Chem
import sqlite3

df = PandasTools.LoadSDF("data/Molecules11.sdf", smilesName="SMILES")

print(df)
print(f"\nDataFrame shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")

# Convert ROMol objects to SMILES strings (which are serializable)
if 'ROMol' in df.columns:
    df = df.drop('ROMol', axis=1)  # Drop the non-serializable ROMol column

# Connect to SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect('data/chem_database.db')
c = conn.cursor()
df.to_sql(name='molecules', con=conn, if_exists='replace', index=False)
conn.commit()
conn.close()
print(f"\nDatabase created and {len(df)} rows inserted successfully.")