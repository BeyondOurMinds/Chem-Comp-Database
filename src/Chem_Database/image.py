from rdkit import Chem
from rdkit.Chem import Draw
import sqlite3

user_input = input("Enter a CdId to visualize the molecule (e.g., 1): ")
con = sqlite3.connect('data/chem_database.db')
query = f"SELECT SMILES FROM molecules WHERE CdId = {user_input}"
smile = con.execute(query).fetchone()[0]

mol = Chem.MolFromSmiles(smile)
img = Draw.MolToImage(mol)
img.show()