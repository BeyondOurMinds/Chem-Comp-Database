from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
import pubchempy as pcp
from io import BytesIO
import sqlite3
import asyncio

df = PandasTools.LoadSDF("data/Molecules11.sdf", smilesName="SMILES")

print(df)
print(f"\nDataFrame shape: {df.shape}")


# Convert ROMol objects to SMILES strings (which are serializable)
if 'ROMol' in df.columns:
    df = df.drop('ROMol', axis=1)  # Drop the non-serializable ROMol column

print(f"Columns: {df.columns.tolist()}")

# Fetch IUPAC names with rate limiting
async def get_compound_async(smiles, semaphore, delay=0.3):
    async with semaphore:  # Limit concurrent requests
        try:
            await asyncio.sleep(delay)  # Add delay between requests
            # Run the blocking PubChem call in a thread pool
            chem = await asyncio.to_thread(pcp.get_compounds, smiles, 'smiles') # Fetch compound info from PubChem
            return chem[0].iupac_name if chem and chem[0].iupac_name else "N/A" # Return IUPAC name or "N/A" if not found for first compound (if multiple compounds are returned, we take the first one)
        except Exception as e:
            print(f"Error fetching IUPAC name for SMILES: {smiles} - {e}")
            return "N/A"

async def fetch_all_iupac_names():
    smiles_list = df['SMILES'].tolist() # Extract SMILES from DataFrame
    semaphore = asyncio.Semaphore(4)  # Only 4 concurrent requests at a time to avoid hitting PubChem rate limits
    return await asyncio.gather(*[get_compound_async(smiles, semaphore) for smiles in smiles_list])

print("Fetching IUPAC names from PubChem (this may take a few minutes)...")
iupac_names = asyncio.run(fetch_all_iupac_names())
print("IUPAC names fetched!")

# Connect to SQLite database (or create it if it doesn't exist)
con = sqlite3.connect('data/chem_database.db')
cur = con.cursor()
# Create the molecules table
cur.execute('''
    CREATE TABLE IF NOT EXISTS molecules (
        CdId INTEGER PRIMARY KEY,
        Structure BLOB,
        SMILES TEXT,
        IUPAC_NAME TEXT,
        Mol_Weight REAL,
        LogP REAL,
        H_Bond_Donors INTEGER,
        H_Bond_Acceptors INTEGER,
        Rotatable_Bonds INTEGER,
        Ring_Count INTEGER
    );
''')
for idx, row in enumerate(df.itertuples()):
    id = row.CdId
    smiles = row.SMILES
    iupac_name = iupac_names[idx]  # Get the corresponding IUPAC name
    mol = Chem.MolFromSmiles(smiles)
    bond = rdMolDescriptors.CalcNumRotatableBonds(mol)
    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    bonddonor = rdMolDescriptors.CalcNumHBD(mol)
    bondacceptor = rdMolDescriptors.CalcNumHBA(mol)
    ringCount = rdMolDescriptors.CalcNumRings(mol)
    molWeight = rdMolDescriptors.CalcExactMolWt(mol)
    img = Draw.MolToImage(mol)
    img_bytes = BytesIO()
    img.save(img_bytes, format='PNG')
    img_data = img_bytes.getvalue()
    cur.execute("INSERT INTO molecules (CdId, Structure, SMILES, IUPAC_NAME, Mol_Weight, LogP, H_Bond_Donors, H_Bond_Acceptors, Rotatable_Bonds, Ring_Count) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (id, img_data, smiles, iupac_name, molWeight, logp, bonddonor, bondacceptor, bond, ringCount))

#df.to_sql(name='molecules', con=con, if_exists='replace', index=False)
con.commit()
cur.close()
con.close()
print(f"\nDatabase created and {len(df)} rows inserted successfully.")