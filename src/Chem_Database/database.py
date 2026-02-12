from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
import pubchempy as pcp
from io import BytesIO
import sqlite3
import asyncio

class Database:
    def __init__(self, sdf_file, db_file='data/chem_database.db'):
        self.sdf_file = sdf_file
        self.db_file = db_file
        self.df = None

    def load_sdf(self):
        """Load SDF file and prepare DataFrame"""
        self.df = PandasTools.LoadSDF(self.sdf_file, smilesName="SMILES")
        
        print(self.df)
        print(f"\nDataFrame shape: {self.df.shape}")
        
        # Convert ROMol objects to SMILES strings (which are serializable)
        if 'ROMol' in self.df.columns:
            self.df = self.df.drop('ROMol', axis=1)  # Drop the non-serializable ROMol column
        
        print(f"Columns: {self.df.columns.tolist()}")

    async def _get_compound_async(self, smiles, semaphore, delay=0.3):
        """Fetch IUPAC name for a single SMILES string with rate limiting"""
        async with semaphore:  # Limit concurrent requests
            try:
                await asyncio.sleep(delay)  # Add delay between requests
                # Run the blocking PubChem call in a thread pool
                chem = await asyncio.to_thread(pcp.get_compounds, smiles, 'smiles')
                return chem[0].iupac_name if chem and chem[0].iupac_name else "N/A"
            except Exception as e:
                print(f"Error fetching IUPAC name for SMILES: {smiles} - {e}")
                return "N/A"

    async def _fetch_all_iupac_names(self):
        """Fetch IUPAC names for all molecules in the DataFrame"""
        smiles_list = self.df['SMILES'].tolist()
        semaphore = asyncio.Semaphore(4)  # Only 4 concurrent requests at a time to avoid hitting PubChem rate limits
        return await asyncio.gather(*[self._get_compound_async(smiles, semaphore) for smiles in smiles_list])

    def create_database(self):
        """Create SQLite database and insert all molecule data"""
        if self.df is None:
            raise ValueError("DataFrame not loaded. Call load_sdf() first.")
        
        print("Fetching IUPAC names from PubChem (this may take a few minutes)...")
        iupac_names = asyncio.run(self._fetch_all_iupac_names())
        print("IUPAC names fetched!")
        
        # Connect to SQLite database (or create it if it doesn't exist)
        con = sqlite3.connect(self.db_file)
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
        
        for idx, row in enumerate(self.df.itertuples()):
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
        
        con.commit()
        cur.close()
        con.close()
        print(f"\nDatabase created and {len(self.df)} rows inserted successfully.")
        return self.db_file