from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
import pubchempy as pcp
from io import BytesIO
import sqlite3
import asyncio
from tkinter import messagebox

class Database:
    def __init__(self, sdf_file, db_file='data/chem_database.db'):
        self.sdf_file = sdf_file
        self.db_file = db_file
        self.df = None
        self.skipped_entries = []

    def load_sdf(self):
        """Load SDF file and prepare DataFrame"""
        self.skipped_entries = []
        supplier = Chem.SDMolSupplier(self.sdf_file, sanitize=False, strictParsing=False)
        for idx in range(len(supplier)):
            mol = supplier[idx]
            cdid = None
            if mol is not None and mol.HasProp("CdId"):
                cdid = mol.GetProp("CdId")
            else:
                try:
                    block = supplier.GetItemText(idx)
                    if block:
                        lines = block.splitlines()
                        for i, line in enumerate(lines):
                            if line.strip() == "> <CdId>" and i + 1 < len(lines):
                                cdid = lines[i + 1].strip()
                                break
                except Exception:
                    cdid = None

            if mol is None:
                if cdid:
                    self.skipped_entries.append(cdid)
                continue

            # Try sanitization; any non-zero result indicates a sanitization error.
            sanitize_result = Chem.SanitizeMol(mol, catchErrors=True)
            if sanitize_result != 0 and cdid:
                self.skipped_entries.append(cdid)

        self.df = PandasTools.LoadSDF(self.sdf_file, smilesName="SMILES")
        messagebox.showwarning("Skipped Entries", f'Skipped {len(self.skipped_entries)} entries due to sanitization errors. CdIds: {self.skipped_entries}')
        
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

    def create_database(self, filter_list=None):
        """Create SQLite database and insert all molecule data"""
        if self.df is None:
            raise ValueError("DataFrame not loaded. Call load_sdf() first.")
        
        messagebox.showinfo("Fetching IUPAC Names", "Fetching IUPAC names from PubChem (this may take a few minutes)...")
        iupac_names = asyncio.run(self._fetch_all_iupac_names())
        print("IUPAC names fetched!")
        
        # Connect to SQLite database (or create it if it doesn't exist)
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        
        # Create the molecules table
        cur.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                CdId INTEGER PRIMARY KEY,
                EntryOrder INTEGER,
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
        
        entry_order = 0
        for idx, row in enumerate(self.df.itertuples()):
            id = row.CdId
            smiles = row.SMILES
            iupac_name = iupac_names[idx]  # Get the corresponding IUPAC name
            mol = Chem.MolFromSmiles(smiles)
            bond = rdMolDescriptors.CalcNumRotatableBonds(mol)
            logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
            if filter_list and "LogP" in filter_list and logp >= 5:
                continue
            bonddonor = rdMolDescriptors.CalcNumHBD(mol)
            if filter_list and "H_Bond_Donors" in filter_list and bonddonor >= 5:
                continue
            bondacceptor = rdMolDescriptors.CalcNumHBA(mol)
            if filter_list and "H_Bond_Acceptors" in filter_list and bondacceptor >= 10:
                continue
            ringCount = rdMolDescriptors.CalcNumRings(mol)
            molWeight = rdMolDescriptors.CalcExactMolWt(mol)
            if filter_list and "Mol_Weight" in filter_list and molWeight >= 500:
                continue
            img = Draw.MolToImage(mol)
            img_bytes = BytesIO()
            img.save(img_bytes, format='PNG')
            img_data = img_bytes.getvalue()
            entry_order += 1
            cur.execute("INSERT INTO molecules (CdId, EntryOrder, Structure, SMILES, IUPAC_NAME, Mol_Weight, LogP, H_Bond_Donors, H_Bond_Acceptors, Rotatable_Bonds, Ring_Count) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        (id, entry_order, img_data, smiles, iupac_name, molWeight, logp, bonddonor, bondacceptor, bond, ringCount))
        
        query = "SELECT COUNT(*) FROM molecules WHERE IUPAC_NAME = 'N/A'"
        missing_iupac_count = con.execute(query).fetchone()[0]
        if missing_iupac_count > 0:
            messagebox.showwarning("IUPAC Name Warning", f"Could not fetch IUPAC names for {missing_iupac_count} molecules and are marked as 'N/A'. This may be due to your internet connection, rate limits, unrecognized SMILES, or other issues. Please check your internet connection and try again later.")

        con.commit()
        cur.close()
        con.close()
        print(f"\nDatabase created and {entry_order} rows inserted successfully.")
        return self.db_file
    
    # Changing display_order to be a method of the Database class so it can access the database connection and perform queries based on the selected option from the dropdown menu in the GUI
    def display_order(self, selected_option):
        column_mapping = {
            "CdId": "CdId",
            "Molecular weight": "Mol_Weight",
            "LogP": "LogP",
            "H-bond Donors": "H_Bond_Donors",
            "H-bond Acceptors": "H_Bond_Acceptors"
        }
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        option = column_mapping.get(selected_option, "CdId")  # Default to CdId if option not found
        query = f"SELECT * FROM molecules ORDER BY {option} DESC"
        results = cur.execute(query).fetchall()
        
        con.close()
        print(f"\nResults ordered by {option}:")
