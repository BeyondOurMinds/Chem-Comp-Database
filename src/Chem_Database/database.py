from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
import pubchempy as pcp
from io import BytesIO
import sqlite3
import asyncio
import os
from tkinter import messagebox
import logging
import certifi
import ssl


os.environ['SSL_CERT_FILE'] = certifi.where() # certifi used to ensure that the application can verify SSL (Secure Sockets Layer) certificates when making requests to PubChem, which is necessary for fetching IUPAC names based on SMILES strings. This sets the environment variable to point to certifi's certificate bundle, allowing secure HTTPS connections. This wasn't needed when run directly, but when turned into an executable with pyinstaller, the SSL certificate verification was failing without this, likely because the bundled Python environment in the executable couldn't find the default certificates. By explicitly setting this, we ensure that SSL requests work correctly in the packaged application.
ssl._create_default_https_context = ssl.create_default_context


logging.basicConfig(
    filename="database_errors.log",
    level=logging.WARNING,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


class Database:
    def __init__(self, sdf_file, db_file=os.path.join(os.getcwd(), "chem_database.db")):
        self.sdf_file = sdf_file
        self.db_file = db_file
        self.df = None
        self.skipped_entries = []

    def load_sdf(self):
        """Load SDF file and prepare DataFrame"""
        self.skipped_entries = []
        supplier = Chem.SDMolSupplier(self.sdf_file, sanitize=False, strictParsing=False) # Disable sanitization and strict parsing to allow loading of all molecules, even those with issues. We'll handle sanitization manually to identify and skip problematic entries while still loading the rest of the data.
        
        # Check if CdId property exists in the SDF file
        has_cdid = False
        if len(supplier) > 0:
            first_mol = supplier[0]
            if first_mol is not None and first_mol.HasProp("CdId"):
                has_cdid = True
        
        for idx in range(len(supplier)):
            ''' Try to get the CdId for logging purposes, even if the molecule fails to load or sanitize. This way we can keep track of which entries are being skipped due to issues. We use CdId if available, otherwise use the entry position (idx + 1) as the identifier. This helps us maintain a record of all skipped entries for debugging and user feedback. '''
            mol = supplier[idx]
            cdid = None
            if has_cdid:
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
            
            # Use entry position as fallback identifier
            identifier = cdid if cdid else str(idx + 1)

            if mol is None:
                self.skipped_entries.append(identifier) # Log the identifier of the entry that failed to load as a molecule. This way we can keep track of which specific entries are causing issues in the SDF file.
                continue

            # Try sanitization; any non-zero result indicates a sanitization error.
            sanitize_result = Chem.SanitizeMol(mol, catchErrors=True)
            if sanitize_result != 0:
                self.skipped_entries.append(identifier)

        self.df = PandasTools.LoadSDF(self.sdf_file, smilesName="SMILES") # Load the SDF file into a DataFrame, which will include all entries. The sanitization errors will be handled separately, and we will keep track of which entries were skipped due to these errors. This allows us to load as much data as possible while still providing feedback on any issues encountered during loading and sanitization.
        
        # If CdId column doesn't exist, create it using sequential numbering
        if 'CdId' not in self.df.columns:
            self.df['CdId'] = range(1, len(self.df) + 1)
            messagebox.showinfo("CdId Column Created", "CdId column not found in SDF file. Created sequential IDs (1, 2, 3...) for database compatibility.")
        else:
            # Ensure CdId is numeric (convert from string if necessary)
            self.df['CdId'] = self.df['CdId'].astype(int)
        
        if len(self.skipped_entries) > 0:
            logging.warning(f"Skipped {len(self.skipped_entries)} entries due to sanitization errors. Identifiers: {self.skipped_entries}")
            messagebox.showwarning("Skipped Entries", f'Skipped {len(self.skipped_entries)} entries due to sanitization errors. Identifiers: {self.skipped_entries}')
        
        print(self.df)
        print(f"\nDataFrame shape: {self.df.shape}")
        
        # Convert ROMol objects to SMILES strings (which are serializable)
        if 'ROMol' in self.df.columns:
            self.df = self.df.drop('ROMol', axis=1)  # Drop the non-serializable ROMol column
        
        print(f"Columns: {self.df.columns.tolist()}")
        print(f"First 5 rows:\n{self.df.head()}")

    async def get_compound_async(self, smiles, semaphore, delay=0.3):
        """Fetch IUPAC name for a single SMILES string with rate limiting"""
        async with semaphore:  # Limit concurrent requests
            try:
                await asyncio.sleep(delay)  # Add delay between requests
                # Run the blocking PubChem call in a thread pool
                chem = await asyncio.to_thread(pcp.get_compounds, smiles, 'smiles')
                return chem[0].iupac_name if len(chem) > 0 and chem[0].iupac_name else "N/A"
            except Exception as e:
                logging.error(f"SMILES: {smiles} | Error: {e}")
                return "N/A"

    async def fetch_all_iupac_names(self):
        """Fetch IUPAC names for all molecules in the DataFrame"""
        if self.df is None or 'SMILES' not in self.df.columns:
            raise ValueError("DataFrame not loaded or SMILES column not found. Call load_sdf() first.")
        smiles_list = self.df['SMILES'].tolist()
        semaphore = asyncio.Semaphore(4)  # Only 4 concurrent requests at a time to avoid hitting PubChem rate limits
        return await asyncio.gather(*[self.get_compound_async(smiles, semaphore) for smiles in smiles_list])

    def create_database(self, filter_list=None):
        """Create SQLite database and insert all molecule data"""
        if self.df is None:
            raise ValueError("DataFrame not loaded. Call load_sdf() first.")
        
        messagebox.showinfo("Fetching IUPAC Names", "Fetching IUPAC names from PubChem (this may take a few minutes)...")
        iupac_names = asyncio.run(self.fetch_all_iupac_names())
        print("IUPAC names fetched!")
        
        # Connect to SQLite database (or create it if it doesn't exist)
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        
        # Create the molecules table
        cur.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                CdId INTEGER PRIMARY KEY,
                EntryOrder INTEGER,
                ID TEXT,
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
        
        entry_order = 0 # used to keep track of the order in which the entries are inserted into the database. This was initially used for sorting and calling purposes, but is now no longer called directly. However, it is still being inserted into the database as a column, which allows for potential future use in sorting or other operations. It also serves as a record of the original order of the entries as they were processed and inserted into the database, which could be useful for debugging or analysis purposes.
        for idx, row in enumerate(self.df.itertuples()):
            """
            Iterate over each row in the DataFrame and insert the molecule data into the database. For each molecule, we extract the necessary information such as CdId, SMILES, and IUPAC name. We then use RDKit to calculate various molecular descriptors like molecular weight, LogP, number of hydrogen bond donors and acceptors, number of rotatable bonds, and ring count. We also generate an image of the molecule and convert it to binary data for storage in the database. The filter_list parameter allows us to apply certain filters (e.g., excluding molecules with LogP >= 5) before inserting them into the database.
            """
            cd = row.CdId
            smiles = row.SMILES
            id = row.ID
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
            cur.execute("INSERT INTO molecules (CdId, EntryOrder, ID, Structure, SMILES, IUPAC_NAME, Mol_Weight, LogP, H_Bond_Donors, H_Bond_Acceptors, Rotatable_Bonds, Ring_Count) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        (cd, entry_order, id, img_data, smiles, iupac_name, molWeight, logp, bonddonor, bondacceptor, bond, ringCount))
        
        query = "SELECT COUNT(*) FROM molecules WHERE IUPAC_NAME = 'N/A'"
        missing_iupac_count = con.execute(query).fetchone()[0]
        if missing_iupac_count > 0:
            logging.warning(f"Could not fetch IUPAC names for {missing_iupac_count} molecules, marked as 'N/A'")
            messagebox.showwarning("IUPAC Name Warning", f"Could not fetch IUPAC names for {missing_iupac_count} molecules and are marked as 'N/A'. This may be due to your internet connection, rate limits, unrecognized SMILES, or other issues. Please check your internet connection and try again later.")

        con.commit()
        cur.close()
        con.close()
        print(f"\nDatabase created and {entry_order} rows inserted successfully.")
        return self.db_file
