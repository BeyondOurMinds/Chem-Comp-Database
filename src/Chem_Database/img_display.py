import sqlite3
import logging

# Configure logging
logging.basicConfig(
    filename="database_errors.log",
    level=logging.WARNING,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

class InfoHandler:
    """
    Handles fetching molecule information from the SQLite database that will be displayed in the GUI. This class is separate from ImageHandler to maintain a clear separation of concerns.
    """
    def __init__(self, db_file='data/chem_database.db'):
        self.db_file = db_file

    def get_info_by_offset(self, offset, sort_column='CdId', sort_direction='ASC'):
        """Fetch molecule information from the database using CdId and return as a dictionary"""
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        query = f'SELECT CdId, SMILES, IUPAC_NAME, Mol_Weight, LogP, H_Bond_Donors, H_Bond_Acceptors, Rotatable_Bonds FROM molecules ORDER BY {sort_column} {sort_direction} LIMIT 1 OFFSET ?'
        cur.execute(query, (offset,))
        result = cur.fetchone()
        con.close()

        if result is None:
            error_msg = f"No information found at offset {offset} (sort: {sort_column} {sort_direction})"
            logging.error(error_msg)
            raise ValueError(error_msg)
        
        info = {
            'CdId': result[0],
            'SMILES': result[1],
            'IUPAC_NAME': result[2],
            'Mol_Weight': result[3],
            'LogP': result[4],
            'H_Bond_Donors': result[5],
            'H_Bond_Acceptors': result[6],
            'Rotatable_Bonds': result[7]
        }
        
        return info