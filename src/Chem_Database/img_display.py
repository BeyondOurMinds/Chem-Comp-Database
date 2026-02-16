import sqlite3

class InfoHandler:
    def __init__(self, db_file='data/chem_database.db'):
        self.db_file = db_file

    def get_info_by_offset(self, offset, sort_column='CdId', sort_direction='ASC'):
        """Fetch molecule information from the database using CdId and return as a dictionary"""
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        query = f'SELECT CdId, SMILES, IUPAC_NAME, Mol_Weight, LogP, H_Bond_Donors, H_Bond_Acceptors FROM molecules ORDER BY {sort_column} {sort_direction} LIMIT 1 OFFSET ?'
        cur.execute(query, (offset,))
        result = cur.fetchone()
        con.close()

        if result is None:
            raise ValueError(f"No information found at offset {offset}")
        
        info = {
            'CdId': result[0],
            'SMILES': result[1],
            'IUPAC_NAME': result[2],
            'Mol_Weight': result[3],
            'LogP': result[4],
            'H_Bond_Donors': result[5],
            'H_Bond_Acceptors': result[6]
        }
        
        return info