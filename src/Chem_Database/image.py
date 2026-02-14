import sqlite3
from io import BytesIO
from PIL import Image

class ImageHandler:
    def __init__(self, db_file='data/chem_database.db'):
        self.db_file = db_file

    def get_image_by_offset(self, offset, sort_column='CdId', sort_direction='ASC'):
        """Fetch image data from the database using CdId and return a PIL Image object"""
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        query = f'SELECT Structure FROM molecules ORDER BY {sort_column} {sort_direction} LIMIT 1 OFFSET ?'
        cur.execute(query, (offset,))
        result = cur.fetchone()
        con.close()

        if result is None:
            raise ValueError(f"No image found at offset {offset}")
        
        img_data = result[0]
        img = Image.open(BytesIO(img_data))
        return img

