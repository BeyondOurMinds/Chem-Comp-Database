import sqlite3
from io import BytesIO
from PIL import Image
import logging

# Configure logging
logging.basicConfig(
    filename="database_errors.log",
    level=logging.WARNING,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

class ImageHandler:
    """
    Handles fetching molecule images from the SQLite database that will be displayed in the GUI. This class is separate from InfoHandler to maintain a clear separation of concerns.
    """
    def __init__(self, db_file='data/chem_database.db'):
        self.db_file = db_file

    def get_image_by_offset(self, offset, sort_column='CdId', sort_direction='ASC'):
        """Fetch image data from the database using offset and return a PIL Image object"""
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        query = f'SELECT Structure FROM molecules ORDER BY {sort_column} {sort_direction} LIMIT 1 OFFSET ?'
        cur.execute(query, (offset,))
        result = cur.fetchone()
        con.close()

        if result is None:
            error_msg = f"No image found at offset {offset} (sort: {sort_column} {sort_direction})"
            logging.error(error_msg)
            raise ValueError(error_msg)
        
        img_data = result[0]
        img = Image.open(BytesIO(img_data))
        return img

