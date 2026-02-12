import sqlite3
from io import BytesIO
from PIL import Image, ImageTk

class ImageHandler:
    def __init__(self, db_file='data/chem_database.db'):
        self.db_file = db_file

    def get_image_by_cd_id(self, cd_id):
        """Fetch image data from the database using CdId and return a PIL Image object"""
        con = sqlite3.connect(self.db_file)
        query = f"SELECT Structure FROM molecules WHERE CdId = {cd_id}"
        img_data = con.execute(query).fetchone()[0]
        con.close()
        
        img = Image.open(BytesIO(img_data))
        return img

