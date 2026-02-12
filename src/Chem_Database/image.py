import sqlite3
from io import BytesIO
from PIL import Image

user_input = input("Enter a CdId to visualize the molecule (e.g., 1): ")
con = sqlite3.connect('data/chem_database.db')
query = f"SELECT Structure FROM molecules WHERE CdId = {user_input}"
img_data = con.execute(query).fetchone()[0]
con.close()

img = Image.open(BytesIO(img_data))
img.show()