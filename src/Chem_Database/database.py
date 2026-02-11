import pandas as pd
from rdkit.Chem import PandasTools

df = PandasTools.LoadSDF("data/Molecules11.sdf")

print(df)