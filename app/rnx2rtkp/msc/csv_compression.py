import numpy as np
import pandas as pd

# df = pd.read_csv("C:/capstone/3d_model_gnss_capstone/app/rnx2rtkp/msc/CombinedRaster1_Complete.csv", header=None)

# df = pd.read_csv("C:---PUT-YOUR-CSV-FILE-PATH-HERE---/CombinedRaster1_Complete.csv", header=None)


uint16_data = df.astype(np.uint16)
uint16_data.to_numpy().tofile("data.bin")