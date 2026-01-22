import os
import re
import pandas as pd

dir_path = "C:\\capstone\\dsm_tiles\\DSM_CGY_5x5km_Res1m_SnappedProperly"

pattern = re.compile(r"_(-?\d+(?:\.\d+)?)E_(-?\d+(?:\.\d+)?)N\.bin$")

rows = []

for fname in os.listdir(dir_path):
    match = pattern.search(fname)
    if match:
        east = float(match.group(1))
        north = float(match.group(2))
        rows.append({
            "filename": fname,
            "East": east,
            "North": north
        })


df = pd.DataFrame(rows)

print(df)

outfile = open("DTMCombinedNorth_5672949.28_-15989.47_5N_5E.bin", "wb")


while not df.empty:
    df_max_north = df[df["North"] == df["North"].max()]
    df_max_north = df_max_north.sort_values(by="East", ascending=True)


    starting_index = 0
    for starting_index in range(0, 25000000, 5000):
        ending_index = starting_index + 4999
        for _, row in df_max_north.iterrows(): 
            
            file_path = dir_path + "\\" + row["filename"]
            print("Calling index in range: " + str(starting_index) + " " + str(ending_index))
            print(row)

            with open(file_path, "rb") as src:
                # Seek to the first uint16
                src.seek(starting_index * 2)

                # Read raw bytes
                data = src.read(10000)

                # Write to destination
                outfile.write(data)
    df = df.drop(df_max_north.index)
