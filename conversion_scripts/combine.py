import os
import re
import pandas as pd
 
dir_path = "C:\\capstone\\dsm_tiles\\DSM_CGY_5x5km_2mResolution_North"
 
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
 
outfile = open("DTMCombinedNorth2_2m.bin", "wb")
 
resolution = 2
base_size = 5000
 
pixel_length = base_size / resolution
print(pixel_length)
input()

num_bytes_read = 0
north_count =0
while not df.empty:
    df_max_north = df[df["North"] == df["North"].max()]
    df_max_north = df_max_north.sort_values(by="East", ascending=True)
    north_count += 1
    print(north_count)

 
    starting_index = 0
    for starting_index in range(0, int(pixel_length*pixel_length), int(pixel_length)):
        ending_index = starting_index + pixel_length - 1
        for _, row in df_max_north.iterrows():
           
            file_path = dir_path + "\\" + row["filename"]
            print("Calling index in range: " + str(starting_index) + " " + str(ending_index))
            print(row)
 
            with open(file_path, "rb") as src:
                # Seek to the first uint16
                src.seek(starting_index * 2)
 
                # Read raw bytes
                data = src.read(int(pixel_length) * 2)
                num_bytes_read += (int(pixel_length) * 2)

                # Write to destination
                outfile.write(data)
    df = df.drop(df_max_north.index)

    print(num_bytes_read)
 
 