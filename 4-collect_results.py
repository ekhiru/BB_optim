#!/usr/bin/env python3
import os
import os.path
from glob import glob
import re
import pandas as pd
import sys

if len(sys.argv) < 2:
    res_dir = "./results/"
else:
    res_dir = sys.argv[1]

if not os.path.exists(res_dir):
    sys.stderr.write(f'ERROR: directory {res_dir} was not found!')
    sys.exit(1)

# traverse root directory, and list directories as dirs and files as files
for root, dirs, files in os.walk(res_dir):
    print(root)
    unique_files = []
    for file in files:
        result = re.match(r"(.+)-r([0-9]+)\.csv\.xz", file)
        if result:
            unique_files.append(result.group(1))
    unique_files = set(unique_files)
    for file in unique_files:
        runs = []
        run_data = []
        for run_file in sorted(glob(root + os.sep + file + "-r*.csv.xz")):
            result = re.match(r".+-r([0-9]+)\.csv\.xz", run_file)
            runs.append(result.group(1))
            try:
                run_data.append(pd.read_csv(run_file))
            except:
                print(f"Error: reading {run_file}")
                raise

        res_file = root + os.sep + file + ".csv.xz"
        if os.path.isfile(res_file):
            df = pd.read_csv(res_file)
            # Remove old data
            common_runs = np.intersect1d(df["seed"].unique(),runs)
            print(f'Replacing {len(common_runs)} runs {common_runs} ')
            df = df[~df["seed"].isin(runs)]
            run_data.insert(0, df)
        df = pd.concat(run_data, ignore_index = True)
        print(f"Writing runs {runs} in {res_file}")
        df.to_csv(res_file, index=False, compression = "xz")
        for run_file in glob(root + os.sep + file + "-r*.csv.xz"):
            print(f"Deleting {run_file}")
            os.remove(run_file)
