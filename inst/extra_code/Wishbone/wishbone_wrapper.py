import wishbone
import os
import sys
import json
import pandas as pd

temp_folder = sys.argv[1]

# Load params
p = json.load(open(temp_folder + "/params.json", "r"))

# Load sample data
scdata = wishbone.wb.SCData.from_csv(temp_folder + "/counts.tsv", data_type='sc-seq', normalize=p["normalize"])
scdata.run_pca()
scdata.run_diffusion_map(knn=p["knn"], epsilon=p["epsilon"], n_diffusion_components=p["n_diffusion_components"], n_pca_components=p["n_pca_components"], markers=p["markers"])

wb = wishbone.wb.Wishbone(scdata)
wb.run_wishbone(start_cell=p["start_cell_id"], components_list=p["components_list"], num_waypoints=p["num_waypoints"], branch=p["branch"], k=p["k"])

wb.trajectory.to_json(temp_folder + "/trajectory.json")
if p["branch"]:
    wb.branch.to_json(temp_folder + "/branch.json")
else:
    pd.Series([1 for i in range(len(wb.trajectory))], index=wb.trajectory.index).to_json(temp_folder + "/branch.json")
wb.scdata.diffusion_eigenvectors.to_csv(temp_folder + "/dm.csv")
