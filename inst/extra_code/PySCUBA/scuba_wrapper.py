import PySCUBA
import json
import sys

temp_folder = sys.argv[1]
log_mode = bool(sys.argv[2])
rigorous_gap_stats = bool(sys.argv[3])

cell_IDs, data, markers, cell_stages, data_tag, output_directory = PySCUBA.Preprocessing.RNASeq_preprocess(temp_folder + "/counts.tsv", pseudotime_mode=True)
centroid_coordinates, cluster_indices, parent_clusters = PySCUBA.initialize_tree(data, cell_stages)
centroid_coordinates, cluster_indices, parent_clusters = PySCUBA.refine_tree(data, centroid_coordinates, cluster_indices, parent_clusters, cell_stages)

json.dump({
        "tree":parent_clusters,
        "labels":cluster_indices.tolist()
    },
          open(temp_folder + "/output.json", "w")
 )
