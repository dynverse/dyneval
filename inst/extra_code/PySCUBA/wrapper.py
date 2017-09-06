import PySCUBA
import json
import sys
import numpy

temp_folder = sys.argv[1]
rigorous_gap_stats = bool(sys.argv[3].capitalize())
N_dim = int(sys.argv[4])
low_gene_threshold = float(sys.argv[5])
low_gene_fraction_max = float(sys.argv[6])
min_split = int(sys.argv[7])
min_percentage_split = float(sys.argv[8])


cell_IDs, data, markers, cell_stages, data_tag, output_directory = PySCUBA.Preprocessing.RNASeq_preprocess(temp_folder + "/counts.tsv", pseudotime_mode=True, log_mode=False, N_dim=N_dim, low_gene_threshold=low_gene_threshold, low_gene_fraction_max=low_gene_fraction_max)
centroid_coordinates, cluster_indices, parent_clusters = PySCUBA.initialize_tree(data, cell_stages, rigorous_gap_stats=rigorous_gap_stats, min_split=min_split, min_percentage_split=min_percentage_split)
centroid_coordinates, cluster_indices, parent_clusters, new_tree = PySCUBA.refine_tree(data, centroid_coordinates, cluster_indices, parent_clusters, cell_stages, output_directory=temp_folder)

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)


json.dump({
        "tree":parent_clusters,
        "labels":cluster_indices.tolist(),
        "new_tree":new_tree
    },
          open(temp_folder + "/output.json", "w")
, cls=MyEncoder
 )
