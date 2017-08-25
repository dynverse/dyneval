import pandas as pd
import numpy as np

from GPfates import GPfates

import sys

temp_folder = sys.argv[1]
log_expression_cutoff = float(sys.argv[2])
min_cells_expression_cutoff = float(sys.argv[3])
nfates = int(sys.argv[4])

etpm = pd.read_table(temp_folder + 'expression.csv', index_col=0)
etpm = etpm[(etpm > log_expression_cutoff).sum(1) >min_cells_expression_cutoff]
logexp = np.log10(etpm + 1)

cellinfo = pd.read_table(temp_folder + 'cellinfo.csv', index_col=0)

m = GPfates.GPfates(cellinfo, logexp)

print("Dimensionality reduction--------------------------------------")
m.dimensionality_reduction()

print("Story DR------------------------------------------------------")
m.store_dr() # store the dr in the sample table (m.s), so it can be used in the gplvm

print("Infer pseudotime----------------------------------------------")
m.infer_pseudotime(s_columns=['bgplvm_0', 'bgplvm_1']) # use the first two components to infer pseudotime

print("Model cell fates----------------------------------------------")
m.model_fates(X=['bgplvm_1'], C=nfates)

print("Saving--------------------------------------------------------")
m.s.pseudotime.to_csv(temp_folder + "pseudotimes.csv")
pd.DataFrame(m.fate_model.phi, index=m.s.pseudotime.index).to_csv(temp_folder + "phi.csv")

print("Finished------------------------------------------------------")
