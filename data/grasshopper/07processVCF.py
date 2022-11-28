# conda create -n allel -c bioconda -c conda-forge matplotlib pandas numpy sckikit-allel
# conda activate allel


import allel, zarr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy

zarr_path = "mitos.zarr"

allel.vcf_to_zarr("mitoVarsFreeb12dupsMarked.vcf", zarr_path, rename_fields={"NUMALT":"numalt2"}, fields="*",overwrite=True, )

callset = zarr.open_group(zarr_path, mode='r')
#print(callset.tree())



gtDf = pd.DataFrame(callset["calldata/GT"][:][:,:,0])

gtDf.columns = callset["samples"]

gtDf["POS"] = callset["variants/POS"][:]

gtDf.to_csv("genotypes.csv", "\t")

alleleCounts = np.dstack((callset["calldata/RO"][:], callset["calldata/AO"][:]))


alleleCounts[alleleCounts==-1] = 0
sumPerSampleAndSite = np.sum(alleleCounts, axis=2)
acDf = pd.DataFrame(np.vstack([alleleCounts[:,:,i] for i in range(4)]))
acDf["position"] = np.tile(callset["variants/POS"][:], 4)
acDf["allele"] =np.repeat([1,2,3,4], sumPerSampleAndSite.shape[0])
acDf.to_csv("alleleCounts.csv", "\t")
