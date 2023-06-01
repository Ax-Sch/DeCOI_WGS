 
#%%

import pandas as pd
import pandas_plink
from xarray import DataArray
import shutil
import dask.array

# %%
mask_def_file = snakemake.input.regenie_mask_def
aaf_file = snakemake.input.data_anno_aaf_file
anno_file = snakemake.input.regenie_anno_file

af_cutoffs = list(map(float, snakemake.params.af_cutoffs))
plink_in_prefix = snakemake.params.plink_in_prefix

out_prefix = snakemake.params.plink_out_folder
set_file = snakemake.output.setFile_rvtest

debug_joined_out = snakemake.output.debug_joined_out
#%%

# read plink with pandas_plink
(bim, fam, bed) = pandas_plink.read_plink(plink_in_prefix, verbose=False)
bim.cm = bim.cm.astype(int)
#%%
# Read regenie aaf_file, anno_file, mask_def
aafs = pd.read_csv(aaf_file, sep=' ', names=['ID', 'aaf'])

mask_def = pd.read_csv(mask_def_file, sep=' ', names=['mask', "annots"])
mask_def.annots = mask_def.annots.map(lambda a: a.split(','))
mask_def = mask_def.set_index("mask")

anno = pd.read_csv(anno_file, sep=' ', names=['ID', 'gene', 'annot'])
anno = anno.dropna()
# %%
# Create M0-M4 True/False columns in table
for mask, annots in mask_def.annots.items():
    anno[mask] = anno.annot.isin(annots)

# Join AF info into table
anno_aaf = anno.join(aafs.set_index("ID"), on="ID", how="outer")

# %%

# join bim (variant table) with regenie anno
joined = bim.join(anno_aaf.set_index("ID"), on="snp", how="inner")
joined = joined.sort_values(by=['gene', 'pos'])

joined.to_csv(debug_joined_out, sep="\t")

#%%

for mask in mask_def.index.unique():
    for af_cutoff in af_cutoffs:
        mask_out_prefix = f"{out_prefix}{mask}_{af_cutoff}"
        # subset variant table
        subset = joined[joined[mask] & (joined["aaf"] < af_cutoff)]
        # export genotype data (subset["i"] contains positions of variants in table)
        pandas_plink.write_plink1_bin(
            DataArray(bed[subset["i"]], dims=["variant", "sample"]),
            bed=mask_out_prefix+".bed"
            )
        # write variant subset to bim file
        subset[["gene","snp","cm","pos","a0","a1"]].to_csv(mask_out_prefix+".bim", sep=" ", header=False, index=False)
        # reuse original fam file
        shutil.copy(plink_in_prefix + ".fam", mask_out_prefix +".fam")

#%%

# create setfile, each gene has its own "chromosome"/contig
set_str = ""
for gene in joined.gene.unique():
    set_str += f"{gene}\t{gene}:0-{joined.pos.max() + 100}\n"

with open(set_file, "wt") as f:
    f.write(set_str)

