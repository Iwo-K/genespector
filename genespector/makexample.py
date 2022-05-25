# Scripts for making sall example adata file, based on the
# LKLSK_smallexample_compressed.h5ad file

import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read("./data/LKLSK_smallexample_compressed.h5ad")

# Making the file smaller, so that it does not increase the repo size too much
sc.pp.subsample(adata, fraction=0.05, random_state=123)

# Adding a field in .obs slot with both negative and positive values for testing the diverging colourscale
adata.obs["posneg"] = np.random.randn(adata.shape[0])

# Adding specific colours in the .uns slot
adata.uns["louvain_colors"] = [
    "#569072",
    "#4f2a98",
    "#7ae26f",
    "#bd51c0",
    "#55a735",
    "#756de1",
    "#a8d358",
    "#642876",
    "#cec33d",
    "#54509c",
    "#d29933",
    "#6886d6",
    "#d34e30",
    "#5fda98",
    "#dc5196",
    "#46863d",
    "#c883d4",
    "#cfcc72",
    "#32315c",
    "#d07839",
    "#6dc5dc",
    "#d7425e",
    "#71ddca",
    "#962f68",
    "#aecc98",
    "#5d2847",
    "#7e8233",
    "#96acdf",
    "#94342d",
    "#3a879c",
    "#cf7f74",
    "#344b26",
    "#e0a4c7",
    "#7f5d30",
    "#4a658e",
    "#d4af81",
    "#652e26",
    "#a0688a",
]

adata.uns["random4_colors"] = ["#b25c4d", "#64acaf", "#8b5aa5", "#91ad58"]

adata.write("./data/tiny_example1.h5ad", compression="lzf")
