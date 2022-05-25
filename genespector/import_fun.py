import scanpy as sc
import pandas as pd
import numpy as np
import scipy as scipy
import re


def import_adata(adata, use_raw=True, make_var_unique=True):
    """
    Function for extracting data from adata objects
    required for plotting by genespector
    Consider using backed adata to save memory?
    """

    if make_var_unique:
        adata.var_names_make_unique()

    # Extracting categorical/string data
    catcols = [str(x) in ["category"] for x in adata.obs.dtypes]
    catinfo = adata.obs.loc[:, catcols]
    # The app now takes only the categorical columns from the adata.obs slot
    # Previous code to import both and covert everything into string
    # catcols = [str(x) in ['object', 'category'] for x in adata.obs.dtypes]
    # for x in catinfo:
    #     catinfo[x] = catinfo[x].astype(str)
    # Might need a warning here for columns which have too many levels

    # Extracting coordinate data (all the obsm fields)
    coords = pd.DataFrame(index=adata.obs.index)
    for i in list(adata.obsm.keys())[::-1]:
        colnames = [i + "_" + str(x) for x in range(1, adata.obsm[i].shape[1] + 1)]
        coordsI = pd.DataFrame(adata.obsm[i], columns=colnames, index=adata.obs.index)
        coords = pd.concat((coords, coordsI), axis=1)

    # Extracting gene expression data
    if use_raw:
        genexpr = adata.raw.X
    else:
        genexpr = adata.X

    # Using scipy sparse matrix (csr) type for efficiency. If the array is not sparse,
    # it is converted
    if type(genexpr) != scipy.sparse.csr.csr_matrix:
        genexpr = scipy.sparse.csr_matrix(genexpr)
    if use_raw:
        genekeys = (
            adata.raw.var.index.values
        )  # Names of the columns in the gene expression array
    else:
        genekeys = adata.var.index.values

    # Extracting numerical data from the .obs slot
    numcols = [
        str(x) in ["float64", "float32", "in32", "int64"] for x in adata.obs.dtypes
    ]
    nummeta = adata.obs.loc[:, numcols]
    numkeys = nummeta.columns.values

    numinfo = np.empty((adata.obs.shape[0], 0))
    for i in nummeta:
        numI = adata.obs[i].astype("float64")
        numI = numI.values.reshape((numI.shape[0], 1))
        numinfo = np.concatenate((numinfo, numI), axis=1)

    # extracting colour data
    color_catinfo = {}
    for i in adata.uns.keys():
        i = str(i)
        check = re.sub("(.*)(_)(.*)", "\\3", i)
        if check == "colors":
            category = re.sub("(.*)(_)(.*)", "\\1", i)
            color_catinfo[category] = {
                x: adata.uns[i][n]
                for n, x in enumerate(adata.obs[category].cat.categories)
            }

    return {
        "catinfo": catinfo,
        "color_catinfo": color_catinfo,
        "numkeys": numkeys,
        "numinfo": numinfo,
        "genexpr": genexpr,
        "genekeys": genekeys,
        "coords": coords,
    }
    # Returns a dictionary (graphdata) with 5 dataframes/arrays
