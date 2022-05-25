.. role:: pyth(code)
  :language: python

Genespector
===========

Genespector is a Dash/Plotly app for interactive visualisation of
large scRNA-Seq datasets analysed using the scanpy/AnnData packages. Genespector allows responsive plotting of datasets up to 100,000 cells.
It remains stable at least up to 250,000 cells but there is a visible delay in plotting.

Genespector is based on the Dash/Plotly framework: https://dash.plot.ly/


Installation
------------

Python > 3.4 and pip are required. To install the package:

1. Clone the repository:

.. code-block:: text

    git clone https://github.com/Iwo-K/genespector

2. Install the dependencies

.. code-block:: text

    pip install -r ./genespector/requirements.txt

3. Install the package

.. code-block:: text

    pip install -e ./genespector/

Examples
--------

A simple example of scRNA-Seq dataset displayed with Genespector can be found here: http://128.232.224.252/example/

For anyone interested in haematopoiesis, the Gottgens lab is hosting several Genespector apps displaying following scRNA-Seq datasets:

1. Mouse haematopoietic stem and progenitor cells:
  a. http://128.232.224.252/sfdata/ - from the paper: Nestorowa et al. Blood 2016, PMID: 27365425
  b. http://128.232.224.252/LKdata/ - from the paper: Dahlin et al. Blood 2018, PMID: 29588278
2. Human Bone Marrow cells (from Human Cell Atlas): http://128.232.224.252/HCABM/


AnnData objects pre-processing
------------------

Genespector simply requires an AnnData object to run. Most adata object should work straight away, but you may want to select features to be shown in the app by removing unnecessary information from the AnnData object. An example jupyter notebook with data preprocessing can be found here: https://nbviewer.jupyter.org/github/Iwo-K/genespector/blob/master/example_genespector.ipynb

The adata objects are imported as follows:
  - Gene expression values are imported from the :pyth:`adata.X` (or :pyth:`adata.raw.X`, depending on the :pyth:`use_raw` argument)
  - Gene names are imported from the :pyth:`adata.var.index` (or :pyth:`adata.raw.var.index`, depending on the :pyth:`use_raw` argument) and become available from the **Color points** menu.
  - Numeric values present in the :pyth:`adata.obs` slot are imported and become available from the **Color points** menu.
  - Categorical data present in the :pyth:`adata.obs` slot are imported and become available from the **Color points** menu. Plots are rendered in layers according to the selected value in **Split data by** menu. Note: **only categorical columns are imported, columns containing strings are ignored**. To change the column type use the method: :pyth:`.astype('category')`
  - All available coordinates are imported from the :pyth:`adata.obs` slot and become available in the **Choose the X and Y axis of the plot** menu
  - User-defined color scales for categorical data are imported from the :pyth:`adata.uns` slot, as long as the names are matching. For instance :pyth:`adata.obs['louvain']` matches the entry :pyth:`adata.uns['louvain_colors']`.

When dealing with large datasets, to save memory try using sparse matrices with normalised/log-transformed data and avoid using full matrices containing scaled data.


App initialisation
------------------

In python, where adata is the AnnData object of choice:

.. code-block:: text

    import genespector as gp
    gp.make_app(adata)

Example output:

.. code-block:: text

    .........
    .........
    * Serving Flask app "app" (lazy loading)
    .........
    * Running on http://127.0.0.1:8050/ (Press CTRL+C to quit)

Copy the address to your browser or click the link.

make_app() accepts the following arguments:
  - :pyth:`adata` - an AnnData object
  - :pyth:`use_raw` - logical, whether .raw.X slot should be used instead of .X (default: False)
  - :pyth:`make_var_unique` - logical, whether .var.index should be converted to unique values using the var_names_make_unique() function from scanpy
  - :pyth:`main_title` and :pyth:`subtitle`  - string, titles displayed above the App
  - :pyth:`layout` - function, creating a dictionary which controls the website layout, needs to contain necessary elements for the app. For an example see the layout.py file
  - :pyth:`server` - logical or name of the server used for deploying the app
  - :pyth:`url_base_pathname` - string, specifies the url address for the app (default: \'/\')
  - :pyth:`assets_folder` - string, path to the folder containing static files, e.g. the .css file
  - :pyth:`run` - logical, whether the app should be run (if True) or return a dash.Dash (if False) instance of Flask app, useful for deploying the app on a server.

Usage
-----

Interface is quite simple. The plot area allows zooming, selection of points and exporting to png files
(in this cases removing axes may be useful).

Specific subsets (plotted as layers) of the data can be selected by clicking the legend
(double-click to isolate a specific subset)

Subsets are specified from the dropdown menu. To colour-code by subsets select 'categorical data' from the 'Colour points' menu.

Gene expression values can be chosen from the **Colour points** menu.

The app can also be embedded within an existing website using <iframe> or integrated into another Flask application.
In the latter case the arguments: :pyth:`server`, :pyth:`url_base_pathname`, :pyth:`assets_folder` are passed into the dash.Dash() call. Setting :pyth:`run` to FALSE will return the app object.
For details see:
https://dash.plot.ly/integrating-dash


Planned features
---------------
1. Displaying multiple adata files (selection from a dropdown menu)
2. 3d coordinate system
3. Comparison of gene expression levels across clusters - violin plots
