# Importing modules
import dash
from dash import html
from dash import dcc
import pandas as pd
import plotly.graph_objs as go
import numpy as np
from .import_fun import import_adata
import sys
from .layout import make_layout
import scanpy as sc
from pathlib import Path
from anndata import AnnData

#Function for running the genespector

def make_app(adata = None,
    use_raw = False,
    make_var_unique = True,
    main_title = "Genespector",
    subtitle = None,
    layout = None,
    server = True,
    url_base_pathname = '/',
    run = True):

    #Specifying the current path to access the example dataset
    HERE = Path(__file__).parent

    #adata argument checking and processing
    if adata is None:
        adata_path = HERE / 'data/tiny_example2.h5ad'
        adata = sc.read(adata_path)

    elif isinstance(adata, str):
        adata = sc.read(adata)

    elif isinstance(adata, AnnData):
        pass

    else:
        print('Error: adata argument needs to be either an AnnData object or path ')
        return

    #Importing data
    graphdata = import_adata(adata, use_raw, make_var_unique)
    #Setting up dash server
    app = dash.Dash(__name__,
        server = server,
        url_base_pathname = url_base_pathname)

    #App Layout
    if subtitle is None:
        subtitle = 'Data: LK/LSK cells (Dahlin et al. 2018)'

    if layout is None:
        app.layout = make_layout(graphdata, main_title, subtitle)
    else:
        app.layout = layout(graphdata, main_title, subtitle)

    ''' App callback
    Specifies inputs and outputs from sliders and dropdown menus
    and generates graph objects, which are updated respectively
    '''
    @app.callback(
        dash.dependencies.Output('maingraph', 'figure'),
        [dash.dependencies.Input('split_dropdown', 'value'),
        dash.dependencies.Input('color_dropdown', 'value'),
        dash.dependencies.Input('Xaxis_choice', 'value'),
        dash.dependencies.Input('Yaxis_choice', 'value'),
        dash.dependencies.Input('pointsize_slider', 'value'),
        dash.dependencies.Input('Axis_switch', 'value'),
         dash.dependencies.Input('cmin_input', 'value'),
         dash.dependencies.Input('cmax_input', 'value'),
         dash.dependencies.Input('opacity_slider', 'value')])
    def update_figure(splitbyvalue, colorvalue, Xaxis, Yaxis, pointsize, axiswitchvalue, cmin_input, cmax_input, opacity):

        if colorvalue == 'categorical data':
            pass
        elif colorvalue in graphdata['numkeys']:
            colors = np.squeeze(np.array(graphdata['numinfo'][:,graphdata['numkeys'] == colorvalue]))
        else:
            colors = np.squeeze(np.array(graphdata['genexpr'][:,graphdata['genekeys'] == colorvalue].todense()))

        traces = []

        category = graphdata['catinfo'][splitbyvalue].cat.categories

        for i in category:

            boolSUB = graphdata['catinfo'][splitbyvalue] == i

            if colorvalue == 'categorical data':

                graph = go.Scattergl(
                    x = graphdata['coords'].loc[boolSUB,:][Xaxis].values,
                    y = graphdata['coords'].loc[boolSUB,:][Yaxis].values,
                    mode = 'markers',
                    marker = dict(opacity = opacity/100,
                        size = pointsize),
                    name = str(i))

                #Checking if there is matchig color code for the category
                if splitbyvalue in graphdata['color_catinfo'].keys():
                    marker_color = graphdata['color_catinfo'][splitbyvalue][i]
                    graph['marker']['color'] = marker_color

                traces.append(graph)

            else:
                colorSUB = colors[boolSUB]
                catinfoSUB = graphdata['catinfo'].loc[boolSUB,:]

                #Added possibility of having diverging colorscale
                if any(colors < 0) and any(colors > 0):
                # colorscale = 'RdBu'
                    colorscale = [[0, 'rgba(31, 58, 147, 1)'],
                    [-min(colors)/(max(colors)-min(colors)), '#DCDCDC'],
                    [1, 'rgba(207, 0, 15, 1)']]
                else:
                    colorscale = 'Viridis'#'Blackbody'

                # Specifying the markers
                cmin = min(colors)
                cmax = max(colors)
                if cmin_input != '':
                    cmin = float(cmin_input)
                if cmax_input != '':
                    cmax = float(cmax_input)
                marker = dict(color = colorSUB,
                    cmin = cmin,
                    cmax = cmax,
                    opacity = opacity/100,
                    colorscale = colorscale,
                    colorbar = dict(thickness=14, x = 1.15, title = colorvalue, titleside = 'bottom'), #perhaps make optional?, then just pass dict() to switch off
                    size = pointsize)

                traces.append(go.Scattergl(
                    x = graphdata['coords'].loc[boolSUB,:][Xaxis].values,
                    y = graphdata['coords'].loc[boolSUB,:][Yaxis].values,
                    mode = 'markers',
                    marker = marker,
                    name = str(i),
                    text = [catinfoSUB.index.values[x] + '<br>' +
                     splitbyvalue + ' ' + catinfoSUB[splitbyvalue].values[x] +
                       '<br>' + str(colorSUB[x]) for x in range(0, len(catinfoSUB.index))],
                        hoverinfo = 'text')
                            )

        if axiswitchvalue == 'ON':
            axiscolor = '#BBBBBB'
        else: axiscolor = '#ffffff'

        return {
            'data': traces,
            'layout': {
                'xaxis' : {'title': Xaxis, 'color' : axiscolor},
                'yaxis' : {'title': Yaxis, 'color' : axiscolor},
                'height' : 800,
                'width' : 800,
                'hovermode' : 'closest',
                'uirevision' : '{}{}'.format(Xaxis, Yaxis),
            }
        }
    if run == True:
        app.run_server(debug=False)
    else:
        return(app)
