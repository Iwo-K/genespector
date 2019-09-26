import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import numpy as np

import sys

def make_layout(graphdata, main_title, subtitle):

    layout = html.Div([
        html.H4(children=main_title,
            style = {'textAlign': 'center'}),
        html.H6(children= subtitle,
            style = {'textAlign': 'center'}),
        html.Div([
            html.Div([
                dcc.Graph(id='maingraph',
                    config = {
                        'toImageButtonOptions' : {
                        'format' : 'png',
                        'filename' : 'newplot',
                        'height' : 1500,
                        'width' : 1500,
                        'scale' : 1.5
                        }
                    }
                ),
            ], className="two columns"),
            html.Div([
                html.Label(children='Split data by (categorical data):'),
                dcc.Dropdown(
                    id='split_dropdown',
                    options=[{'label' : x, 'value' : x} for x in graphdata['catinfo'].columns.values],
                    value=graphdata['catinfo'].columns.values[0]
                    ),
                html.Label(children='Colour points by (numeric data):'),
                dcc.Dropdown(
                    id='color_dropdown',
                    options=[{'label' : x, 'value' : x} for x in np.concatenate((np.array(['categorical data']), graphdata['numkeys'], graphdata['genekeys']))],
                    value='categorical data'
                    ),
                    html.Label(children='Choose the X and Y axis of the plot'),
                dcc.Dropdown(
                    id='Xaxis_choice',
                    options=[{'label' : x, 'value' : x} for x in graphdata['coords'].columns.values],
                    value = graphdata['coords'].columns.values[0]
                    ),
                dcc.Dropdown(
                    id='Yaxis_choice',
                    options=[{'label' : x, 'value' : x} for x in graphdata['coords'].columns.values],
                    value = graphdata['coords'].columns.values[1]
                    ),
                html.Label(children='Axis'),
                dcc.RadioItems(
                    id = 'Axis_switch',
                    options=[
                        {'label': 'ON', 'value': 'ON'},
                        {'label': 'OFF', 'value': 'OFF'}],
                    value='ON'),
                html.Label(children='Color scale: Min/Max'),
                dcc.Input(
                    id = 'cmin_input',
                    type = 'text',
                    value=''),
                dcc.Input(
                    id = 'cmax_input',
                    type = 'text',
                    value=''),
                html.Label(children='Choose point size'),
                dcc.Slider(
                    id='pointsize_slider',
                    min=2,
                    max=14,
                    value=7,
                    marks={str(x): str(x) for x in range(2, 14, 2)}),
                html.Br(),
                html.Label(children='Choose point opacity (%)'),
                dcc.Slider(
                    id='opacity_slider',
                    min=0,
                    max=100,
                    value=50,
                    marks={x : str(x) for x in range(0,110, 10)})
                    ], style = {'width' : 400, 'marginLeft': 200, 'marginTop' : 100, 'float' : 'left'}, className="two columns")
                    #], style = {'float': 'right'}, className="two columns")
                ], className="row"),
        html.Label(children='Genespector by: Iwo Kucinski, Gottgens lab',
                            style = {'textAlign': 'left', 'color' : '#CCCCCC'}),
        html.Label(children='Code: https://github.com/Iwo-K/genespector',
                            style = {'textAlign': 'left', 'color' : '#CCCCCC'}),
        html.Label(children='Based on Dash/Plotly framework (https://plot.ly/products/dash/)',
                            style = {'textAlign': 'left', 'color' : '#CCCCCC'})
])

    return(layout)
