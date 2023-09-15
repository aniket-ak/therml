import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(
    html.Title(["A big title"]),
    html.Div(className='row', children=[
        
        html.Div(className="four columns", children=[
                html.Label(['Dataset:'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown-dataset',
                    options=[
                        {'label': 'Diabetes', 'value': 'diabetes'},
                        {'label': 'Boston Housing', 'value': 'boston'},
                        {'label': 'Sine Curve', 'value': 'sin'}

                    ],
                    value='diabetes',
                    searchable=False,
                    clearable=False,
                ),
                html.Label('Slider', style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Slider(
                    min=0,
                    max=9,
                    marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
                    value=5,
                ),
            ]),

        
        html.Div(className="four columns", children=[
                html.Label(['Model Type'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown-select-model',
                    options=[
                        {'label': 'Linear Regression', 'value': 'linear'},
                        {'label': 'Lasso', 'value': 'lasso'},
                        {'label': 'Ridge', 'value': 'ridge'},
                        {'label': 'Polynomial', 'value': 'polynomial'},
                        {'label': 'elastic-net', 'value': 'elastic-net'},

                    ],
                    value='linear',
                    searchable=False,
                    clearable=False
                ),
                html.Label('Slider', style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Slider(
                    min=0,
                    max=9,
                    marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
                    value=5,
                ),
            ]),

        html.Div(className="four columns", children=[
                html.Label(['Add data'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown-custom-selection',
                    options=[
                        {'label': 'Add Training Data', 'value': 'training'},
                        {'label': 'Add Test Data', 'value': 'test'},
                        {'label': 'Remove Data point', 'value': 'remove'},
                    ],
                    value='training',
                    clearable=False,
                    searchable=False
                ),
                html.Label('Slider', style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Slider(
                    min=0,
                    max=9,
                    marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
                    value=5,
                ),
            ]),
    ],style=dict(display='flex')),
)


if __name__ == '__main__':
    app.run(debug=True)