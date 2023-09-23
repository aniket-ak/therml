import dash
import dash_bootstrap_components as dbc
from dash import dcc
from dash import Input, Output, State, html
from dash.dependencies import Input, Output

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
# app.stylesheets.serve_locally=True
app.scripts.serve_locally = True

header = html.H1(
    "MarlinSim TherML"
)

left_accordion = dbc.Accordion([
    dbc.AccordionItem([
        html.Label("Mold"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number",
                    required=True, size="sm")], width=4),
        ]),

        html.Label("Die"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number",
                    required=True, size="sm")], width=4),
        ]),

        html.Label("Underfill"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number",
                    required=True, size="sm")], width=4),
        ]),

        html.Label("Bumps"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number",
                    required=True, size="sm")], width=4),
        ]),

        html.Label("Substrate"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number",
                    required=True, size="sm")], width=4),
        ]),

        html.Label("Solder balls"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number",
                    required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number",
                    required=True, size="sm")], width=4),
        ]),

        dbc.Card([
            html.Label("Package layout"),
            dbc.CardImg(src="http://localhost:8000/package_simplified.png")
        ], style={"margin-top": "7px"})


    ], title="Geometry and Model"),
    dbc.AccordionItem([
        html.Label("Top BC"),
        dcc.Dropdown(id="top-bc", options=[
            {"label": "Insulated", "value": "insulated"},
            {"label": "Const. T", "value": "const_T"},
            {"label": "Const. HTC", "value": "const_h"},
            {"label": "Const. heat flux", "value": "const_q"}]),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="BC Value", type="number"),]),
            dbc.Col([dbc.Input(placeholder="Reference temperature", type="number")])
        ]),

        html.Label("Bottom BC"),
        dcc.Dropdown(id="bottom-bc", options=[
            {"label": "Insulated", "value": "insulated"},
            {"label": "Const. T", "value": "const_T"},
            {"label": "Const. HTC", "value": "const_h"},
            {"label": "Const. heat flux", "value": "const_q"}]),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="BC Value", type="number"),]),
            dbc.Col([dbc.Input(placeholder="Reference temperature", type="number")])
        ]),
    ], title="Boundary Conditions"),
    dbc.AccordionItem([
        html.Label("Ambient temperature"),
        dbc.Input(placeholder="Ambient temperature",
                  type="number", required=True, value=25),

        dbc.Row([
            dbc.Col([
                html.Label("Start time"),
                dbc.Input(placeholder="Start time", type="number", value=0.0)]),

            dbc.Col([html.Label("End time"),
                     dbc.Input(placeholder="End time", type="number")])]),
    ], title="Problem Setup"),
])

settings_modal = dbc.Modal([
    dbc.ModalHeader(dbc.ModalTitle(
        "Model settings")),
    left_accordion,
    dbc.ModalFooter(
        dbc.Button(
            "Submit", id="close", className="ms-auto", n_clicks=0
        )
    ),
], id="modal", is_open=False,)

# bottom_console = dbc.

navbar = dbc.NavbarSimple(
    children=[
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("New Project"),
                dbc.DropdownMenuItem("Exit"),
            ],
            nav=True,
            in_navbar=True,
            label="File",
        ),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Units"),
            ],
            nav=True,
            in_navbar=True,
            label="Settings",
        ),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Documentation"),
                dbc.DropdownMenuItem("Raise SR"),
                dbc.DropdownMenuItem("Version"),
            ],
            nav=True,
            in_navbar=True,
            label="Help",
        ),
    ],
    brand="MarlinSim therML",
    brand_href="#",
    color="darkblue",
    dark=True,
)

upload_component = dcc.Upload(
    id='upload-data',
    children=html.Div([
        'Drag and Drop or ',
        html.A('Select Files')
    ]),
    style={
        'width': '100%',
        'height': '40px',
        'lineHeight': '40px',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'borderRadius': '5px',
        'textAlign': 'center',
        'margin': '1px'
    },
    multiple=True
)

app.layout = dbc.Container([
    navbar,
    html.Hr(),
    dbc.Tabs([
        dbc.Tab(
            dbc.Row([
                    dbc.Col([
                            dbc.Button("Settings", id="open", n_clicks=0),
                            settings_modal,
                            dbc.Col(["Upload Power dissipation files"]),
                            upload_component,
                            ], width=3),
                    dbc.Col([
                        dbc.Card([
                            dbc.Row([
                                dbc.Col([
                                ], width=6, id="scenario"),
                                dbc.Card([dbc.Col(["Visualization"])],)]),
                            dbc.Card(
                                [dbc.Col([dbc.Pagination(id="pagination", max_value=5),]),])
                        ])
                    ], width=9),

                    ]), label="Thermal Simulation"),
        dbc.Tab("ML flow", label="Machine Learning")
    ]),
])

# @app.callback(Output('output-data-upload', 'children'),
#               Input('upload-data', 'contents'),
#               State('upload-data', 'filename'),
#               State('upload-data', 'last_modified'))
# def update_output(list_of_contents, list_of_names, list_of_dates):
#     print(list_of_names)
#     return list_of_names

# @app.callback(
#     Output("modal", "is_open"),
#     [Input("open", "n_clicks"), Input("close", "n_clicks")],
#     [State("modal", "is_open")],
# )
# def toggle_modal(n1, n2, is_open):
#     if n1 or n2:
#         return not is_open
#     return is_open


@app.callback(
    Output("scenario", "children"),
    [Input("pagination", "active_page")],
)
def change_page(page):
    if page:
        return f"Page selected: {page}"
    return "Select a page"


if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
