import dash
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
from dash.dependencies import Input, Output

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

header = html.H1(
    "MarlinSim TherML"
)

left_accordion = dbc.Accordion([
    dbc.AccordionItem([
        html.Label("Mold"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number", required=True, size="sm")], width=4),
        ]),

        html.Label("Die"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number", required=True, size="sm")], width=4),
        ]),

        html.Label("Underfill"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number", required=True, size="sm")], width=4),
        ]),

        html.Label("Bumps"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number", required=True, size="sm")], width=4),
        ]),

        html.Label("Substrate"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number", required=True, size="sm")], width=4),
        ]),

        html.Label("Solder balls"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Y", type="number", required=True, size="sm")], width=4),
            dbc.Col([dbc.Input(placeholder="Z", type="number", required=True, size="sm")], width=4),
        ]),

        dbc.CardImg(src="https://drive.google.com/file/d/1-auPFqRIR9s0WCYqPaAplXvzNfQpqBCW/view?usp=sharing")

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
        dbc.Input(placeholder="Ambient temperature", type="number", required=True, value=25),

        dbc.Row([
        dbc.Col([
            html.Label("Start time"),
            dbc.Input(placeholder="Start time", type="number", value=0.0)]),

        dbc.Col([html.Label("End time"),
        dbc.Input(placeholder="End time", type="number")])]),
    ], title="Problem Setup"),
])

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

app.layout = dbc.Container(
    [
        navbar,
        html.Hr(),
        dbc.Col([
            dbc.Card([
                dbc.Col([
                    left_accordion
                ], width=12)], body=True)
        ], width=3),
        # dbc.Col([
        #     dbc.Card([
        #         db
        #     ])
        # ])
    ]
)


if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
