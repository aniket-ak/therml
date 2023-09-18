import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Col([dbc.Card([dbc.Col([
            dbc.Accordion([
                dbc.AccordionItem([
                    html.Label("Top BC"),
                    dcc.Dropdown(id="top-bc", options=[{"label": "Option 1", "value": "option1"}, {"label": "Option 2", "value": "option2"}]),
                    
                    html.Label("Bottom BC"),
                    dcc.Dropdown(id="bottom-bc", options=[{"label": "Option 1", "value": "option1"}, {"label": "Option 2", "value": "option2"}]),
                ], title= "Boundary Conditions")
            ]),


        ], width=12)], body=True)], width=3)


if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
