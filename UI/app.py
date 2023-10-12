import dash
import dash_bootstrap_components as dbc
from dash import dcc
from dash import Input, Output, State, html, ctx, ALL
from dash.dependencies import Input, Output
import pandas as pd
import plotly.express as px
import numpy as np
import os
import subprocess
from dash.exceptions import PreventUpdate
import random


first_names = ["Blue", "Red", "Green", "Yellow", "Orange", "Purple", "Silver", "Golden", "Ruby", "Sapphire", "Emerald", "Diamond",
                "Amber", "Topaz", "Jade", "Pearl", "Opal", "Crimson", "Lavender", "Azure", "Cobalt", "Indigo", "Violet", "Magenta",
                  "Turquoise", "Coral", "Platinum", "Bronze", "Titanium", "Graphite", "Quartz", "Onyx", "Ivory", "Ebony", "Obsidian",
                    "Sardonyx", "Moonstone", "Peridot", "Aquamarine", "Citrine", "Garnet", "Agate", "Hematite", "Zircon", "Malachite", 
                    "Rhodonite", "Carnelian", "Aventurine", "Selenite", "Beryl", "Chrysoprase", "Serpentine", "Fluorite", "Kyanite", 
                    "Labradorite", "Amazonite", "Rhodochrosite", "Sunstone", "Turritella", "Unakite", "Flint", "Obsidian", "Lapis", "Apatite",
                      "Azurite", "Calcite", "Jasper", "Moonstone", "Larimar", "Celestite", "Sugilite", "Charoite", "Pietersite", "Thulite", 
                      "Howlite", "Dendritic", "Variscite", "Angelite", "Prehnite", "Septarian", "Pyrite", "Amethyst", "Kunzite", "Selenite", 
                      "Cavansite", "Uvarovite", "Nuummite", "Chrysocolla", "Aragonite"]

last_names = ["Lion", "Tiger", "Eagle", "Wolf", "Dragon", "Bear", "Falcon", "Fox", "Hawk", "Leopard", "Panther",
               "Lynx", "Cheetah", "Cobra", "Serpent", "Phoenix", "Viper", "Raven", "Heron", "Hound", "Stallion", "Turtle",
                 "Jaguar", "Puma", "Coyote", "Hedgehog", "Gazelle", "Kangaroo", "Zebra", "Ocelot", "Vulture", "Raccoon", "Penguin", 
                 "Badger", "Dolphin", "Mongoose", "Koala", "Platypus", "Gorilla", "Ostrich", "Polarbear", "Piranha", "Jellyfish", 
                 "Octopus", "Squirrel", "Seagull", "Shark", "Crab", "Pelican", "Seahorse", "Owl", "Peacock", "Sloth", "Walrus", "Komodo", 
                 "Anaconda", "Armadillo", "Bison", "Buffalo", "Camel", "Cockatoo", "Dalmatian", "Dingo", "Elephant", "Ferret", "Giraffe", 
                 "Hamster", "Iguana", "Kookaburra", "Llama", "Meerkat", "Narwhal", "Panda", "Rattlesnake", "Salamander", "Sloth", "Tapir", 
                 "Wombat", "Yak", "Zorilla", "Bobcat", "Chameleon", "Firefly", "Jaguarundi", "Lemur", "Mantis", "Ocelot", "Platypus", "Scorpion",
                   "Squid", "Tarsier"]

random_first_name = random.choice(first_names)
random_last_name = random.choice(last_names)
random_name = random_first_name + "_" + random_last_name

df_tables = pd.DataFrame()
df_tables["Scenario"] = ["None"]
df_tables["Status"] = ["None"]
df_tables["Progress"] = [dbc.Progress(value=0, striped=True, id="progress0")]
df_tables["Simulate"] = [dbc.Button("Simulate", color="primary", className="me-1", id="simulate0",disabled=True)]

app = dash.Dash(
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.1"},
    ],
)
# app.stylesheets.serve_locally=True
app.scripts.serve_locally = True
app._favicon = "favicon.ico"
app.title = "therML"

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
        dbc.Select(id="top-bc", options=[
            {"label": "Insulated", "value": "insulated"},
            {"label": "Const. T", "value": "const_T"},
            {"label": "Const. HTC", "value": "const_h"},
            {"label": "Const. heat flux", "value": "const_q"}]),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="BC Value", type="number"),]),
            dbc.Col([dbc.Input(placeholder="Reference temperature", type="number")])
        ]),

        html.Label("Bottom BC"),
        dbc.Select(id="bottom-bc", options=[
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
    color="#0c2340",
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
                            dbc.Row([dbc.Col([dbc.Label("Model Settings")]), 
                                     dbc.Col([dbc.Button("Settings", id="open", n_clicks=0),])
                                     ],style={"margin-top": "30px"}),
                            settings_modal,
                            dbc.Label("Select the working directory"),
                            dbc.Input(placeholder="Working Directory", type="text", id="working_dir", required=True),
                            dbc.Col(["Upload Power dissipation files"]),
                            upload_component,
                            dbc.Label("Enter a name for simulation"),
                            dbc.Input(placeholder=random_name, type="text", id="simulation_name")
                            ], width=3),
                    dbc.Col([
                            dbc.Row([dbc.Col(html.H4("Scenarios") , width=6), 
                                     dbc.Col([
                                    dbc.Select(options=[{"label":i,"value":i} for i in ["None"]], id="viz_dropdown")
                                    ], width=6, id="dropdowns"),],style={"margin-top": "30px"}),
                            
                            dbc.Row([
                                dbc.Col([dbc.Table.from_dataframe(df_tables, striped=True, bordered=True, hover=True)], width=6, id="scenario"),
                                dbc.Col([dcc.Graph(id="contour")], width=6, id="visualization"),
                            ],style={"margin-top": "30px"}),
                            html.Hr(),
                            dbc.Col([dbc.Pagination(id="pagination", max_value=5),]),
                        
                    ], width=9),

                    ]), label="Thermal Simulation"),
        dbc.Tab("ML flow", label="Machine Learning")
    ]),
])


@app.callback(Output('scenario', 'children'),
              Output('dropdowns', 'children'),
              Output('pagination', 'max_value'),
              Input('upload-data', 'filename'),
              State('upload-data', 'filename'),
              Input('pagination', 'active_page'))
def update_output(list_of_names, list_of_names1, active_page):
    global df_tables
    if list_of_names is not None:
        df_tables = pd.DataFrame()
        df_tables["Scenario"] = [dbc.Col([name], id={"type":"scenario", "index":i}) for (i,name) in enumerate(list_of_names)]
        df_tables["Status"] = [dbc.Badge("Complete", color="white",text_color="success") for i in list_of_names]
        # TODO - Save the dataframe state to a file and read it here for the progress to be overwritten every time pagination click happens
        df_tables["Progress"] = [dbc.Progress(value=10, striped=True, id={"type":"progress","index":i}) for (i,name) in enumerate(list_of_names)]
        df_tables["Simulate"] = [dbc.Button("Simulate", color="primary", className="me-1", id={"type":"simulate","index":i}) for (i,name) in enumerate(list_of_names)]
        children = [dbc.Select(options=[{"label":i,"value":i} for i in list_of_names],id="viz_dropdown")]
    else:
        df_tables = pd.DataFrame()
        df_tables["Scenario"] = ["None"]
        df_tables["Status"] = ["None"]
        df_tables["Progress"] = [dbc.Progress(value=0, striped=True, id="progress0")]
        df_tables["Simulate"] = [dbc.Button("Simulate", color="primary", className="me-1", id="simulate0",disabled=True)]
        children = [dbc.Select(options=[{"label":i,"value":i} for i in ["None"]],id="viz_dropdown")]
    
    if active_page is None and list_of_names is not None:
        active_page=1
        table_ = dbc.Table.from_dataframe(df_tables[(active_page-1)*10:(active_page)*10], striped=True, bordered=True, hover=True)
    elif active_page is None and list_of_names is None:
        table_ = dbc.Table.from_dataframe(df_tables[0:10], striped=True, bordered=True, hover=True)
    else:
        # TODO - Cleanup logic here
        table_ = dbc.Table.from_dataframe(df_tables[(active_page-1)*10:(active_page)*10], striped=True, bordered=True, hover=True)

    max_value = int(df_tables.shape[0]/10)+1
    
    return table_, children, max_value

@app.callback(Output('contour', 'figure'),
              Input('viz_dropdown', 'value'))
def update_output(name):
    import time
    time.sleep(1)

    if name is not None:
        path = "/Users/aniket/Documents/MarlinSim/04_testing/scenarios"
        csv_ = os.path.join(path, name)
        matrix = pd.read_csv(csv_, header=None).values
        fig = px.imshow(matrix,labels=dict(x="X", y="Y",color="Power"),color_continuous_scale='jet')
    else:
        fig = px.area()
    return fig

@app.callback(
    Output({"type": "progress", "index": ALL}, "value"), 
    Input({"type": "simulate", "index": ALL}, "n_clicks"),
    Input({"type": "progress", "index": ALL}, "value"),
    Input({"type": "scenario", "index": ALL}, "children"),
    State('pagination', 'active_page'),
    State('working_dir', 'value'),
    State('simulation_name', 'value'))
def filter_heatmap(n_clicks, values, children, active_page, working_dir, simulation_name):
    if simulation_name is None:
        simulation_name = random_name
    
    output_dir = os.path.join(working_dir,simulation_name)
    log_dir = os.path.join(output_dir, "Logs")
    sol_dir = os.path.join(output_dir, "Solution")
    temp_dir = os.path.join(output_dir, "Temp")
    os.mkdir(output_dir)
    os.mkdir(log_dir)
    os.mkdir(sol_dir)
    os.mkdir(temp_dir)

    if len(n_clicks) < 1:
        raise PreventUpdate
    n_clicks = ctx.triggered[0]["value"]
    if not n_clicks:
        raise PreventUpdate
    button_id = ctx.triggered_id.index
    if active_page is None:
        active_page=1
    
    mod_values = values
    if isinstance(button_id, int):
        mod_values[button_id-(active_page-1)*10] = 100
    
    power_file = os.path.join(working_dir, children[button_id][0])
    print("Executing power file - ", power_file)
    
    result = subprocess.run(["julia", "/Users/aniket/Documents/MarlinSim/03_code/therml/3d/3d_fvm.jl", 
                             "-t", "4", "-working_dir", working_dir, "-power", power_file, "-run_name", 
                             simulation_name], capture_output=True, text=True)
    # julia --project=./therml_environment ./therml_environment/precompile_.jl -t 4 -working_dir 
    #   /Users/aniket/Documents/MarlinSim/04_testing/scenarios -power file_1.csv -run_name "sim_1"

    # Check for errors
    if result.returncode != 0:
        print("Error:", result.stderr)
        return None

    return mod_values


@app.callback(
    Output("working_dir", "invalid"),
    Input("working_dir", "value")
)
def check_working_dir(dir):
    if dir == "":
        return True
    else:
        return False

@app.callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


# @app.callback(
#     Output("scenario", "children"),
#     [Input("pagination", "active_page")],
# )
# def change_page(page):
#     if page:
#         return f"Page selected: {page}"
#     return "Select a page"


if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
