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
from functions import *
import time
import json
import plotly.graph_objects as go

working_dir_proj = ""
run_name_proj = ""

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
app.title = "Marlinsim therML"

header = html.H1(
    "MarlinSim TherML"
)

material_file = os.path.join("../", "3d/common/materials.json")
with open(material_file, 'r') as file:
    data = json.load(file)

material_names = [list(i.keys())[0] for i in data['materials']]
material_names.sort()

toast = dbc.Toast(
        "Settings saved successfully",id="save_success_toast",header="Notification",
        is_open=False,dismissable=True, icon="success", duration=4000, style={"position": "fixed", 
        "top": 66, "right": 10, "width": 350},
    )

left_accordion = dbc.Accordion([
    dbc.AccordionItem([
        dbc.Row([
            dbc.Col([html.Label("Mold [mm]")],width=2),
            dbc.Col([""],width=2),
            dbc.Col([""],width=2),
            dbc.Col([""],width=2),
            dbc.Col(["Material"],width=2),
            dbc.Col(["Surface Material"],width=2)
        ]),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", size="sm",id="mold_x",value=15)], width=2),
            dbc.Col([dbc.Input(placeholder="Y", type="number", size="sm",id="mold_y",value=15)], width=2),
            dbc.Col([dbc.Input(placeholder="Z", type="number", size="sm",id="mold_z",value=5)], width=2),
            dbc.Col([dbc.Select(material_names, "Epoxy Molding Compound (EMC)", id="mold_material")], width=3),
            dbc.Col([dbc.Select(material_names, "Epoxy Molding Compound (EMC)", id="mold_surface_material")], width=3),
        ]),

        html.Label("Die [mm]"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", size="sm",id="die_x",value=10)], width=2),
            dbc.Col([dbc.Input(placeholder="Y", type="number", size="sm",id="die_y",value=10)], width=2),
            dbc.Col([dbc.Input(placeholder="Z", type="number", size="sm",id="die_z",value=4)], width=2),
            dbc.Col([dbc.Select(material_names, "Silicon", id="die_material")], width=3),
            dbc.Col([dbc.Select(material_names, "FR4", id="die_surface_material")], width=3)
        ]),

        html.Label("Epoxy [mm]"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", size="sm",id="underfill_x",value=12)], width=2),
            dbc.Col([dbc.Input(placeholder="Y", type="number", size="sm",id="underfill_y",value=12)], width=2),
            dbc.Col([dbc.Input(placeholder="Z", type="number", size="sm",id="underfill_z",value=4)], width=2),
            dbc.Col([dbc.Select(material_names, "Epoxy Molding Compound (EMC)", id="underfill_material")], width=3),
            dbc.Col([dbc.Select(material_names, "FR4", id="epoxy_surface_material")], width=3)
        ]),

        html.Label("Bumps [mm]"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", size="sm",id="bumps_x",value=10)], width=2),
            dbc.Col([dbc.Input(placeholder="Y", type="number", size="sm",id="bumps_y",value=10)], width=2),
            dbc.Col([dbc.Input(placeholder="Z", type="number", size="sm",id="bumps_z",value=3)], width=2),
            dbc.Col([dbc.Select(material_names, "SnPb Alloy", id="bumps_materials")], width=3),
            dbc.Col([dbc.Select(material_names, "FR4", id="bumps_surface_materials")], width=3)
        ]),

        html.Label("Substrate [mm]"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", size="sm",id="substrate_x",value=15)], width=2),
            dbc.Col([dbc.Input(placeholder="Y", type="number", size="sm",id="substrate_y",value=15)], width=2),
            dbc.Col([dbc.Input(placeholder="Z", type="number", size="sm",id="substrate_z",value=5)], width=2),
            dbc.Col([dbc.Select(material_names, "FR4", id="substrate_materials")], width=3),
            dbc.Col([dbc.Select(material_names, "FR4", id="substrate_surface_materials")], width=3)
        ]),

        html.Label("Solder balls [mm]"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="X", type="number", size="sm",id="solder_x",value=15)], width=2),
            dbc.Col([dbc.Input(placeholder="Y", type="number", size="sm",id="solder_y",value=15)], width=2),
            dbc.Col([dbc.Input(placeholder="Z", type="number", size="sm",id="solder_z",value=2)], width=2),
            dbc.Col([dbc.Select(material_names, "SnPb Alloy", id="solder_materials")], width=3),
            dbc.Col([dbc.Select(material_names, "FR4", id="solder_surface_materials")], width=3)
        ]),

        dbc.Card([
            html.Label("Package layout"),
            dcc.Graph(id="geometry")
        ], style={"margin-top": "7px"})


    ], title="Geometry and Model"),
    dbc.AccordionItem([
        html.Label("Top BC"),
        dbc.Select([
            {"label": "Insulated", "value": "insulated"},
            {"label": "Const. T [C]", "value": "const_T"},
            {"label": "Const. HTC [W/m2K]", "value": "const_h"},
            {"label": "Const. heat flux [W/m2]", "value": "const_q"}], ["const_h"], id="top_bc_type"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="BC Value", id="top_bc_value", value=5, type="number"),]),
            dbc.Col([dbc.Input(placeholder="Reference temperature", id="top_bc_ref_temp", type="number", value=25)])
        ]),

        html.Label("Bottom BC"),
        dbc.Select([
            {"label": "Insulated", "value": "insulated"},
            {"label": "Const. T [C]", "value": "const_T"},
            {"label": "Const. HTC [W/m2K]", "value": "const_h"},
            {"label": "Const. heat flux [W/m2]", "value": "const_q"}], ["const_h"],id="bottom_bc_type"),
        dbc.Row([
            dbc.Col([dbc.Input(placeholder="BC Value", id="bot_bc_value", value=5, type="number"),]),
            dbc.Col([dbc.Input(placeholder="Reference temperature", id="bot_bc_ref_temp", type="number", value=25)])
        ]),
    ], title="Boundary Conditions"),
    dbc.AccordionItem([
        dbc.Row([
            dbc.Col([
                html.Label("Ambient temperature"),
                dbc.Input(placeholder="Ambient temperature [C]",id="ambient_temp", type="number", required=True, value=25),
            ]),
            dbc.Col([
                html.Label("Start time [s]"),
                dbc.Input(placeholder="Start time",id="start_time", type="number", value=0.0)]),

            dbc.Col([
                html.Label("End time [s]"),
                dbc.Input(placeholder="End time",id="end_time", type="number", value=10.0)])]),
    ], title="Problem Setup"),
])

settings_modal = dbc.Modal([
    dbc.ModalHeader(dbc.ModalTitle(
        "Model and Problem Settings")),
    left_accordion,
    toast,
    dbc.ModalFooter(
        dbc.Button(
            "Submit", id="save_settings", className="ms-auto", n_clicks=0
        )
    ),
], id="modal", is_open=False, size="xl")

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

viz_modal = dbc.Modal([
    dbc.ModalHeader(dbc.ModalTitle(
        "Power and Temperature Visualization"), style={"margin-top": "10px"}),
    dbc.Col([
        dbc.Row([dbc.Col("Select the scenario from dropdown"),dbc.Col([dbc.Select(options=[{"label":i,"value":i} for i in ["None"]], id="viz_dropdown")])])
        ], width=6, id="dropdowns", align="center", style={"margin-top": "10px"}),
    dbc.Col([dcc.Graph(id="contour")], width=12, id="visualization"),
    dbc.Row([
        dbc.Col(["Select time step"],width=2), 
        dbc.Col(dbc.Input(placeholder=0, type="number", size="sm", id="viz_time_value"), width=1), 
        dbc.Col("[s]", width=1),
        dbc.Col("Cut-plane", width=1),
        dbc.Col(dbc.Select(["XY", "YZ", "ZX"],"XY", id="viz_cutplane"), width=1),
        dbc.Col(["Select cut-section"], width =2), 
        dbc.Col([dcc.Slider(0, 1, id='cut_section_slider',tooltip={"placement": "bottom", "always_visible": True})], width=4)]),
    dbc.Col([dcc.Graph(id="temp_contour")], width=12, id="visualization_temperature"),
    dbc.ModalFooter(
        dbc.Button(
            "Close", id="viz_close", className="ms-auto", n_clicks=0
        )
    ),
], id="viz_modal", is_open=False,size="xl")

app.layout = dbc.Container([
    navbar,
    html.Hr(),
    dcc.Interval(id='interval-component', interval=5*1000, n_intervals=0),
    dbc.Tabs([
        dbc.Tab(
            dbc.Row([
                    dbc.Col([
                            dbc.Label("1. Select the working directory"),
                            dbc.Input(placeholder="Working Directory", type="text", id="working_dir", required=True),
                            dbc.Row([dbc.Col([dbc.Label("2. Model Settings")]), 
                                     dbc.Col([dbc.Button("Settings", id="open", n_clicks=0),])
                                     ],style={"margin-top": "30px"}, justify="between"),
                            settings_modal,
                            dbc.Col(["3. Upload Power dissipation files"]),
                            upload_component,
                            dbc.Label("4. Enter a name for simulation"),
                            dbc.Input(placeholder=random_name, type="text", id="run_name")
                            ], width=3),
                    dbc.Col([
                            dbc.Row([dbc.Col(html.H4("5. Scenarios") , width=6), dbc.Col([dbc.Button("6. Visualization", id="viz_open", n_clicks=0),], width="auto"), viz_modal
                                     ],style={"margin-top": "30px"}, justify="between"),
                            dbc.Row([
                                dbc.Col([dbc.Table.from_dataframe(df_tables, striped=True, bordered=True, hover=True)], width=12, id="scenario"),
                                dbc.Toast(
                                        "Simulation submitted",id="success_toast",header="Notification",
                                        is_open=False,dismissable=True, icon="success", duration=4000, style={"position": "fixed", 
                                        "top": 66, "right": 10, "width": 350},
                                    ),
                            ],style={"margin-top": "30px"}),
                            html.Hr(),
                            dbc.Col([dbc.Pagination(id="pagination", max_value=5),]),
                        
                    ], width=9),

                    ]), label="Thermal Simulation"),
        dbc.Tab("ML flow", label="Machine Learning")
    ]),
], fluid=True)

@app.callback(Output('scenario', 'children'),
              Output('dropdowns', 'children'),
              Output('pagination', 'max_value'),
              Input('upload-data', 'filename'),
              State('upload-data', 'filename'),
              Input('pagination', 'active_page'))
def update_output(list_of_names, list_of_names1, active_page):
    global df_tables
    children = [dbc.Col(["Select scenario from dropdown"]),dbc.Select(options=[{"label":i,"value":i} for i in ["None"]],id="viz_dropdown")]
    
    triggered_id = ctx.triggered_id
    cols_to_display = ["Scenario","Status","Progress","Simulate"]

    if triggered_id == "upload-data":
        df_tables = pd.DataFrame()
        df_tables["Scenario"] = list_of_names #[dbc.Col([name], id={"type":"scenario", "index":i}) for (i,name) in enumerate(list_of_names)]
        df_tables["status_internal"] = ["yet_to_start" for i in range(len(list_of_names))]
        df_tables["Status"] = [dbc.Badge("Not started", color="secondary",text_color="white", id={"type":"status", "index":i}) for (i,name) in enumerate(list_of_names)]
        df_tables["Progress"] = [dbc.Progress(value=0, striped=True, id={"type":"progress","index":i}) for (i,name) in enumerate(list_of_names)]
        df_tables["Simulate"] = [dbc.Button("Simulate", color="primary", className="me-1", id={"type":"simulate","index":i}) for (i,name) in enumerate(list_of_names)]
        table_ = dbc.Table.from_dataframe(df_tables[cols_to_display][0:10], striped=True, bordered=True, hover=True)
    else:
        if active_page is None and list_of_names is not None:
            active_page=1
            table_ = dbc.Table.from_dataframe(df_tables[cols_to_display][(active_page-1)*10:(active_page)*10], striped=True, bordered=True, hover=True)
        elif active_page is None and list_of_names is None:
            table_ = dbc.Table.from_dataframe(df_tables[cols_to_display][0:10], striped=True, bordered=True, hover=True)
        else:
            # TODO - Cleanup logic here
            table_ = dbc.Table.from_dataframe(df_tables[cols_to_display][(active_page-1)*10:(active_page)*10], striped=True, bordered=True, hover=True)
    children = ["Select scenario from dropdown",dbc.Select(options=[{"label":i,"value":i} for i in df_tables["Scenario"].to_list()],id="viz_dropdown")]
    max_value = int(df_tables.shape[0]/10)+1
    
    return table_, children, max_value

@app.callback(Output('contour', 'figure'),
              Output('temp_contour', 'figure'),
              Input('viz_dropdown', 'value'),
              Input('cut_section_slider', 'value'),
              Input('viz_time_value', 'value'),
              Input('viz_cutplane', 'value'))
def update_output(name, slider_value, viz_time_value, viz_cut_plane):
    time.sleep(1)
    global working_dir_proj

    if name is not None and name != "None":
        csv_ = os.path.join(working_dir_proj, name)
        matrix = pd.read_csv(csv_, header=None).values
        fig = px.imshow(matrix,labels=dict(x="X", y="Y",color="Power"),color_continuous_scale='jet')
    else:
        fig = px.area()
    
    if name is not None and name != "None":
        sol_name = name+"__solution.sol"
        solution_dir = os.path.join(os.path.join(working_dir_proj, run_name_proj), "Solution")
        file = os.path.join(solution_dir, sol_name)
        print("File name", file)
        if slider_value is None:
            slider_value = 0
        if viz_time_value is None:
            viz_time_value = 1
        if os.path.exists(file):
            solution_, time_ = read_solution(file)
            contour_time = time_[np.where(time_==viz_time_value)]
            contour_index = int(np.where(time_==viz_time_value)[0])
            if contour_time is not []:
                print("inside plotting 1")
                if viz_cut_plane == "XY":
                    print("XY")
                    contour = solution_[contour_index, 1:-1, 1:-1, int(slider_value*solution_.shape[3])]
                    print("contour", contour.shape)
                    print(contour)
                    fig_temp = px.imshow(contour,labels=dict(x="X", y="Y",color="Temperature"),color_continuous_scale='jet')
                elif viz_cut_plane == "XZ":
                    print("XZ")
                    contour = solution_[contour_index, int(slider_value*solution_.shape[1]), 1:-1, 1:-1]
                    fig_temp = px.imshow(contour,labels=dict(x="Z", y="X",color="Temperature"),color_continuous_scale='jet')
                else:
                    print("YZ")
                    contour = solution_[contour_index, 1:-1, int(slider_value*solution_.shape[2]), 1:-1]
                    fig_temp = px.imshow(contour,labels=dict(x="Z", y="Y",color="Temperature"),color_continuous_scale='jet')
                
            else:
                print("else1")
                fig_temp = px.area()
        else:
            print("else2")
            fig_temp = px.area()
    else:
        print("else 3")
        fig_temp = px.area()
    return fig, fig_temp

@app.callback(
    Output("success_toast", "is_open"),
    Input({"type": "simulate", "index": ALL}, "n_clicks"),
    # State({"type": "scenario", "index": ALL}, "children"),
    State('pagination', 'active_page'),
    State('working_dir', 'value'),
    State('run_name', 'value'),
    State({"type":"status", "index":ALL}, "children"),
    State({"type": "status", "index": ALL}, "color"),)
def filter_heatmap(n_clicks, active_page, working_dir, run_name, status,colors):
    if run_name is None:
        run_name = random_name
    
    global temp_dir
    global working_dir_proj
    global run_name_proj
    working_dir_proj = working_dir
    run_name_proj = run_name

    if working_dir is not None:
        output_dir = os.path.join(working_dir,run_name)
        log_dir = os.path.join(output_dir, "Logs")
        sol_dir = os.path.join(output_dir, "Solution")
        temp_dir = os.path.join(output_dir, "Temp")
    else:
        return False

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    if not os.path.exists(sol_dir):
        os.mkdir(sol_dir)
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    if len(n_clicks) < 1:
        raise PreventUpdate
    n_clicks = ctx.triggered[0]["value"]
    if not n_clicks:
        raise PreventUpdate
    button_id = ctx.triggered_id.index
    if active_page is None:
        active_page=1
    # print("active_page:", active_page, "button_id", button_id, "children len", len(children))
    if button_id>10:
        scenario_name = df_tables["Scenario"][button_id] #children[button_id-active_page*10][0]
    else:
        scenario_name = df_tables["Scenario"][button_id] #children[button_id][0]
    df_tables.loc[button_id,'status_internal'] = 'wip'
    print("Executing power file - ", scenario_name)
    
    # result = subprocess.run(["julia", "--project=/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment", 
    #                          "/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment/precompile_.jl", 
    #                          "-t", "4", "-working_dir", working_dir, "-power", scenario_name, "-run_name", 
    #                          run_name], capture_output=True, text=True)

    # run_julia(working_dir=working_dir, scenario_name=scenario_name, run_name=run_name)
    output_queue = queue.Queue()
    julia_thread = threading.Thread(target=run_julia, args=(working_dir, scenario_name, run_name, output_queue))
    julia_thread.start()
    print("Captured output:")
    captured_output = output_queue.get()
    print(captured_output)

    julia_thread.join()

    
    # julia --project=./therml_environment /Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment/precompile_.jl -t 4 -working_dir 
    #   /Users/aniket/Documents/MarlinSim/04_testing/scenarios -power file_1.csv -run_name "sim_1"

    return True

@app.callback(
    Output({"type": "progress", "index": ALL}, "value"),
    Output({"type": "progress", "index": ALL}, "label"),
    Output({"type": "status", "index": ALL}, "children"), 
    Output({"type": "status", "index": ALL}, "color"),
    State('pagination', 'active_page'),
    State({"type": "progress", "index": ALL}, "value"),
    #State({"type": "scenario", "index": ALL}, "children"),
    [Input('interval-component', 'n_intervals')],
    State({"type":"status", "index":ALL}, "children"),
    State({"type": "status", "index": ALL}, "color"),
    Input({"type": "simulate", "index": ALL}, "n_clicks"),)
def timer(active_page, values, n, status,colors, n_clicks):
    if active_page is None:
        active_page=1
    children = df_tables.loc[(active_page-1)*10:(active_page)*10,"Scenario"].tolist() #[i[0] for i in children]
    current_progress = {s:v for (s,v) in zip(children, values)} 
    
    if 'temp_dir' in globals():
        if os.path.exists(temp_dir) and len(os.listdir(temp_dir))>0:
            list_of_files = os.listdir(temp_dir)
            for i in list_of_files:
                if i.split("__")[0] in children:
                    f=open(os.path.join(temp_dir,i), "r")
                    content = f.readlines()

                    if len(content)>0:
                        latest_progress = float(content[-1].split("\n")[0])
                        current_progress[i.split("__")[0]] = latest_progress
                    f.close()
    # print("progress ",current_progress)
    n_clicks = ctx.triggered[0]["value"]
    if not n_clicks:
        raise PreventUpdate
    button_id = ctx.triggered_id.index
    mod_status = status
    mod_colors = colors

    for i, (s, v) in enumerate(current_progress.items()):
        row_num = int(df_tables[df_tables["Scenario"] == s].index[0])
        if v > 99.9:
            mod_status[row_num-(active_page-1)*10] = "Complete"
            mod_colors[row_num-(active_page-1)*10] = "success"
            df_tables.loc[row_num,"Status"] = dbc.Badge("Complete", color="success",text_color="white", id={"type":"status", "index":row_num})
        elif df_tables.loc[row_num,"status_internal"] == 'wip':
            mod_status[row_num-(active_page-1)*10] = "In Progress"
            mod_colors[row_num-(active_page-1)*10] = "primary"
            df_tables.loc[row_num,"Status"] = dbc.Badge("In Progress", color="primary",text_color="white", id={"type":"status", "index":row_num})

    # df_tables["Scenario"] = [dbc.Col([name], id={"type":"scenario", "index":i}) for (i,name) in enumerate(list_of_names)]
    
    # df_tables["Progress"] = [dbc.Progress(value=0, striped=True, id={"type":"progress","index":i}) for (i,name) in enumerate(list_of_names)]
    # df_tables["Simulate"] = [dbc.Button("Simulate", color="primary", className="me-1", id={"type":"simulate","index":i}) for (i,name) in enumerate(list_of_names)]
    # df_tables.to_csv('./df_tables1.csv')
    return list(current_progress.values()), list(current_progress.values()), mod_status, mod_colors

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
    Output("save_success_toast", "is_open"),
    [Input("save_settings", "n_clicks")],
    [State('working_dir', 'value'),State("modal", "is_open"), State("open", "n_clicks"),
    State("mold_x","value"), State("mold_y","value"), State("mold_z","value"), State("die_x","value"), State("die_y","value"), State("die_z","value"),
    State("underfill_x","value"), State("underfill_y","value"), State("underfill_z","value"), State("bumps_x","value"), State("bumps_y","value"), State("bumps_z","value"), 
    State("substrate_x","value"), State("substrate_y","value"), State("substrate_z","value"), State("solder_x","value"), State("solder_y","value"), State("solder_z","value"),
    State("mold_material", "value"), State("die_material", "value"), State("underfill_material", "value"), State("bumps_materials", "value"), State("substrate_materials", "value"),
    State("solder_materials", "value"), State("mold_surface_material", "value"), State("die_surface_material", "value"), State("epoxy_surface_material", "value"),
    State("bumps_surface_materials", "value"), State("substrate_surface_materials", "value"), State("solder_surface_materials", "value"), State("ambient_temp", "value"),
    State("start_time", "value"), State("end_time", "value"), State("top_bc_type", "value"), State("top_bc_value", "value"), State("top_bc_ref_temp", "value"),
    State("bottom_bc_type", "value"), State("bot_bc_value", "value"), State("bot_bc_ref_temp", "value"),
])
def toggle_modal(n1, working_dir, is_open, n2, mold_x, mold_y, mold_z, die_x, die_y, die_z, underfill_x, underfill_y, underfill_z,
                 bumps_x, bumps_y, bumps_z, substrate_x, substrate_y, substrate_z, solder_x, solder_y, solder_z, 
                 mold_material, die_material, underfill_material, bumps_materials, substrate_materials, solder_materials,
                 mold_surface_material, die_surface_material, epoxy_surface_material, bumps_surface_materials, substrate_surface_materials, solder_surface_materials,
                 ambient_temp, start_time, end_time, top_bc_type, top_bc_value, top_bc_ref_temp, bottom_bc_type, bot_bc_value, bot_bc_ref_temp):
    global working_dir_proj
    working_dir_proj = working_dir
    if n1:
        if working_dir_proj != "" and working_dir_proj is not None:
            settings_file = os.path.join(working_dir_proj, "settings.json")
            save_to_json(settings_file, (mold_x, mold_y, mold_z, die_x, die_y, die_z, underfill_x, underfill_y, underfill_z,
                    bumps_x, bumps_y, bumps_z, substrate_x, substrate_y, substrate_z, solder_x, solder_y, solder_z, 
                    mold_material, die_material, underfill_material, bumps_materials, substrate_materials, solder_materials,
                    mold_surface_material, die_surface_material, epoxy_surface_material, bumps_surface_materials, substrate_surface_materials, solder_surface_materials,
                    ambient_temp, start_time, end_time, top_bc_type, top_bc_value, top_bc_ref_temp, bottom_bc_type, bot_bc_value, bot_bc_ref_temp))
            return True
        else:
            return False
    else:
        return False

@app.callback(
    Output("geometry", "figure"),
    [Input("mold_x","value"), Input("mold_y","value"), Input("mold_z","value"), Input("die_x","value"), Input("die_y","value"), Input("die_z","value"),
    Input("underfill_x","value"), Input("underfill_y","value"), Input("underfill_z","value"), Input("bumps_x","value"), Input("bumps_y","value"), Input("bumps_z","value"), 
    Input("substrate_x","value"), Input("substrate_y","value"), Input("substrate_z","value"), Input("solder_x","value"), Input("solder_y","value"), Input("solder_z","value"),
])
def toggle_modal(mold_x, mold_y, mold_z, die_x, die_y, die_z, underfill_x, underfill_y, underfill_z,
                 bumps_x, bumps_y, bumps_z, substrate_x, substrate_y, substrate_z, solder_x, solder_y, solder_z):
    fig = go.Figure()
    
    #8dd3c7 #f1f0bc #d9a0ac #a9a1b3 #efb46f #bbdc77 #f0d1e1 #c9a8c9 #c8d3c3 #ffed6f
    if None not in [mold_x, mold_y, mold_z, die_x, die_y, die_z, underfill_x, underfill_y, underfill_z, bumps_x, bumps_y, bumps_z, substrate_x, substrate_y, substrate_z, solder_x, solder_y, solder_z]:
        fig.add_shape(type="rect", x0=-solder_x/2, y0=0, x1=solder_x/2, y1=solder_z, line=dict( color="black", width=2), fillcolor="#8dd3c7",opacity=0.5,label=dict(text="solder", textposition="top left"))
        fig.add_shape(type="rect", x0=-substrate_x/2, y0=solder_z, x1=substrate_x/2, y1=solder_z+substrate_z, line=dict( color="black", width=2), fillcolor="#f1f0bc",opacity=0.5, label = dict(text="substrate",textposition="top left"))
        fig.add_shape(type="rect", x0=-mold_x/2, y0=solder_z+substrate_z, x1=mold_x/2, y1=solder_z+substrate_z+mold_z, line=dict( color="black", width=2), fillcolor="#d9a0ac",opacity=0.5, label=dict(text="mold",textposition="top left"))
        fig.add_shape(type="rect", x0=-underfill_x/2, y0=solder_z+substrate_z, x1=underfill_x/2, y1=solder_z+substrate_z+underfill_z, line=dict( color="black", width=2), fillcolor="#a9a1b3",opacity=0.5,label=dict(text="underfill",textposition="top left"))
        fig.add_shape(type="rect", x0=-bumps_x/2, y0=solder_z+substrate_z, x1=bumps_x/2, y1=solder_z+substrate_z+bumps_z, line=dict( color="black", width=2), fillcolor="#efb46f",opacity=0.5,label=dict(text="bumps",textposition="top left"))
        fig.add_shape(type="rect", x0=-die_x/2, y0=solder_z+substrate_z+bumps_z, x1=die_x/2, y1=solder_z+substrate_z+bumps_z+die_z, line=dict( color="black", width=2), fillcolor="#bbdc77",opacity=0.5,label=dict(text="die",textposition="top left"))
        
        bounds_x = (-mold_x/2-5, mold_x/2+5)
        bounds_y = (0-5, solder_z+substrate_z+mold_z+5)

        fig.update_shapes(dict(xref='x', yref='y'))
        # Set axes properties
        fig.update_xaxes(range=[bounds_x[0], bounds_x[1]], showgrid=False)
        fig.update_yaxes(range=[bounds_y[0], bounds_y[1]])
    return fig

@app.callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks")],
    [State("modal", "is_open")])
def toggle_modal(n1, is_open):
    if n1:
        return not is_open
    return is_open

@app.callback(
    Output("viz_modal", "is_open"),
    [Input("viz_open", "n_clicks"), 
     Input("viz_close", "n_clicks")],
    [State("viz_modal", "is_open")],
)
def viz_toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

if __name__ == "__main__":
    app.run_server(debug=True, port=8051)
