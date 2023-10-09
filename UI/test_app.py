import dash
from dash import dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_table
import pandas as pd
import io
import base64
import os

app = dash.Dash(__name__)

app.layout = html.Div([
    html.Label("Choose Folder:"),
    dcc.Input(id='file-input', type='file', multiple=True)
])


@app.callback(
    Output('store-uploaded-data', 'data'),
    Input('file-input', 'value')
)
def store_files(directory_content):
    # Handle the files from the directory here
    # For simplicity, let's assume these are text files and just store their content.
    # For larger files or a more complex directory structure, you'd adjust this accordingly.
    stored_data = {}
    if directory_content:
        for file in directory_content:
            # Convert file content to string
            stored_data[file.filename] = file.read().decode('utf-8')
    return stored_data


@app.callback(
    Output('output-folder-upload', 'children'),
    Input('process-button', 'n_clicks'),
    State('store-uploaded-data', 'data')
)
def update_output(n_clicks, stored_data):
    if n_clicks and stored_data:
        children = [html.Div([html.H5(filename), html.Pre(content)]) for filename, content in stored_data.items()]
        return children


if __name__ == '__main__':
    app.run_server(debug=True, port=8052)
