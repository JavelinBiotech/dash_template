from distutils.command.upload import upload
from dash import html, dash_table, dcc
import pandas as pd
import dash_bootstrap_components as dbc

#Connect Database
from database import connect

df = connect.load_data()

table = dash_table.DataTable(df.head().to_dict('records'), [{"name": i, "id": i} for i in df.columns],
                             style_table={'height': 'auto', 'overflowY': 'auto'})

upload_card = dcc.Upload(
    id='upload-data',
    children=html.Div([
        'Drag and Drop or ',
        html.A('Select Files')
    ]),
    style={
        'width': '100%',
        'height': '60px',
        'lineHeight': '60px',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'borderRadius': '5px',
        'textAlign': 'center',
        'margin': '10px'
    },
    # Allow multiple files to be uploaded
    multiple=True
),


layout_2 = html.Div(['Hello',
                    html.Div(upload_card),
                    html.Div(id='output-data-upload'),

                    dbc.Button('Save to PostgreSQL',
                                id='save_to_postgres', n_clicks=0),

                    html.Div(id='output-div'),

 ])