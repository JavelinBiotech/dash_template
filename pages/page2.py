from distutils.command.upload import upload
from dash import html, dash_table, dcc
import pandas as pd
import dash_bootstrap_components as dbc

#Connect Database
from database import connect, export

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

# received = html.Div([
#     # html.Div([
#     #     dcc.Input(
#     #         id='adding-rows-name',
#     #         placeholder='Enter a column name...',
#     #         value='',
#     #         style={'padding': 10}
#     #     ),
#     #     html.Button('Add Column', id='adding-columns-button', n_clicks=0)
#     # ], style={'height': 50}),

#     # activated once/week or when page refreshed
#     dcc.Interval(id='interval_pg', interval=86400000*7, n_intervals=0),
#     html.Div(id='postgres_datatable'),

#     # # html.Button('Add Row', id='editing-rows-button', n_clicks=0),
#     # html.Button('Save to PostgreSQL', id='save_to_postgres', n_clicks=0),

#     # Create notification when saving to excel
#     html.Div(id='placeholder', children=[]),
#     # dcc.Store(id="store", data=0),
#     # dcc.Store(id="store", data=0),

#     dcc.Interval(id='interval', interval=1000),

#     # Graph results
#     # dcc.Graph(id='my_graph')

# ])

layout_2 = html.Div(['Hello',
                    html.Div(upload_card),
                    html.Div(id='output-data-upload'),
                    # dcc.Store(id="stored-data", data=None), # placeholder # need to handle error
                    dbc.Button('Save to PostgreSQL',
                                id='save_to_postgres', n_clicks=0),
                    # html.Div(table),
                    html.Div(id='output-div'),
                    # html.Div(received)
 ])