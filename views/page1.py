from dash import Dash, html, dcc, dash_table
from mmap import PAGESIZE
import plotly.express as px
import pandas as pd
import dash_bootstrap_components as dbc



df = pd.read_csv('demo_data.csv')

fig = px.bar(df, x="name", y="molecular_weight")
# df = pd.DataFrame()
table = dash_table.DataTable(df.to_dict('records'), [{"name": i, "id": i} for i in df.columns])

smile_input = html.Div(
        [
            dbc.Label("SMILE Input"),
            dbc.Input(
                id="input_smiles",
                value='OC(=O)CC1=CC=CC=C1NC1=C(Cl)C=CC=C1Cl',
                type="str",
                style={'color': 'black'}
            ),
            dbc.FormText(
            "Sample SMILE String: OC(=O)CC1=CC=CC=C1NC1=C(Cl)C=CC=C1Cl",
            color="secondary",
        )
        ],
        className="mb-3",
    )
    
form = dbc.FormFloating(
        [
            smile_input
        ]
    )



layout_1 = html.Div([
    dbc.Container([
        html.H1(children='Sample Data'),

        html.Div(table),

        html.Div(children='''
        Dash: A web application framework for your data.
        '''),

        dcc.Graph(
        id='example-graph',
        figure=fig
        ),

        html.Div(children=[
        form,
        dbc.Button(id='submit-smile', n_clicks=0, children='Submit')
            ]),

        html.Br(),

        html.H2('Molecular Structures'),#'SMILES'),
        
        html.Div(
            dbc.Col(dbc.CardImg(
				id='smile-string-image'),
			width="18rem", style={"margin-left": "8rem"})
        ),

        html.H3("Generated Descriptors"),
		html.Div(id='generated-descriptors'),
    
    ])
    
])