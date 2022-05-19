from dash import Dash, html, dcc, dash_table
import plotly.express as px
import pandas as pd
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import callbacks

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
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


app.layout = html.Div([
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
    # form,

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



# CALL BACKS
from dash import dcc
import dash_bootstrap_components as dbc
import dash as html
from dash import dash_table
from dash.dependencies import Input, Output, State

# RDKIT to Draw Smiles
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
import base64
from io import BytesIO
import numpy as np 
import PIL

import analysis

# Single STRING and QSAR Prediction Input Generation
@app.callback(
    # Output('user-input-output', 'children'),
    Output('generated-descriptors', 'children'),
    Output('smile-string-image', 'src'),
    Input('submit-smile', 'n_clicks'),
    State('input_smiles', 'value'),
)
def show_smile_string(n_clicks, value):
    print(n_clicks)

    smile_df = analysis.generate_df(value)
    mol = analysis.generate_mol(smile_df)    
    print(smile_df)
    result_list = [Chem.MolFromSmiles(smiles) for smiles in smile_df.smiles]
    img = MolsToGridImage(result_list)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    encoded_image = base64.b64encode(buffered.getvalue())
    src_str = 'data:image/png;base64,{}'.format(encoded_image.decode())



    complete_df = analysis.generate_2Ddescriptors(smile_df, mol)
    # Can add condition if user != Javelin then rename the columns otherwise, pass
    complete_df.columns = [f'Generated Descriptor: {i}' for i in range(complete_df.shape[1])]
    des_table = dash_table.DataTable(data=complete_df.to_dict('records'), columns=[{'name': str(i), 'id': str(i)} for i in complete_df.columns],
								style_cell={'whiteSpace': 'normal','height': 'auto','width': 'auto'},
								style_table={'height': 'auto', 'overflowY': 'auto'}
								)


    # from model import model
    # model.predict(value)


    return des_table, src_str

if __name__ == '__main__':
    app.run_server(debug=True)

