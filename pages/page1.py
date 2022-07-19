from dash import html, dash_table, dcc
import plotly.express as px
import pandas as pd
import dash_bootstrap_components as dbc

df = pd.read_csv('demo_data.csv')

fig = px.bar(df, x="name", y="molecular_weight")
# df = pd.DataFrame()
table = dash_table.DataTable(df.to_dict('records'), [{"name": i, "id": i} for i in df.columns],
                                style_table={'height': 'auto', 'overflowY': 'auto'})
#Div for smile input
smile_input = html.Div(
        [
            dbc.Label("SMILE Input"),
            dbc.Input(
                id="input_smiles",
                value='O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl',
                type="str",
                style={'color': 'black'}
            ),
            dbc.FormText(
            "Sample SMILE String; Diclofenac: O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
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

        html.Br(),

        html.H2(children='Molecular Weights'),

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