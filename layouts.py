from dash import Dash, html, dcc, dash_table
from mmap import PAGESIZE
import plotly.express as px
import pandas as pd
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
# import callbacks

from app import app

# app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

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

# Styles & Colors
##########################################################################

NAVBAR_STYLE = {
	"position": "fixed",
	"top": 0,
	"left": 0,
	"bottom": 0,
	"width": "16rem",
	"padding": "2rem 1rem",
	"background-color": "#f8f9fa",
}

CONTENT_STYLE = {
	"top":0,
	"margin-top":'2rem',
	"margin-left": "18rem",
	"margin-right": "2rem",
}

##########################################################################
# Create Auxiliary Components Here
##########################################################################

# Javelin Logo
image_filename = 'Javelin-Biotech_LOGO.png'

def nav_bar():
	"""
	Creates Navigation bar
	"""
	navbar = html.Div(
	[
		#html.H4("Javelin Apps", className="display-10",style={'textAlign':'center'}),
		#html.Div([dcc.Graph(figure=fig_javelin_logo)]),
		#html.Img(src=image_filename, height="30px"),
		#html.Div([html.Img(src='data:image/png;base64,{}'.format(encoded_logo))]),
		html.Div([html.Img(src=app.get_asset_url(image_filename), width="200px", style={'align': 'center'})]),
		html.Hr(),
		dbc.Nav(
			[
				dbc.NavLink("Sample Page 1", href="/page1",active="exact", external_link=True),
				dbc.NavLink("Sample Page 2", href="/page2", active="exact", external_link=True),
				dbc.NavLink("Sample Page 3", href="/page5", active="exact", external_link=True),
				dbc.NavLink("Sample Page 4", href="/page4", active="exact", external_link=True),
				dbc.NavLink("Sample Page 5", href="/page3", active="exact", external_link=True),

			],
			pills=True,
			vertical=True
		),
	],
	style=NAVBAR_STYLE,
	)  
	return navbar



layout1 = html.Div([
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

layout2 = html.Div('Page 2')
layout3 = html.Div('Page 3')
layout4 = html.Div('Page 4')
layout5 = html.Div('Page 5')


if __name__ == '__main__':
    app.run_server(debug=True)

