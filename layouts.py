"""
Displays main HTML layout of application
"""

from dash import html
import dash_bootstrap_components as dbc
from views import page1

from app import app

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
		html.Div([html.Img(src=app.get_asset_url(image_filename),
							width="200px", style={'align': 'center'})]),
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


layout1 = page1.layout_1

layout2 = html.Div('Page 2')
layout3 = html.Div('Page 3')
layout4 = html.Div('Page 4')
layout5 = html.Div('Page 5')


if __name__ == '__main__':
    app.run_server(debug=True)
