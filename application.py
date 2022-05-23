"""
Instantiates the Dash app and identify the server
https://towardsdatascience.com/callbacks-layouts-bootstrap-how-to-create-dashboards-in-plotly-dash-1d233ff63e30
"""
from dash import dcc, html
from dash.dependencies import Input, Output

from app import app, server
# import your navigation, styles and layouts from layouts.py here
from layouts import nav_bar, layout1, layout2, layout3, layout4, layout5, CONTENT_STYLE
import callbacks

# Define basic structure of app:
# A horizontal navigation bar on the left side with page content on the right.
app.layout = html.Div([
    dcc.Location(id='url', refresh=False), #this locates this structure to the url
    nav_bar(),
    html.Div(id='page-content',style=CONTENT_STYLE) #we'll use a callback to change the layout of this section 
])

# This callback changes the layout of the page based on the URL
# For each layout read the current URL page "http://127.0.0.1:5000/pagename" and return the layout
@app.callback(Output('page-content', 'children'), #this changes the content
              [Input('url', 'pathname')]) #this listens for the url in use
def display_page(pathname):
    if pathname == '/':
        return layout1
    elif pathname == '/page1':
        return layout1
    elif pathname == '/page2':
         return layout2
    elif pathname == '/page3':
         return layout3
    elif pathname == '/page4':
         return layout4
    elif pathname == '/page5':
        return layout5
    else:
        return '404' #If page not found return 404

#Runs the server at http://127.0.0.1:5000/
if __name__ == '__main__':
    app.run_server(port=5000, host= '127.0.0.1',debug=True)