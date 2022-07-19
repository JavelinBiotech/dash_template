from app import app, db
import base64
import datetime
import io

from dash import html, dash_table, dcc
from dash.dependencies import Input, Output, State

import pandas as pd

def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    # return filename, date, df
    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        dash_table.DataTable(
            df.to_dict('records'),
            [{'name': i, 'id': i} for i in df.columns], 
            style_table={'height': 'auto', 'overflowY': 'auto'}
        ),

        html.Hr(),  # horizontal line

        dcc.Store(id='stored-data', data=df.to_dict('records'))

        # For debugging, display the raw contents provided by the web browser
        # html.Div('Raw Content'),
        # html.Pre(contents[0:200] + '...', style={
        #     'whiteSpace': 'pre-wrap',
        #     'wordBreak': 'break-all'
        # })
    ])

@app.callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              [State('upload-data', 'filename'),
              State('upload-data', 'last_modified')])
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
       
        return children

@app.callback([Output('output-div', 'children')],
              [Input('save_to_postgres', 'n_clicks')],
              [State('stored-data', 'data')]
)
def df_to_db(n_clicks, data):
    if n_clicks > 0:
        print('clicked')
        # print(data)
        df = pd.DataFrame(data)
        df.to_sql("sample_drug_data", con=db.engine,
                                  if_exists='replace', index=False)
        print('Success!')
        
        updated = pd.read_sql_table('sample_drug_data', con=db.engine)
        
        return [html.Div([
            html.H2('Updated Table from Database'), 
            dash_table.DataTable(
                updated.to_dict('records'),
                [{'name': i, 'id': i} for i in updated.columns],
                style_table={'height': 'auto', 'overflowY': 'auto'}
            )
        ])]
