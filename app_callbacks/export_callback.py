# from app import app, db
# import dash
# from dash import dash_table, html
# from dash.dependencies import Input, Output, State
# import pandas as pd
# from database import connect

# @app.callback(Output('postgres_datatable', 'children'),
#               [Input('interval_pg', 'n_intervals')])
# def populate_datatable(n_intervals):
#     # df = pd.read_sql_table('SampleData', con=db.engine)
#     df = connect.load_data()

#     return [
#         dash_table.DataTable(
#             id='our-table',
#             columns=[{"name": i, "id": i} for i in df.columns],
#             data=df.head().to_dict('records'),
#             # editable=True,
#             # row_deletable=True,
#             filter_action="native",
#             sort_action="native",  # give user capability to sort columns
#             sort_mode="single",  # sort across 'multi' or 'single' columns
#             # page_action='none',  # render all of the data at once. No paging.
#             style_table={'height': 'auto', 'overflowY': 'auto'},
#             style_cell={'textAlign': 'left', 'minWidth': '100px',
#                         'width': '100px', 'maxWidth': '100px'},
#         ),
#     ]

# ## ADD COLUMNS
# @app.callback(
#     Output('our-table', 'columns'),
#     [Input('adding-columns-button', 'n_clicks')],
#     [State('adding-rows-name', 'value'),
#      State('our-table', 'columns')],
#     prevent_initial_call=True)
# def add_columns(n_clicks, value, existing_columns):
#     if n_clicks > 0:
#         existing_columns.append({
#             'name': value, 'id': value,
#             'renamable': True, 'deletable': True
#         })
#     return existing_columns

# ## ADD Rows
# @app.callback(
#     Output('our-table', 'data'),
#     [Input('editing-rows-button', 'n_clicks')],
#     [State('our-table', 'data'),
#      State('our-table', 'columns')],
#     prevent_initial_call=True)
# def add_row(n_clicks, rows, columns):
#     if n_clicks > 0:
#         rows.append({c['id']: '' for c in columns})
#     return rows

# Display graph
# @app.callback(
#     Output('my_graph', 'figure'),
#     [Input('our-table', 'data')],
#     prevent_initial_call=True)
# def display_graph(data):
#     # df_fig = pd.DataFrame(data)
#     # fig = px.bar(df_fig, x='Phone', y='Sales')

#     pg_filtered = db.session.query(Product.Phone, Product.Sales)
#     phone_c = [x.Phone for x in pg_filtered]
#     sales_c = [x.Sales for x in pg_filtered]
#     fig = go.Figure([go.Bar(x=phone_c, y=sales_c)])

#     return fig


# @app.callback(
#     [Output('placeholder', 'children'),
#      Output("store", "data")],
#     [Input('save_to_postgres', 'n_clicks'),
#      Input("interval", "n_intervals")],
#     [State('our-table', 'data'),
#      State('store', 'data')],
#     prevent_initial_call=True)
# def df_to_csv(n_clicks, n_intervals, dataset, s):
#     output = html.Plaintext("The data has been saved to your PostgreSQL database.",
#                             style={'color': 'green', 'font-weight': 'bold', 'font-size': 'large'})
#     no_output = html.Plaintext("", style={'margin': "0px"})

#     input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[
#         0]

#     if input_triggered == "save_to_postgres":
#         s = 6
#         pg = pd.DataFrame(dataset)
#         pg.to_sql("productlist", con=db.engine,
#                   if_exists='replace', index=False)
#         return output, s
#     elif input_triggered == 'interval' and s > 0:
#         s = s - 1
#         if s > 0:
#             return output, s
#         else:
#             return no_output, s
#     elif s == 0:
#         return no_output, s
