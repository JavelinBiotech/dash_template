from app import app
from dash import dash_table
from dash.dependencies import Input, Output, State

#Used to draw Structure
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
import base64
from io import BytesIO
import PIL

# Generates 2D Descriptors
import analysis

# Single STRING Input Generation
@app.callback(
    # Output('user-input-output', 'children'), # THIS CAN BE WHERE THE MODEL GOES
    Output('generated-descriptors', 'children'),
    Output('smile-string-image', 'src'),
    Input('submit-smile', 'n_clicks'),
    State('input_smiles', 'value'),
)
def show_smile_string(n_clicks, value):
    # print(n_clicks)

    smile_df = analysis.generate_df(value)
    mol = analysis.generate_mol(smile_df)    
    # print(smile_df)

    result_list = [Chem.MolFromSmiles(smiles) for smiles in smile_df.smiles]
    img = MolsToGridImage(result_list)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    encoded_image = base64.b64encode(buffered.getvalue())
    src_str = 'data:image/png;base64,{}'.format(encoded_image.decode())

    complete_df = analysis.generate_2Ddescriptors(smile_df, mol)
    # Can add condition if user != Javelin then rename the columns otherwise, pass
    complete_df.columns = [f'Generated Descriptor: {i}' for i in range(complete_df.shape[1])]
    des_table = dash_table.DataTable(data=complete_df.to_dict('records'), 
                                columns=[{'name': str(i), 'id': str(i)} for i in complete_df.columns],
								style_cell={'whiteSpace': 'normal','height': 'auto','width': 'auto'},
								style_table={'height': 'auto', 'overflowY': 'auto'}
								)

    # MODEL on AWS predict -> send smile string to AWS,
    # AWS loads model and generates prediction, predicted result is sent back to the application

    return des_table, src_str
