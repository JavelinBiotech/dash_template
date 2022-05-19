import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors



def generate_df(smile):
    data = {'smiles': [smile],
           }
    df = pd.DataFrame(data)
    return df

def generate_mol(df):
    mol = df.smiles.apply(Chem.MolFromSmiles)
    #return mole and smile string as a dataframe -> maybe even draw it out

    return mol

def generate_2Ddescriptors(df, mol):

    df['MolLogP'] = mol.apply(Descriptors.MolLogP)
    df['HeavyAtomCount'] = mol.apply(Descriptors.HeavyAtomCount)
    df['HAccept'] = mol.apply(Descriptors.NumHAcceptors)
    df['Heteroatoms'] = mol.apply(Descriptors.NumHeteroatoms)
    df['HDonor'] = mol.apply(Descriptors.NumHDonors)
    df['MolWt'] = mol.apply(Descriptors.MolWt)
    df['RotableBonds'] = mol.apply(Descriptors.NumRotatableBonds)
    df['RingCount'] = mol.apply(Descriptors.RingCount)
    df['Ipc'] = mol.apply(Descriptors.Ipc)
    df['HallKierAlpha'] = mol.apply(Descriptors.HallKierAlpha)
    df['NumValenceElectrons'] = mol.apply(Descriptors.NumValenceElectrons)
    df['SaturatedRings'] = mol.apply(Descriptors.NumSaturatedRings)
    df['AliphaticRings'] = mol.apply(Descriptors.NumAliphaticRings)
    df['AromaticRings'] = mol.apply(Descriptors.NumAromaticRings)

    # Part of top 10 descriptors from Wang et al
    df['PEOE_VSA6'] = mol.apply(Descriptors.PEOE_VSA6 )  # Total negative 6vdw surface area #2
    df['SlogP_VSA10'] = mol.apply(Descriptors.SlogP_VSA10) # Bin 0 SlogP (-10, -0.40)  #3
    df['SlogP_VSA2'] = mol.apply(Descriptors.SlogP_VSA2 ) # Bin 2 SlogP (-0.20, 0.00]) # 6

    # Subdivided Surgace Areas
    df['SlogP_VSA1'] = mol.apply(Descriptors.SlogP_VSA1)
    df['SlogP_VSA2'] = mol.apply(Descriptors.SlogP_VSA2)
    df['SlogP_VSA3'] = mol.apply(Descriptors.SlogP_VSA3)
    df['SlogP_VSA4'] = mol.apply(Descriptors.SlogP_VSA4)
    df['SlogP_VSA7'] = mol.apply(Descriptors.SlogP_VSA7)
    df['SlogP_VSA9'] = mol.apply(Descriptors.SlogP_VSA9)
    df['SMR_VSA1'] = mol.apply(Descriptors.SMR_VSA1)
    df['SMR_VSA2'] = mol.apply(Descriptors.SMR_VSA2)
    df['SMR_VSA3'] = mol.apply(Descriptors.SMR_VSA3)
    df['SMR_VSA5'] = mol.apply(Descriptors.SMR_VSA5)

    # Partial Charge Descriptors, Not the same as Paper, doesn't denonte charge
    df['PEOE_VSA1'] = mol.apply(Descriptors.PEOE_VSA1)
    df['PEOE_VSA2'] = mol.apply(Descriptors.PEOE_VSA2)
    df['PEOE_VSA4'] = mol.apply(Descriptors.PEOE_VSA4)
    df['PEOE_VSA5'] = mol.apply(Descriptors.PEOE_VSA5)

    # Fingerprint Density
    df['FpDensityMorgan1'] = mol.apply(Descriptors.FpDensityMorgan1)
    df['FpDensityMorgan2'] = mol.apply(Descriptors.FpDensityMorgan2)
    df['FpDensityMorgan3'] = mol.apply(Descriptors.FpDensityMorgan3)

    df['MinPartialCharge'] = mol.apply(Descriptors.MaxPartialCharge)
    df['MaxAbsPartialCharge'] = mol.apply(Descriptors.MinPartialCharge)
    df['MaxPartialCharge'] = mol.apply(Descriptors.MaxAbsPartialCharge)

    return df