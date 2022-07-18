import os
import sqlalchemy
from dotenv import load_dotenv
import pandas as pd


def load_data():
	# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

	# Database connection
	#Load Database connections
    load_dotenv()
    database_url = os.getenv('DATABASE_URL')

	# Creating Database engine
    engine = sqlalchemy.create_engine(database_url)

	# assume you have a "long-form" data frame
	# see https://plotly.com/python/px-arguments/ for more options

    df = pd.read_sql(
    '''SELECT
	drug.drug_id, name AS "Name",
	cl_value AS "Clearance",
    vd_value AS "Volume of Distribution", molecular_weight, molecular_species, cx_logp,
    hba AS "HBA", hbd AS "HBD",
	polar_surface_area, rotatable_bonds, lipinski_ro5_violations,
	canonical_smiles
	FROM molecular_property
	JOIN drug
	ON drug.drug_id = molecular_property.drug_id
	JOIN clinical_clearance
	ON molecular_property.drug_id = clinical_clearance.drug_id
	JOIN clinical_vd
	ON molecular_property.drug_id = clinical_vd.drug_id
	JOIN reference
	ON clinical_vd.ref_id = reference.ref_id


	WHERE
		(clinical_clearance.cl_type = 'total'
		AND clinical_clearance.condition = 'control')
		AND ("clinical_vd"."condition" = 'control')
		AND (clinical_vd.vd_type = 'steady-state' OR clinical_vd.vd_type = 'apparent')
		AND (molecular_species, molecular_weight, cx_logp, HBA, HBD,
		polar_surface_area, rotatable_bonds, lipinski_ro5_violations, canonical_smiles) IS NOT NULL
		'''
		, con=engine)

	#close connection
    engine.dispose()
    return df
