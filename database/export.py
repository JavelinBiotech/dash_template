# Not needed
# from app import app
# import os
# import sqlalchemy
# from dotenv import load_dotenv
# import pandas as pd

# from flask_sqlalchemy import SQLAlchemy

# db = SQLAlchemy(app.server)


# class SampleData(db.Model):
#     __tablename__ = 'sample_drug_data'

#     drug_id = db.Column(db.Integer, nullable=False, primary_key=True)
#     name = db.Column(db.VARCHAR(255), nullable = False)
#     canonical_smiles = db.Column(db.Float, nullable=False)
#     molecular_weight = db.Column(db.Float, nullable=False)

#     def __init__(self, drug_id, name, canonical_smiles, molecular_weight):
#         self.drug_id = drug_id
#         self.name = name
#         self.canonical_smiles = canonical_smiles
#         self.molecular_weight = molecular_weight

# def export_data():
# 	# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
# 	# print(f'loading')
# 	# Database connection
# 	#Load Database connections
# 	load_dotenv()
# 	database_url = os.getenv('DATABASE_URL')

# 	# Creating Database engine
# 	engine = sqlalchemy.create_engine(database_url)
	
# 	pass
