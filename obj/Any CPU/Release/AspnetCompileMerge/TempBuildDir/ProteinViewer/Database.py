from argparse import FileType
import base64
import json
import os
import asyncio
from pickle import TRUE
from numpy import median
import numpy
import pandas as pd
import pyodbc
from sqlalchemy import false, true, values
from sqlalchemy.engine import row
from urllib3 import connection
import fast_to_sql
from Uniprot_Items import UniprotItems
import numpy
class Database:
    @staticmethod
    def _get_connection_string_from_secrets(file_path):
        with open(file_path, 'r') as file:
            data = json.load(file)
            db = data['database']
            
            connection_string = (f"Driver={{{db['driver']}}};"
                                 f"Server={db['server']};"
                                 f"Database={db['database_name']};"
                                 f"Uid={db['uid']};"
                                 f"Pwd={db['pwd']};"
                                 f"Encrypt={db['encrypt']};"
                                 f"TrustServerCertificate={db['trust_server_certificate']};"
                                 f"Connection Timeout={db['connection_timeout']};")
            return connection_string
    
    @classmethod
    def get_connection(cls):
        path = os.path.join(os.getcwd(), "static", "content")
        file_path = os.path.join(path, "secrets.json")
        connection_string = cls._get_connection_string_from_secrets(file_path)
        return pyodbc.connect(connection_string)
    @staticmethod
    def _sync_insert(uniprot_id, img_encoded, pdb_encoded, connection_string):
        # Your usual connection string for SQL Server
        conn_str = connection_string;
        with pyodbc.connect(conn_str) as conn:
            with conn.cursor() as cursor:
                # Construct your INSERT SQL statement here. This is just a placeholder.
                sql = """INSERT INTO YourTableName (UniProtID, ImageData, PDBData)
                         VALUES (?, ?, ?)"""
                cursor.execute(sql, uniprot_id, img_encoded, pdb_encoded)
                conn.commit()

    @staticmethod
    async def insert_into_db_async(uniprot_id, img_encoded, pdb_encoded):
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(None, Database._sync_insert, uniprot_id, img_encoded, pdb_encoded)
    @staticmethod
    def fast_insert(df, table_name, if_exists="append"):
        connection = Database.get_connection()
        fast_to_sql.frame_to_sql(df, table_name, connection, if_exists=if_exists)
        connection.close()
      

    def File_from_db(connection_str, uniprot_id, output_directory):
        with pyodbc.connect(connection_str) as conn:
            cursor = conn.cursor()
            query = """SELECT [Generated_Image], [Generated_Model], [FileName], [FileType] 
                       FROM [ProteinData].[ViewableData] 
                       WHERE [Uniprot_Id] = ?"""
            cursor.execute(query, uniprot_id)
            row = cursor.fetchone()
        
            # Decode and save to the respective file type
            if row:
                file_data = row.Generated_Image if row.FileType == 'png' else row.Generated_Model
                decoded_data = base64.b64decode(file_data)
            
                output_path = os.path.join(output_directory, row.FileName)
                with open(output_path, "wb") as f:
                    f.write(decoded_data)
    
    def insert_into_db(uniprot_id, encoded_data, file_type):
        with Database.get_connection() as conn:
            cursor = conn.cursor()
            query = """INSERT INTO [ProteinData].[ViewableData] 
                   ([Uniprot_Id], [File], [FileName], [FileType]) 
                   VALUES (?, ?, ?, ?)"""
            Name = uniprot_id + file_type
            cursor.execute(query, uniprot_id, encoded_data, Name, file_type)
                   
            conn.commit()
    @staticmethod
    def execute_non_query(sql):
        connection = Database.get_connection()
        cursor = connection.cursor()
        cursor.execute(sql)
        connection.commit()
        connection.close()

    @staticmethod
    def execute_query(sql, dataframe = false):
     connection = Database.get_connection()
     if(dataframe):
         result = pd.read_sql(sql, connection)
     else:
         result = connection.execute(sql)
     return result
 
    @staticmethod
    def get_uniprot_list():
        sql = "SELECT Distinct Uniprot_id from ProteinData.Proteins" 
        cursor = Database.execute_query(sql, False)
        return [row[0] for row in cursor.fetchall()]


    @staticmethod
    def get_data_tuples(uniprot_id):
        sql = f"""
        SELECT Alternate, Position, PathogenicityScore
        FROM ProteinData.Pathogenicity 
        WHERE Uniprot_Id = '{uniprot_id}'
    """
    
       
        cursor = Database.execute_query(sql, False)
        pos_to_val = []

        for row in cursor.fetchall():
            pos_to_val.append((row[1], row[0], row[2]))
        return pos_to_val
   

    
      # This captures the target amino acid of the mutation.

        val = round(float(row['PathogenicityScore']), 4)

        return pos_to_val

    @staticmethod
    def get_protein_data(uniprot_Id):
        sql = f"SELECT Uniprot_Id, ProteinName, Ensembl_Id, Sequence, Length FROM ProteinData.Proteins WHERE Uniprot_Id='{uniprot_Id}'"
        df = Database.execute_query(sql)
        if not df.empty:
            row = df.iloc[0]
            return {
                "Uniprot_Id": row["Uniprot_Id"],
                "ProteinName": row["ProteinName"],
                "Ensembl_Id": row["Ensembl_Id"],
                "Sequence": row["Sequence"],
                "Length": row["Length"]
            }
        return None

    @staticmethod
    def get_protein_name_map():
        sql = "SELECT ProteinName, Uniprot_Id FROM ProteinData.Proteins"
        temps = Database.execute_query(sql, False)

        return dict(temps.fetchall())
        

        connection.close()
        return result
 
    @staticmethod
    def get_median_pathogenicity_data(uniprot_id):
        # Complete set of amino acids
        ALL_AMINO_ACIDS = set(['C', 'F', 'I', 'N', 'R', 'V', 'D', 'G', 'K', 'P', 'S', 'W', 'A', 'E', 'H', 'L', 'Q' ,'T', 'Y', 'M'])

        # Query to fetch data from the Pathogenicity table
        sql = f"""
                  SELECT PA.Position, PA.Alternate, PA.PathogenicityScore
                  FROM ProteinData.Pathogenicity PA
                  WHERE PA.Uniprot_Id = '{uniprot_id}' ORDER BY PA.Position"""
     
        data = Database.execute_query(sql, False).fetchall()
    
        # Dictionary to hold amino acids and their scores by position
        amino_data_by_position = {}

        for position, alternate, score in data:
            if position not in amino_data_by_position:
                amino_data_by_position[position] = {
                    'alternates': set(),
                    'scores': {}
                }
        
            amino_data_by_position[position]['alternates'].add(alternate)
            amino_data_by_position[position]['scores'][alternate] = score

        # Identifying the original amino acid for each position
        protein_sequence = ""
        amino_scores = {}
        for position, data in sorted(amino_data_by_position.items()):
            missing_amino = list(ALL_AMINO_ACIDS - data['alternates'])[0]  # The missing amino acid from the set of 20
            protein_sequence += missing_amino
            amino_scores[position] = data['scores']

        return protein_sequence, amino_scores
