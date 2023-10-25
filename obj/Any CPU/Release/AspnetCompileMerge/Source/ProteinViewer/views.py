from datetime import datetime
from email.policy import default
from math import trunc
from pickle import NONE
from unittest.mock import DEFAULT
from flask import Flask, jsonify, make_response, render_template,redirect, request,  Response, url_for, send_from_directory
import os
from matplotlib import gridspec, pyplot
import requests
import tempfile
import shutil
import argparse
import asyncio
import base64
import gzip
from importlib.metadata import files
from io import BytesIO
from multiprocessing import process
import os
import re
import shutil
import tempfile
from tkinter.font import BOLD
import urllib
from urllib.request import urlretrieve
import sys
import aiofiles
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
import requests
from matplotlib import pyplot, gridspec
from sqlalchemy import null
import sys
import json
from io import BytesIO
from flask import send_file
from sqlalchemy import JSON
from ProteinViewer import app
from missense import missense as ms
from os.path import relpath
import requests
import logging
from flask import url_for
import pyodbc
from Uniprot_Items import UniprotItems
from Database import Database
import time
import numpy

CACHE_TIMEOUT = 3600  # Cache timeout in seconds (1 hour)

protein_name_map_cache = {'data': None, 'timestamp': 0}

protein_name_map = Database.get_protein_name_map();
def output_file():
  cwd = os.getcwd() # Construct a relative path to the content folder.
  content_path = os.path.join(cwd, "static", "content")   # Get the current working directory of the application.
  with open(content_path, "w") as f:
    f.write("This is the content of the file.")
  f.close()


def gen_image(pos_to_val):
    positions, targets, values = zip(*pos_to_val)

    positions = numpy.array(positions)
    targets = numpy.array(targets)
    values = numpy.array(values)
    
    unique_to, unique_indices = numpy.unique(targets, return_inverse=True)
    
    img = numpy.zeros((len(unique_to), positions.max() + 1))
    img[unique_indices, positions] = values
    
    return img
def compute_segments(total_length):
    base_increment = 25

    target_segment_size = (total_length // 4) // base_increment * base_increment
    segments = []

    # Calculate the first three segments
    for i in range(3):
        start = i * target_segment_size
        end = (i + 1) * target_segment_size
        segments.append((start, end))

    # Last segment takes the remainder
    segments.append((segments[-1][1], total_length))

    return segments


def make_and_save_plot(pos_to_val, out_file="plot.png"):
    # Determine the total number of positions to visualize
    total_length = max(p[0] for p in pos_to_val)
    x_label_list = numpy.unique([p[1] for p in pos_to_val])[::-1]
    yticks = [y + 0.5 for y in range(len(x_label_list))]
    img = gen_image(pos_to_val)
    median_img =numpy.median(img, axis=0)
    # Accounting for very small proteins that can't fit in a quad structure nicely.
    if total_length <= 100:
        fig, axes = pyplot.subplots(2, 1, figsize=(8, 12))
    
        median_img = numpy.median(img, axis=0)
        # The top plot
        ax = axes[0]
        ax.imshow(img, aspect='auto', interpolation='none', cmap="bwr")
        ax.set_ylim(20, 0)
        ax.set_xlim(0, total_length)
        ax.set_yticks(yticks)
        ax.set_yticklabels(x_label_list, fontsize=13, weight='bold')
        ax.set_ylabel("Alternate amino acid", fontsize=14, weight='bold')
        ax.set_xlabel("Residue sequence number", fontsize=14, weight='bold')
        ax.tick_params(axis='both', which='major', labelsize=12)

        # The bottom plot
        #ax2 = axes[1]
        #ax2.plot(median_img, linewidth=2)
        #ax2.set_ylim(0, 1.1)
        #ax2.set_xlim(0, total_length)
        #ax2.set_ylabel("Mean Pathogenicity", fontsize=14, weight='bold')
        #ax2.set_xlabel("Residue sequence number", fontsize=14, weight='bold')
        #ax2.tick_params(axis='both', which='major', labelsize=13)
    else:
        segments = compute_segments(total_length)

        fig = pyplot.figure(figsize=(16,12))
        gs = gridspec.GridSpec(4, 2, height_ratios=[5, 5, 5,5], width_ratios=[11, 10])
        for idx, (segment_min, segment_max) in enumerate(segments):
            filtered_data = [p for p in pos_to_val if segment_min <= p[0] < segment_max]
            img = gen_image(filtered_data)
            
            ax = fig.add_subplot(gs[idx, 0])
            ax.imshow(img, aspect='auto', interpolation='none', cmap="bwr")
            ax.set_ylim(20, 0)
            ax.set_xlim(segment_min, segment_max)
            ax.set_yticks(yticks)
            ax.set_yticklabels(x_label_list, fontsize=1, weight='bold', rotation="90")
            ax.set_ylabel("Alternate amino acid", fontsize=12, weight='bold')
            ax.set_xlabel("Residue sequence number", fontsize=12, weight='bold')
            ax.tick_params(axis='both', which='major', labelsize=12)
            
            #ax2 = fig.add_subplot(gs[idx, 1])
            #ax2.plot(median_img, linewidth=2)
            #ax2.set_ylim(0, 1.1)
            #ax2.set_xlim(segment_min, segment_max)
            #ax2.set_ylabel("Mean Pathogenicity", fontsize=12, weight='bold')
            #ax2.set_xlabel("Residue sequence number", fontsize=12, fontweight = 'bold')
            #ax.tick_params(axis='both', which='major', labelsize=12)


    pyplot.tight_layout()
    pyplot.savefig(out_file, dpi=300, bbox_inches='tight')  
    img_io = BytesIO()
    pyplot.savefig(img_io, format='png', bbox_inches='tight')
    pyplot.close()
    img_io.seek(0)
    return img_io
     
    if uniprot_Id == None:
        return jsonify("No Protein with that Name or Id could be found double check or try that Protein's Uniprot_Id")
     
    # Assuming you have the PNG file, you can send it as a response
    # You should set the appropriate content type in the response headers
    response = make_response()
@app.route('/fetchImage/<proteinName>', methods=['GET'])
def fetch_image(proteinName):
    try:
        uniprot_Id = Get_Protein_Id(proteinName)
        pos_to_val = Database.get_data_tuples(uniprot_Id)

        if uniprot_Id is None:
            return jsonify("No Protein with that Name or Id could be found double check or try that Protein's Uniprot_Id"), 404
        
        # Make a plot and get it as a BytesIO object
        img = make_and_save_plot(pos_to_val)

        return send_file(img, mimetype='image/png')

    except Exception as e:
        # Handle any unexpected exceptions
        print(f"An error occurred: {e}")
        return jsonify({"error": "An error occurred generating the image."}), 500
        

def create_pathogenic_pdb(img, uniprot_Id):

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_Id.upper()}"
    response = requests.get(api_url, timeout=20)
    response.raise_for_status()  # Raise exception for bad requests

    r = response.json()
    target_pdb_data = requests.get(r[0]['pdbUrl']).content

    mean_per_pos = np.median(img, axis=0)

    # Create a BytesIO buffer to hold the modified PDB data
    pdb_buffer = BytesIO(target_pdb_data)
    modified_pdb_buffer = BytesIO()

    for line in pdb_buffer.read().decode('utf-8').splitlines():
        if line.startswith("ATOM "):
            pos = int(line[22:26]) - 1 
            value_str = f"{mean_per_pos[pos]:.2f}".rjust(6) 
            edit_line = f"{line[:60]}{value_str}{line[67:]}"
            modified_pdb_buffer.write(edit_line.encode() + b'\n')
        else:
            modified_pdb_buffer.write(line.encode() + b'\n')

    # Reset buffer position to the beginning for subsequent reads
    modified_pdb_buffer.seek(0)

    return modified_pdb_buffer

@app.route('/fetchPDB/<proteinName>', methods=['GET'])
def fetch_pdb(proteinName):
    try:
        uniprot_Id = Get_Protein_Id(proteinName)
        if uniprot_Id == None:
            return jsonify("No Protein with that Name or Id could be found. Double check or try that Protein's Uniprot_Id")
        
        pos_to_val = Database.get_data_tuples(uniprot_Id)
        img = gen_image(pos_to_val)
    
        
        PDBFile = create_pathogenic_pdb(img, uniprot_Id)
        
        # Assuming you have the PDB file, you can send it as a response
        # You should set the appropriate content type in the response headers
        response = make_response(PDBFile.getvalue())
        response.headers['Content-Type'] = 'chemical/x-pdb'
        return response

    except Exception as e:
        try:
            protein_id = UniprotItems.verify_protein_id(proteinName)
            return generate_response(protein_id)
        except:
            return jsonify({'Searches for valid Uniprot_Id, Protein_Id and Protein_Name failed'}), 500
def Get_Protein_Id(proteinName):
    uniprot_Id = None
    current_time = time.time()
    
    try:
        # Check if cache is valid
        if (protein_name_map_cache.get('data') and  (current_time - protein_name_map_cache.get('timestamp', 0)) < CACHE_TIMEOUT):
            Protein_Name_Map = protein_name_map_cache['data']
        else:
            # Refresh cache
            Protein_Name_Map = Database.get_protein_name_map()
            protein_name_map_cache['data'] = Protein_Name_Map
            protein_name_map_cache['timestamp'] = current_time
        
        # Attempt to get UniProt ID using protein name
        uniprot_Id = Protein_Name_Map.get(proteinName)

        # Check if proteinName starts with "ENSG0" - assuming it might be an Ensembl ID
        if not uniprot_Id and proteinName.startswith("ENSG0"):
            uniprot_Id = UniprotItems.get_uniprot_id_from_ensemble(proteinName)
        elif not uniprot_Id and proteinName in set(Protein_Name_Map.values()):
            holder =  [key for key, value in Protein_Name_Map.items() if value == proteinName][0] #someone searched usinga UniprotId
            uniprot_Id = Protein_Name_Map.get(holder)
        elif not uniprot_Id:
            raise ValueError("No match found for protein name or ID")
        
        return uniprot_Id
    
    except Exception as e:
        # Handle exception (e.g., logging)
        print(f"An error occurred: {str(e)}")
        return None

    except Exception as e:
        # Log the exception for debugging
        print(f"An error occurred: {str(e)}")
        return None
def Get_PNG(Uniprot_Id):
    try:
        pos_to_vals = Database.get_data_tuples(Uniprot_Id)
        finalImage = make_and_save_plot(pos_to_vals, "Plot.png")
        return send_file(finalImage, mimetype='image/png')
    except:
        return None
def Get_PDB(Uniprot_Id):
    try:
        pos_to_vals = Database.get_data_tuples(Uniprot_Id)
        intermediate_image = ms.gen_image(pos_to_vals)
        pdb_data = ms.create_pathogenic_pdb(intermediate_image, Uniprot_Id)
        return pdb_data
    except:
        return None

    

@app.route('/api/protein_name_map')
async def protein_name_map():
    data = await Database.get_protein_name_map_async()
    return jsonify(data)

@app.route('/fetchIndividual/<proteinName>', methods=['GET'])
def fetchIndividualProteinMap(proteinName):
    try:
        uniprot_Id = Get_Protein_Id(proteinName)
        data = Database.get_median_pathogenicity_data(uniprot_Id)
  
        return jsonify(data)
    except Exception as e:
        return jsonify({'Searches for valid Uniprot_Id, Protein_Id and Protein_Name failed'}), 500


def fetch_image(proteinName):
    try:
        uniprot_Id = Get_Protein_Id(proteinName)
        pos_to_val = ... # somehow get or calculate your pos_to_val

        if uniprot_Id is None:
            return jsonify("No Protein with that Name or Id could be found double check or try that Protein's Uniprot_Id"), 404
        
        # Make a plot and get it as a BytesIO object
        img = make_and_save_plot(pos_to_val)

        return send_file(img, mimetype='image/png')

    except Exception as e:
        # Handle any unexpected exceptions
        print(f"An error occurred: {e}")
        return jsonify({"error": "An error occurred generating the image."}), 500

def generate_response(protein_id):
    try:
        # Construct the paths based on the protein_id
        content_path = os.path.join(app.static_folder, "content")
        output_path = os.path.join(content_path, protein_id)
        tsvpath = os.path.join(content_path, "alpha.tsv")
        pdbpath = os.path.join(output_path, f'{protein_id}-edit.pdb')

        # Ensure output directory exists
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        # Run Pymissense function
        maxacid = 200  # Or any default value suitable
        ms._run(protein_id, output_path, tsvpath, pdbpath, maxacid)

        # Construct the URL paths for the generated files
        image_url = url_for('static', filename=f'content/{protein_id}/{protein_id}.pdf')
        pdb_url = url_for('static', filename=f'content/{protein_id}/{protein_id}-edit.pdb')

        return jsonify({
            'pdbPath': pdb_url,
            'imgPath': image_url
        })

    except Exception as e:
        # Handle any errors that may arise
        print(e)  # You can also use app.logger.error(e) for better logging
        return jsonify({'error': 'An error occurred while processing the protein.'}), 500
@app.route('/', methods=['GET', 'POST'])
def index():
    db =  Database()
    message = None
    protein_name_map_cache['data'] = Database.get_protein_name_map()
    protein_name_map_cache['timestamp'] = time.time()
    protein_name_map = protein_name_map_cache['data']
    Protein_Name = "SARM1"
    pdb_file_path = ""
    img_file_path = ""
    maxacid = 200  # Or any default value suitable
    message = "Data Sent correct"
    Uniprot_ID = "Q6SZW1"
    data = db.get_protein_data(Uniprot_ID)
    if data:
        Uniprot_ID = data["Uniprot_Id"]
        sequence= data["Sequence"]
        length = data["Length"]
        maxacid = data["Length"]  if data["Length"] == None else int(data["Length"])   
        message = "No data found for the provided protein name/Ensemble ID."
    ProteinName = "SARM1"
    #content_path = os.path.join(os.getcwd(), "static", "content")

    #tsvpath = content_path + "\\alpha.tsv"
    #pdbpath = None

    #ms._run(proteinName, output_path, tsvpath, content_path, maxacid)
    
     
    return render_template('index.html', Protein_Name_Map= protein_name_map, ProteinList= protein_name_map_cache['data'].keys, Uniprot_Id = Uniprot_ID, message = message)

if __name__ == "__main__":
    app.run(debug=True)
@app.route('/contact')
def contact():
    """Renders the contact page."""
    return render_template(
        'contact.html',
        title='Contact',
        year=datetime.now().year,
        message='Your contact page.'
    )

@app.route('/about')
def about():
    """Renders the about page."""
    return render_template(
        'about.html',
        title='About',
        year=datetime.now().year,
        message='Your application description page.'
    )

