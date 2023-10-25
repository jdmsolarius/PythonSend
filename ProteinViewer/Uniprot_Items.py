import os
import requests
from flask import jsonify, url_for

class UniprotItems: 

    @staticmethod
    def get_uniprot_number_from_protein_name(protein_name="sarm1"):
        return UniprotItems._fetch_primary_accession(
            f"https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(organism_id:9606)%20AND%20(protein_name={protein_name})"
        )
    
    @staticmethod
    def verify_protein_id(protein_id):
        return UniprotItems._fetch_primary_accession(
            f"https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(organism_id:9606)%20AND%20(protein_id={protein_id})"
        )

    @staticmethod
    def _fetch_primary_accession(url):
        headers = {"Accept": "application/json"}
        response = requests.get(url, headers=headers)
        
        if response.status_code != 200:
            raise Exception(f"Failed to fetch PDB entries for URL {url}")
        
        return response.json()["results"][0].get("primaryAccession")

    @staticmethod
    def fetch_pdb_and_image(proteinName):
        try:
            protein_id = UniprotItems.get_uniprot_id(proteinName)
            return UniprotItems._generate_response(protein_id)

        except Exception as e:
            protein_id = UniprotItems.verify_protein_id(proteinName)
            if protein_id:
                return UniprotItems._generate_response(protein_id)
            return jsonify({'searchs for valid Uniprot_Id, Protein_Id and Protein_Name failed'}), 500

    @staticmethod
    def get_uniprot_id(proteinName):
        if proteinName.startswith("ENS"):
            return UniprotItems.get_uniprot_id_from_ensemble(proteinName)
        return UniprotItems.get_uniprot_number_from_protein_name(proteinName)

    @staticmethod
    def _generate_response(protein_id):
        # Define app.static_folder, url_for, and ms._run in your context or import them
        content_path = os.path.join(app.static_folder, "content")
        output_path = os.path.join(content_path, protein_id)

        if not os.path.exists(output_path):
            tsv_path = os.path.join(content_path, "alpha.tsv")
            ms._run(protein_id, app.static_folder, tsv_path, None, 200)

        image_url = url_for('static', filename=f'content/{protein_id}/{protein_id}.pdf')
        pdb_url = url_for('static', filename=f'content/{protein_id}/{protein_id}-edit.pdb')
        
        return jsonify({
            'pdbPath': pdb_url,
            'imgPath': image_url
        })

    @staticmethod
    def get_uniprot_id_from_ensemble(ensembl_id):
        url = 'https://rest.uniprot.org/idmapping/run'
        params = {
            'ids': ensembl_id,
            'from': 'Ensembl',
            'to': 'UniProtKB'
        }
        response = requests.post(url, params=params)
        if response.status_code == 200:
            data = response.json()
            return data['results'][0]['to']
        raise Exception('Failed to get UniProt ID for Ensembl ID {}'.format(ensembl_id))