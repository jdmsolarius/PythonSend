from fast_to_sql import fast_to_sql as fs
class Protein:
    def __init__(self, uniprot_id, ensemble_id=None, protein_name=None, sequence=None, length=None, insert_date=None):
        self.id = None
        self.uniprot_id = uniprot_id
        self.protein_name = protein_name
        self.ensemble_id = ensemble_id
        self.sequence = sequence
        self.length = length
        self.insert_date = insert_date
        self.pathogenicities = []

    def add_pathogenicity(self, pathogenicity):
        """Add a ProteinPathogenicity instance to the Protein."""
        if isinstance(pathogenicity, ProteinPathogenicity):
            self.pathogenicities.append(pathogenicity)
        else:
            raise ValueError("Only ProteinPathogenicity instances can be added.")
    
    def generate_protein_insert(self):
        """Generate SQL INSERT statement for the Proteins table."""
        sql = f"""INSERT INTO ProteinData.Proteins (Uniprot_Id, Ensemble_Id, ProteinName, Sequence, Length)
                  VALUES ('{self.uniprot_id}', '{self.ensemble_id}', '{self.protein_name}', '{self.sequence}', {self.length});"""
        return sql
    def Create_Protein_Map(self):
        """Generate SQL INSERT statement for the Proteins table."""
        sql = f"""INSERT INTO ProteinData.Proteins (Uniprot_Id, Ensemble_Id, ProteinName, Sequence, Length)
                  VALUES ('{self.uniprot_id}', '{self.ensemble_id}', '{self.protein_name}', '{self.sequence}', {self.length});"""
        return sql

    def __str__(self):
        return f"Protein(Id: {self.id}, Uniprot_Id: {self.uniprot_id}, Ensemble_Id: {self.ensemble_id}, ProteinName: {self.protein_name}, Sequence Length: {self.length}, InsertDate: {self.insert_date})"


class ProteinPathogenicity:
    def __init__(self, protein, reference, position, alternate, pathogenicity_score, pathogenic):
        self.protein = protein  # Linking to the Protein object
        self.id = None
        self.reference = reference
        self.position = position
        self.alternate = alternate
        self.variant = f"{reference}{position}{alternate}"
        self.pathogenicity_score = pathogenicity_score
        self.pathogenic = pathogenic
        self.pathogenic_insert_date = None

    def generate_pathogenicity_insert(self):
        """Generate SQL INSERT statement for the Pathogenicity table."""
        sql = f"""INSERT INTO ProteinData.Pathogenicity (ProteinId, Reference, Position, Alternate, PathogenicityScore, Pathogenic)
                  VALUES ({self.protein.id}, '{self.reference}', {self.position}, '{self.alternate}', {self.pathogenicity_score}, '{self.pathogenic}');"""
        return sql

    def __str__(self):
        return f"Pathogenicity(Id: {self.id}, Protein({self.protein.uniprot_id}), Variant: {self.variant}, Score: {self.pathogenicity_score}, Pathogenic: {self.pathogenic}, InsertDate: {self.pathogenic_insert_date})"