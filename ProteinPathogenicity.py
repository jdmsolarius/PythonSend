class Protein:
    def __init__(self, uniprot_id, protein_name=None, sequence=None, length=None, insert_date=None):
        self.id = None
        self.uniprot_id = uniprot_id
        self.protein_name = protein_name
        self.sequence = sequence
        self.length = length
        self.insert_date = insert_date

    def __str__(self):
        return f"Protein(Id: {self.id}, Uniprot_Id: {self.uniprot_id}, ProteinName: {self.protein_name}, Sequence Length: {self.length}, InsertDate: {self.insert_date})"


class ProteinPathogenicity(Protein):
    def __init__(self, uniprot_id, reference, position, alternate, pathogenicity_score, pathogenic,
                 protein_name=None, sequence=None, length=None, insert_date=None):
        super().__init__(uniprot_id, protein_name, sequence, length, insert_date)
        self.pathogenicity_id = None  # This will be set when the pathogenicity data is stored in the database
        self.reference = reference
        self.position = position
        self.alternate = alternate
        self.variant = f"{reference}{position}{alternate}"  # This is derived from the other attributes
        self.pathogenicity_score = pathogenicity_score
        self.pathogenic = pathogenic
        self.pathogenic_insert_date = None  # This will be set when the pathogenicity data is stored in the database

    def generate_pathogenicity_insert(self):
        """Generate SQL INSERT statement for the Pathogenicity table."""
        # Assuming protein_pathogenicity.id is the ID in the Proteins table after it's inserted
        sql = f"""INSERT INTO ProteinData.Pathogenicity (ProteinId, Reference, Position, Alternate, PathogenicityScore, Pathogenic)
                  VALUES ({self.id}, '{self.reference}', {self.position}, '{self.alternate}', {self.pathogenicity_score}, '{self.pathogenic}');"""
        return sql

    def __str__(self):
        base_str = super().__str__()
        return f"{base_str}, Pathogenicity(Id: {self.pathogenicity_id}, Variant: {self.variant}, Score: {self.pathogenicity_score}, Pathogenic: {self.pathogenic}, InsertDate: {self.pathogenic_insert_date})"