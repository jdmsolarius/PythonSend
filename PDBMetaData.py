from symbol import argument
from typing import Self
from DiscoveryMethod import DiscoveryMethod;

class PdbMetaData(DiscoveryMethod):
    def __init__(self, method, pdb_Id, resolutionA, uniprot_id, chains, fileName):
        self.method = DiscoveryMethod(method)
        self.pdb_id = int(pdb_Id)
        self.chains = str(chains)
        self.resolutionA = float(resolutionA.replace("A", "").trim())
        self.uniprot_id = int(uniprot_id)
        self.file_name = "default.pdb"
        
        
        
    


    
    