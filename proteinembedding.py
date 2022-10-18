import warnings
import os
import torch
from deeprank.features import FeatureClass
from deeprank import config
import pdb2sql
from pprint import pprint

class ProteinEmbedding(FeatureClass):
    def __init__(self, pt_path=None, pdb_file=None,chain1='A', chain2='B'):
        super().__init__("Residue")
        
        self.pt_path= pt_path
        self.pdb_file = pdb_file
        self.chain1 = chain1
        self.chain2 = chain2

        self.feature_data = {}
        self.feature_data_xyz = {}
        res_names = config.AA_codes_pssm_ordered
        self.feature_names = tuple(['embedding_' + n for n in res_names])
        for name in self.feature_names:
                self.feature_data[name] = {}
                self.feature_data_xyz[name] = {}

    def read_pt_data(self,cutoff=5.5):
        sql = pdb2sql.interface(self.pdb_file)

        xyz_info, xyz = self.get_residue_center(sql)
        xyz_dict = {}
        for pos, info in zip(xyz, xyz_info):
            xyz_dict[tuple(info)] = pos

        ctc_res = sql.get_contact_residues(cutoff=cutoff,
                            chain1=self.chain1, chain2=self.chain2)
        sql._close()
        ctc_res = ctc_res[self.chain1] + ctc_res[self.chain2]
        total_len = len(ctc_res)


        if total_len  == 0:
            raise ValueError(
                f" Failed to calculate the features of protein embedding.")
        elif total_len  < 5:  # this is an empirical value
            warnings.warn(
                f"Only {total_len } interface residues found"
                f" with cutoff {cutoff}Å. Be careful with")
            
        #pt_path = os.path.join(pt_path, str(pdb_file).split()[-1])
        #need to change after the embedding generation code
        result = torch.load(self.pt_path)["representations"][33]
        #print(result.size())
        

        self.embed = {}
        for i in range(len(ctc_res)):
            self.embed[ctc_res[i]] = result[i]


        for res in ctc_res:
            #print(res)
            chain = {self.chain1: 0, self.chain2: 1}[res[0]]
            key = tuple([chain] + xyz_dict[res])
            #print(key)
            for name, value in zip(self.feature_names, self.embed):
                self.feature_data[name][res] = self.embed[res]
                self.feature_data_xyz[name][key] = self.embed[res]


########################################################################
#
#   THE MAIN FUNCTION CALLED IN THE INTERNAL FEATURE CALCULATOR
#
########################################################################


def __compute_feature__(pdb_file, featgrp, featgrp_raw, chain1, chain2, out_type='protein_embedding'):
    """Main function called in deeprank for the feature calculations.
    Args:
        pdb_file (list(bytes)): pdb information
        featgrp (str): name of the group where to save xyz-val data
        featgrp_raw (str): name of the group where to save human readable data
        chain1 (str): First chain ID
        chain2 (str): Second chain ID
        out_type (str): which feature to generate, 'protein_embedding'
    """

    if config.PATH_PE_SOURCE is None:
        raise FileExistsError(f"No available PE source, "
                    f"check 'config.PATH_PSSM_SOURCE'")
    else:
        path = config.PATH_PSSM_SOURCE 

    pe = ProteinEmbedding(pt_path, pdb_file, chain1=chain1, chain2=chain2, out_type=out_type)

    # read and get the raw data
    pe.read_pt_data()


    # export in the hdf5 file
    pe.export_dataxyz_hdf5(featgrp)
    pe.export_data_hdf5(featgrp_raw)

'''
########################################################################
#
#   IF WE JUST TEST THE CLASS
#
########################################################################
if __name__ == '__main__':
    pe = ProteinEmbedding(pt_path='Protein_test.pt', pdb_file='7cei.pdb')
    pe.read_pt_data()
    pprint(pe.feature_data)

'''
