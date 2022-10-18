import pdb2sql
import os

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def write_fastafile(cutoff=5.5,pdb_file=None,chain1='C', chain2='D'):
    sql = pdb2sql.interface(pdb_file)
    ctc_res = sql.get_contact_residues(cutoff=cutoff,chain1=chain1, chain2=chain2)
    sql._close()
    ctc_A = ctc_res[chain1]
    ctc_B = ctc_res[chain2]
    sequence = [i[-1] for i in ctc_A] + [i[-1] for i in ctc_B]
    sequence = [d[i] for i in sequence]
    fasta = ''.join(sequence)
    return fasta

path = '/Users/xiaotongxu/Downloads/deep_ranktests'
pdb_files = os.listdir(path)
file_paths = []
for file in pdb_files:
    if file.endswith('.pdb'):
        file_path = os.path.join(path, file)
        file_paths.append(file_path)
sequence = write_fastafile(cutoff=5.5,pdb_file=file_paths[0],chain1='C', chain2='D')

fastas = []
for i in file_paths:
    sequence = write_fastafile(cutoff=5.5,pdb_file=i,chain1='C', chain2='D')
    fastas.append(sequence)

s = ''
for i in range(len(fastas)):
    file_name = file_paths[i].split('/')[-1].split('.')[0]
    lineone = f"\n>{file_name}\n"
    s  += lineone
    s += fastas[i]

with open('/Users/xiaotongxu/Downloads/deep_ranktests/deep_rank_test.fasta', 'w') as f:
    f.write(s)
print(s)

'''
python /Users/xiaotongxu/tests/protein_language_model/esm/scripts/extract.py esm2_t33_650M_UR50D /Users/xiaotongxu/Downloads/deep_ranktests/deep_rank_test.fasta /Users/xiaotongxu/Downloads/deep_ranktests/1AK4_embeddings --repr_layers 0 32 33 --include mean per_tok
'''