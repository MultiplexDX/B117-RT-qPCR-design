import pandas as pd
from Bio import SeqIO

fasta_sequences = SeqIO.parse(open('sars2_8thApril2021/msa_0406/msa_0406.fasta'),'fasta')
meta_data = pd.read_csv('sars2_8thApril2021/metadata.tsv', delimiter="\t")

strains = tuple(meta_data['Virus name'])
epi = tuple(meta_data['Accession ID'])
host = tuple(meta_data['Host'])
col_date = tuple(meta_data['Collection date'])
country = tuple(meta_data['Location'])
pango_lin = tuple(meta_data['Pango lineage'])
clade = tuple(meta_data['Clade']) 

for fasta in fasta_sequences:
    faEpi = fasta.description.split("|")
    if faEpi[1] in epi:        
        myIndex = epi.index(faEpi[1])
# array starts with 0. position
# position of 6 bp deletion (HV) in the spike + a bit around: 24212. - 24234. bp
# position of 3 bp deletion (Y) in the spike + a bit around: 24482. - 24504. bp
        print(strains[myIndex], epi[myIndex], host[myIndex], col_date[myIndex], country[myIndex], pango_lin[myIndex], clade[myIndex], str(fasta.seq)[24212:24234], str(fasta.seq)[24482:24504])
