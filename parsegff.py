# python3 -m pip install numpy
import numpy as np
from Bio import SeqIO

# Open the GFF3 file


# Open the genome FASTA file
fasta_file = "genome.fa"

# Parse the genome FASTA file

genome = SeqIO.read(fasta_file, "fasta")
gff_table = gff3_parser.parse_gff3(gff_file,parse_attributes=True,verbose=False)
print(genome)

genome_id = gff_table.loc[0,"Seqid"]

assert(genome.id == genome_id) #sanity check
# print(type(gff_table))
for index,row in gff_table.iterrows():
  if index==0:
    continue
  start,end = int(row.Start),int(row.End) # 1-indexed apparently
  print(genome.seq[start-1:end])
# # Parse the GFF3 file
# for record in SeqIO.parse(gff_file, "gff"):
#     print(f"Record ID: {record.id}")
    
#     # Get the corresponding genome sequence
#     genome_seq = genome_dict[record.id].seq
    
#     # Iterate over the features in the record
#     for feature in record.features:
#         print(f"Feature Type: {feature.type}")
#         print(f"Location: {feature.location}")
        
#         # Extract the feature sequence from the genome
#         feature_seq = feature.extract(genome_seq)
#         print(f"Sequence: {feature_seq}")
        
#         print(f"Qualifiers:")
#         for key, value in feature.qualifiers.items():
#             print(f"  {key}: {value}")
#         print("---")