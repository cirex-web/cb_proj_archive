# Generates the test data

import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import subprocess

def reverse_complement(dna:str):
  ans = ""
  nucleotide_map = {'A':'T','C':'G','G':'C','T':'A'}
  for c in dna[-1::-1]:
    ans+=nucleotide_map[c]
  return ans
    
    

def make_path(path):
  if not os.path.exists(path):
    os.makedirs(path)

annotation_file_name = "gencode.v45.basic.annotation.gff3"
transcript_file_name = "gencode.v45.transcripts.fa"
human_genome_file_name = "GRCh38.primary_assembly.genome (1).fa"

# db = gffutils.create_db(annotation_file_name, dbfn='mini.db', force=True, keep_order=False,
# merge_strategy='merge', sort_attribute_values=False)

db = gffutils.FeatureDB('full.db', keep_order=True)
print("Loaded in db!")
transcript_sequences = SeqIO.parse(open(transcript_file_name),'fasta')
transcript_id_to_sequence = {seq_record.id.split("|")[0]:seq_record.seq for seq_record in SeqIO.parse(transcript_file_name,"fasta") }
chromosome_id_to_sequence = {seq_record.id.split(" ")[0]:seq_record.seq for seq_record in SeqIO.parse(human_genome_file_name,"fasta") }
print("Loaded in transcript sequences!")

for gene in db.features_of_type('gene',order_by='start'):
  if gene['level'][0] != "1":
    continue
  
  for transcript in db.children(gene,featuretype='transcript',order_by='start'):
    transcript_id = transcript.id
    if transcript.start is None or transcript.end is None:
      continue
    if transcript_id in transcript_id_to_sequence:
      
      fakeCDNA = "" # we'll just make this as reference to compare against actual transcript... just in case (but we'll still use the actual transcript)
      
      
      reference_cDNA = transcript_id_to_sequence[transcript_id]
      sequence = chromosome_id_to_sequence[transcript.seqid][transcript.start-1:transcript.end]
      exon_count = 0
      exon_intervals = []
      assert(transcript['level'][0]=="1")
      
      for exon in db.children(transcript,featuretype='exon',order_by='start'):
        assert(exon['level'][0]=="1")
        
        assert(exon.start is not None and exon.end is not None)
        startI = exon.start-transcript.start
        endI = exon.end-transcript.start
        if not set(sequence[startI:endI+1]) <= set("ACTG"):
          print("VIOLATION",transcript_id,sequence[startI:endI+1])
          # no N's or lowercase letters allowed (basically we don't want soft-masking/uncertain data)
        fakeCDNA+=sequence[startI:endI+1]
        exon_count+=1
        exon_intervals.append(f"{startI+1}/{endI+1}") # return it back to 1-indexing
      if exon_count <= 1:
        assert(exon_count!= 0)
        continue # really boring transcript
      assert(fakeCDNA==reference_cDNA or fakeCDNA==reverse_complement(reference_cDNA)) 
      
      exon_intervals.append("+" if fakeCDNA==reference_cDNA else "-") # directionality (template vs. non-template strand)
      
        
      # make_path(f"tests/{transcript_id}")
      # SeqIO.write(SeqRecord(sequence,transcript_id,description=gene.id),tfn:=f"tests/{transcript_id}/genome.fa","fasta")
      # SeqIO.write(SeqRecord(reference_cDNA,f"cDNA_{transcript_id}",description=gene.id),cdna_fn:=f"tests/{transcript_id}/cdna.fa","fasta")
      # open(f"tests/{transcript_id}/ans.txt", "w").write("\n".join(exon_intervals))
      # subprocess.run(f"/Users/work/Documents/SD_Projects/gmap-2024/bin/gmap {cdna_fn} -g {tfn} -E genomic -f 3 > tests/{transcript_id}/result.gmap",check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    else:
      print(f"Skipping {transcript_id} cuz not found")
      