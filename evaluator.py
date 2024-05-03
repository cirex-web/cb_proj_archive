#! path/to/venv/bin/python3
# no joke I actually called the path this lol (yeah I'm lazy)
from BCBio import GFF
from Bio import SeqIO
import os
import subprocess
from BCBio.GFF import GFFExaminer



def sort_exons_and_check_invariant(exons):
  assert(len(exons)>0) # you never know lol, esp cuz pairagon is so scuffed
  exons.sort()
  for i in range(len(exons)):
    assert (exons[i][0]<=exons[i][1])
    assert(i==len(exons)-1 or exons[i][1]<exons[i+1][0])
    
"""
Returns exon array, 0-indexed (coords are in reference to the reference genome), and sorted in ascending order (even if reference strand is template strand rather than coding strand)
"""
def extract_exons_from_gmap_file(gmap_file_path):
  exons = [] # 0-indexed coords, both inclusive
  reverse:bool|None = None
  with open(gmap_file_path) as gmap_file:
    recs = []
    for rec in GFF.parse(gmap_file): # wait I really don't know if this is right or not also do we really need the sequence itself is the real question
      recs.append(rec)
    assert(len(recs)==1) # surely ...
    for match in rec.features:
      assert (match.type == "cDNA_match")
      exons.append([int(match.location.start),match.location.end-1])
      reverse = bool(match.location.strand==-1)

  sort_exons_and_check_invariant(exons)
  assert( reverse is not None)
  return exons,reverse

def extract_exons_from_ans(ans_file_path):
  lines = [line for line in open(ans_file_path)]
  is_reversed = lines.pop()=='-'
  exons = []
  for line in lines:
    start,end = line.split("/")
    exons.append([int(start)-1,int(end)-1])
  sort_exons_and_check_invariant(exons)
  return exons, is_reversed

def extract_exons_from_pairagon_file(pairagon_gtf_path):
  is_reversed = None
  exons = []

  with open(pairagon_gtf_path) as gtf_file:
    for line in gtf_file: #ya know what screw BioPython
      row_data = line.split("\t") 
      is_reversed = row_data[6]=='-'
      start,end = int(row_data[3]),int(row_data[4]) # already 0-indexed because pairagon is pairagon
      exons.append([start,end])
    sort_exons_and_check_invariant(exons)
  assert(is_reversed is not None)
  return exons,is_reversed

def score_on_nucleotide_lvl(true_exons:list[tuple[int,int]],true_reversed:bool,found_exons:list[tuple[int,int]],found_reversed:bool):
  if true_reversed != found_reversed:
    return 0 # There's no hope lol
  assert(len(true_exons)>=1 and len(found_exons)>=1)
  startNucleotide = min(true_exons[0][0],found_exons[0][0])
  endNucleotide = max(true_exons[-1][1],found_exons[-1][1])
  
  i = j = 0 # what exon interval we're currently looking at
  
  false_positives = false_negatives = true_positives = 0 # true negatives apparently don't matter
  for pos in range(startNucleotide,endNucleotide+1):
    while(i<len(true_exons)and true_exons[i][1]<pos):
      i+=1
    while j<len(found_exons) and found_exons[j][1]<pos:
      j+=1
    
    is_true_exon = i!=len(true_exons) and pos>=true_exons[i][0]
    is_found_exon = j!=len(found_exons) and pos >= found_exons[j][0]
    # print(pos,is_true_exon,is_found_exon)
    if is_true_exon:
      if is_found_exon:
        true_positives+=1
      else:
        false_negatives+=1
    else:
      if is_found_exon:
        false_positives+=1
  sensitivity = true_positives/(true_positives+false_negatives) 
  accuracy = true_positives/(true_positives+false_positives)
  # print(sensitivity,accuracy)
  if sensitivity == 0 or accuracy == 0:
    return 0 # when you suck so much 
  return 2/(1/sensitivity+1/accuracy) # HM of the two
    
# Assumes that exon lists are sorted and disjoint 
# move tru exon pointer forwards until true end is after current start
# after attempted match, move true exon ptr until true start is > cur end (not needed)
def score_on_exon_lvl(true_exons:list[tuple[int,int]],true_reversed:bool,found_exons:list[tuple[int,int]],found_reversed:bool):
  if true_reversed != found_reversed: # I mean okay if you messed this up you actually kinda suck
    return 0
  true_i = 0
  false_positives = false_negatives = true_positives = 0
  for [start,end] in found_exons:
    while true_i<len(true_exons) and true_exons[true_i][1]<start:
      true_i+=1
      false_negatives+=1 # we just never matched this true exon :(
    if true_i<len(true_exons) and true_exons[true_i] == [start,end]:
      true_positives+=1
      true_i+=1
    else:
      false_positives+=1 # found exon is wrong lol
  false_negatives+=len(true_exons)-true_i # everything that we didn't get to
  sensitivity = true_positives/(true_positives+false_negatives) 
  accuracy = true_positives/(true_positives+false_positives)
  # print(sensitivity,accuracy)
  if sensitivity == 0 or accuracy == 0:
    return 0 # when you suck so much 
  return 2/(1/(sensitivity)+1/(accuracy)) # HM of the two
  
def get_scores():
  print("fasdf")
  transcript_to_score_map = {}
  # read in files 
  for path,_,files in os.walk("./derek/simulatedmRNAs"):
    if path == "./derek/simulatedmRNAs":
      continue
    transcript_id = path.replace("./derek/simulatedmRNAs/","")
    assert ("gmap.gff" in files) # GMAP is op (100% coverage)
    assert("pairagon_seeded.gtf" in files)

    # Just temp code to transfer pairagon alignments to non-corrupt .gtf format
    # subprocess.run(f"alignmentConvert.pl -i {path}/pairagon.pair -o gtf -q {path}/cdna.fa -t {path}/genome.fa > {path}/pairagon.gtf",check=True,shell=True)

    # if transcript_id != "ENST00000397517.6":
    #   continue
    
    # print(transcript_id)
    # seq_dict = SeqIO.to_dict(SeqIO.parse(f"{path}/genome.fa", "fasta")) 

    # orientation 
    gmap_exons, gmap_is_reversed = extract_exons_from_gmap_file(f"{path}/gmap.gff")
    pairagon_exons, pairagon_is_reversed = extract_exons_from_pairagon_file(f"{path}/pairagon_seeded.gtf")
    true_exons,true_is_reversed = extract_exons_from_ans(f"{path}/ans.txt")
    
    
    # print(gmap_is_reversed,gmap_exons)
    # print(pairagon_is_reversed,pairagon_exons)
    # print(true_is_reversed,true_exons)

    # print(score_on_exon_lvl([[1,1],[2,2],[3,9],[10,12]],False,[[1,1],[2,2],[3,3],[4,4],[10,12]],False)) # type: ignore
    gmap_exon_score = score_on_exon_lvl(true_exons,true_is_reversed,gmap_exons,gmap_is_reversed)
    pairagon_exon_score = score_on_exon_lvl(true_exons,true_is_reversed,pairagon_exons,pairagon_is_reversed)

    gmap_nucleotide_score = score_on_nucleotide_lvl(true_exons,true_is_reversed,gmap_exons,gmap_is_reversed)
    pairagon_nucleotide_score = score_on_nucleotide_lvl(true_exons,true_is_reversed,pairagon_exons,pairagon_is_reversed)
    transcript_to_score_map[transcript_id] = {
      "gmap_exon":gmap_exon_score,
      "gmap_nucleotide":gmap_nucleotide_score,
      "pairagon_exon":pairagon_exon_score,
      "pairagon_nucleotide":pairagon_nucleotide_score
    }
    # print(transcript_id,gmap_exon_score,gmap_nucleotide_score,pairagon_exon_score,pairagon_nucleotide_score)
  
  return transcript_to_score_map
    # print(gmap_exons,gmap_is_reverse)
  # convert everything over to exon segments 

T = get_scores()
for key, value in T.items(): 
  print(key)
  print(value)