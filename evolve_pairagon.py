
import os
from pickle import FALSE
from tokenize import group
import numpy as np
import subprocess
from multiprocessing import Pool
import time
import shutil


def sort_exons_and_check_invariant(exons):
  assert(len(exons)>0) # you never know lol, esp cuz pairagon is so scuffed
  exons.sort()
  for i in range(len(exons)):
    assert (exons[i][0]<=exons[i][1])
    assert(i==len(exons)-1 or exons[i][1]<exons[i+1][0])
    

def create_offspring(cur_gen):
  for i in range(5):
    with open(f"pairagon_generations/{cur_gen}/{i}.zhmm") as current_file:
      with open(f"pairagon_generations/{cur_gen}/{i+5}.zhmm","w") as new_file:
        sum_map = {}
        freq_map = {}
        in_state_transitions = False
        probability_lines = []
        for line in current_file:
          if line=="<STATE_TRANSITIONS>\n":
            in_state_transitions = True
            
          if line == "<STATE_DURATIONS>\n":
            in_state_transitions = False
            for p_line in probability_lines:
              if sum_map[p_line[0]]<=.000000000001:
                p_line[2] = 1/freq_map[p_line[0]] #just balance it out - lol (cuz everything is 0)
              else:
                p_line[2]/=sum_map[p_line[0]]
              p_line[2] = str(p_line[2])
              # print(p_line)
              new_file.write("\t".join(p_line)+"\n")
          
          if in_state_transitions and len(line.split("\t")) == 3:
              a,b,prob = line.split("\t")
              prob = float(prob)
              prob = max(0,(prob+np.random.normal(0,.05))) # keep it in there lol
              sum_map[a] = sum_map[a]+prob if a in sum_map else prob
              freq_map[a] = freq_map[a]+1 if a in freq_map else 1
              probability_lines.append([a,b,prob])
          else:
            new_file.write(line)

def run_pairagon(path: str,child_id:int,generation):
    print("running", path,generation,child_id)
    subprocess.run(
        f"Pairagon.pl -o ./{path}/ -t ./{path}/genome.fa -q ./{path}/cdna.fa -gtf -p '{child_id}'_ -alignment -a GMap -param pairagon_generations/{generation}/{child_id}.zhmm", # prefix is pairagon and also output alignment in case the gff is cancer
        check=True,
        shell=True,
        stdout=subprocess.DEVNULL,
    )
    time.sleep(1)
    print("done!", path,child_id)
  
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
  specificity = true_positives/(true_positives+false_positives)
  # print(sensitivity,accuracy)
  return sensitivity,specificity

    
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
  specificity = true_positives/(true_positives+false_positives)
  # print(sensitivity,accuracy)
  return sensitivity,specificity

def test_and_propogate_offspring(cur_gen,top_5_scores, with_proportionate_selection=False):
  print(f"Beggining Generation {cur_gen}")
  children_ids = range(5,10) # the only new ones
  test_paths = [test_path[0] for test_path in os.walk("mini_test") if test_path[0] != "mini_test"]
  # test_paths = test_paths[0:1]
  with Pool(32) as p:
    p.starmap(run_pairagon,[(path,child,cur_gen) for path in test_paths for child in children_ids])
  scores = [[top_5_scores[i],i] for i in range(len(top_5_scores))]
  
  # score what we just ran
  for child in children_ids:
    total_score = 0
    for path in test_paths:
      true_exons,reversed = extract_exons_from_ans(f"{path}/ans.txt")
      
      pairagon_exons,pairagon_reversed = extract_exons_from_pairagon_file(f"{path}/{child}_.gtf")
      # if reversed != pairagon_reversed:
      #   continue # 0 score
      se,sp = score_on_nucleotide_lvl(true_exons,reversed,pairagon_exons,pairagon_reversed)
      total_score+=se+sp
    scores.append([total_score,child])
  unsorted_scores = [str(score[0]) for score in scores]
  scores.sort(reverse=True,key=lambda a:a[0]) #don't sort based on child_id. only score

  # create directory for next generation
  if not os.path.exists(f"./pairagon_generations/{cur_gen+1}"): # !!! can just use os.makedirs with exists_ok = True
    os.mkdir(f"./pairagon_generations/{cur_gen+1}")

  if with_proportionate_selection:

    # determine total "fitness" of the group
    aggregate_group_fitness = sum(score - 115 if score > 115 else 0 for score, _ in scores) # subtracting 115 from every score to heighten signal. Kind of hackey, but if you don't do this probabilities will all be essentially uniform.

    # assign probabilities for each individual to be selected
    group_probabilities = [(score - 115) / aggregate_group_fitness if score > 115 else 0 for score, _ in scores] # set individual probability of reproduction equal to relative weight of group fitness
    assert(abs(1 - sum(group_probabilities)) <= 0.001) # just make sure it's close enough

    # probablistically select individuals
    children_ids_bank = range(0, 10)
    selected_children_ids = np.random.choice(children_ids_bank, size=5, replace=False, p=group_probabilities)
    new_children = [scores[i] for i in selected_children_ids]


  else: # -- truncated child selection --
    latest_best = 4 # we'll at least take the first 5
    for i in range(len(scores)): # this just finds the most recent best score
      if scores[i][0]!=scores[0][0]: 
        break
      latest_best = max(latest_best,i)

    # get the new children we will be bringing to the next generation
    # this is just the most recent top5 best scores
    new_children = [scores[i] for i in range(latest_best-5+1,latest_best+1)]
      

  # new_child_id = 0
  # for i in range(latest_best-5+1,latest_best+1):
  #   child_id = scores[i][1]
  #   print("yanking",child_id)
  #   shutil.copy(f"pairagon_generations/{cur_gen}/{child_id}.zhmm",f"pairagon_generations/{cur_gen+1}/{new_child_id}.zhmm")
  #   new_child_id+=1
  #   new_5.append(scores[i][0])

  new_5 = [] # renamed this to "new_5" because that's what makes sense in conjunction with proportional child selection. Really it shouldn't be hardcoded to 5 but its whatever
  # idx will be what this child will be in the new generation. new generation contains 5 from the old generation
  for idx, child in enumerate(new_children):
    child_id = scores[idx][1]
    print("yanking", child_id)
    shutil.copy(f"pairagon_generations/{cur_gen}/{child_id}.zhmm",f"pairagon_generations/{cur_gen+1}/{idx}.zhmm")
    new_5.append(scores[idx][0])

  return unsorted_scores,new_5

gen = None
best_5 = None

# when continuing from existing scoring
with open("scores.tsv") as score_file:
  last_line = ""
  for line in score_file:
    if len(line)<=1:
      continue # sucky line
    last_line = line
  generation,*scores = last_line.strip().split("\t")
  gen = int(generation) + 1 #next gen baby!
  scores_numerical = [float(score) for score in scores]
  scores_numerical.sort(reverse=True)
  best_5 = scores_numerical[:5]

# gen = 26
# best_5 = [115.77101904002731, 115.77101904002731, 115.77101904002731, 115.75221670715541, 115.75221670715541] # original run (3 pairagon.zhmm and 2 pairagonx.zhmm)
# best_5 = [115.77451816107964, 115.77451816107964, 115.77451816107964, 115.77101904002731, 115.77101904002731]
# best_5 = [115.78322179930521,115.78322179930521,115.78322179930521,115.78322179930521,115.78322179930521]

print(f"Starting Generation: {gen} \n Best 5: {best_5}")

with open("scores.tsv","a") as score_file:
  while True:
    create_offspring(gen) # fill in last 5
    scores,new_best_5 = test_and_propogate_offspring(gen,best_5,with_proportionate_selection=True) # score and propogate best 5 to next gen (cache best-5 scores in variable)
    score_txt = "\t".join(scores)
    score_file.write(f"{gen}\t{score_txt}\n")
    score_file.flush()
    print(scores,new_best_5) # scores are unsorted so they directly correspond with 0,1,2...new_best_5 corresponds with 0,1,2,3,4 of next gen, respectively (cuz we just copied it over)
    best_5 = new_best_5
    gen+=1
# nohup python3 -u evolve_pairagon.py &
# Patch1 (right before gen 26): Change min from .1 to .000001
# patch2 (right before gen 35): Genetic drift enabled (if generation has more than 5 optimal scoring organisms, kill off the oldest)
# also if sum_map[p_line[0]]<=.000000000001:
      #   p_line[2] = 1/freq_map[p_line[0]] #just balance it out - lol (cuz everything is 0)
      # else:
      #   p_line[2]/=sum_map[p_line[0]]
      
# gens 33 and below are uncorrupted

# patch3 enabled for probabilistic child selection as opposed to randomized selection