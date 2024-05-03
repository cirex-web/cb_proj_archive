#! /home/work/cb-proj/path/to/venv/bin/python3
import gtfparse
import polars as pl
import matplotlib.pyplot as plt
import random
import os
from Bio import SeqIO


###################
def load_fasta(file_path = 'GRCm39.genome.fa'):
    # This will load a FASTA file and return a dictionary of sequences.
    # The keys are the sequence identifiers, and the values are the sequences.
    sequences = {}
    for seq_record in SeqIO.parse(file_path, "fasta"):
        sequences[seq_record.id] = seq_record

    print("Successfully Loaded Fasta File")
    return sequences
############################################################

def load_cDNA(file_path):
    sequences = {}
    for seq_record in SeqIO.parse(file_path, "fasta"):
        sequences[seq_record.id.split("|")[0]]  = seq_record

    print("Successfully Loaded Transcript Fasta File")
    return sequences

###################
def extract_exon_sequences(sequences, exons):
    # Exons should be a list of tuples (start, end) for the exon positions.
    # Sequences is a dictionary from load_fasta.
    extracted_sequences = {}
    for gene_id, sequence in sequences.items():
        gene_exons = []
        for start, end in exons:
            # Biopython uses 0-based indexing, subtract 1 from start positions
            exon_seq = sequence[start-1:end]
            gene_exons.append(exon_seq)
        extracted_sequences[gene_id] = gene_exons
    return extracted_sequences
############################################################

###################
def analyzeGenome(transcript_group: dict):
    # Initialize a dictionary to store the count of exons for each gene
    exon_counts = {}

    # Loop through each grouped DataFrame to count the exons
    for transcript_name, group_df in transcript_group.items():
        # Count the number of rows in the DataFrame for each gene, which represents the number of exons
        exon_count = group_df.shape[0]
        exon_counts[transcript_name] = exon_count  # Correctly update exon_counts instead of transcript_group
        print(f"{transcript_name} has {exon_count} exons")


    # Count the frequency of each exon count
    exon_frequency = {}
    for exon in exon_counts.values():
        if exon in exon_frequency:
            exon_frequency[exon] += 1
        else:
            exon_frequency[exon] = 1

    # Data for plotting
    exon_numbers = list(exon_frequency.keys())  # Number of exons (x-axis)
    frequency = list(exon_frequency.values())  # Frequency of each exon count (y-axis)

    # Create bar chart
    plt.bar(exon_numbers, frequency, color='blue')

    # Adding title and labels
    plt.title('Frequency of Exon Counts per Transcript')
    plt.xlabel('Number of Exons')
    plt.ylabel('Frequency')

    # Show plot
    plt.yscale('log')
    plt.xlim(0, 1000)  # Adjust this based on your data

    plt.savefig("ExonCountPerTranscript.png")
############################################################
    
###################
def loadingGenome(): # returns transcript_group
    # Load the GTF file
    gtf_file = 'gencode.vM34.chr_patch_hapl_scaff.basic.annotation.gtf'
    gtf_data = gtfparse.read_gtf(gtf_file)


    # filter to get only exons and only lab verified
    exons = gtf_data.filter(gtf_data['feature'] == 'exon')
    # exons = gtf_data.filter(gtf_data['level'] == 1)
    

    exon_details = exons.select(['seqname', 'start', 'end', 'feature', 'transcript_name', 'gene_id', 'transcript_id'])
    exon_details = exon_details.with_columns((exons['end'] - exons['start'] + 1).alias('exon_length'))
    exon_details = exon_details.sort(['transcript_name', 'start'])

    # sort exons by transcript
    transcript_group = {name: group for name, group in exon_details.group_by("transcript_id")}

    print("Successfully Loaded Annotated Genome")
    return transcript_group
############################################################

###################
def mutate_dna(sequence, mutation_rate=0.01):
    """
    Mutate a given percentage of nucleotides in the DNA sequence.
    """
    # Calculate the total number of mutations needed
    n_mutations = int(len(sequence) * mutation_rate)
    
    # Create a set to store mutation indices (to avoid duplicating mutations)
    mutation_indices = set()
    while len(mutation_indices) < n_mutations:
        mutation_indices.add(random.randint(0, len(sequence) - 1))
    
    # Convert sequence to a list for mutability
    sequence_list = list(sequence)
    
    # Possible nucleotides
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Perform mutations
    for index in mutation_indices:
        original_nucleotide = sequence_list[index]
        new_nucleotide = random.choice([n for n in nucleotides if n != original_nucleotide])
        sequence_list[index] = new_nucleotide
    
    return ''.join(sequence_list)
############################################################

###################
def createTranscripts(transcript_group, new_dir, count = 100, mutation_rate = None):

    if mutation_rate:
        ogloc = new_dir + '-' + str(mutation_rate)
    else: 
        ogloc = new_dir
        assert(False)

    os.makedirs(ogloc, exist_ok=True)

    # load the sequences
    sequences = load_fasta()

    transcript_id_to_sequence = load_cDNA(file_path = 'gencode.vM34.transcripts.fa')

    # loop through every gene in transcript_group
    for transcript in list(transcript_group.keys()):
        count -= 1
        if count <= 0: 
            print("Successfully Created Artificial cDNA strands")
            break

        if transcript not in transcript_id_to_sequence:
            continue

        
        # only do transcripts more than 3
        exon_count = transcript_group[transcript].shape[0]
        if exon_count <= 2:
            continue 

        # make a folder to store the mRNA strings for this gene
        loc = os.path.join(ogloc, transcript)  # Use os.path.join for OS-independent path handling
        os.makedirs(loc, exist_ok=True)  # Use exist_ok here as well


        # randomly choose how many exons will be present in this transcript
        number_of_exons = exon_count
        # number_of_exons = random.randint(1, exon_count)

        # randomly choose a subset of size number_of_exons
        exon_subset_indices = random.sample(range(1, exon_count), number_of_exons - 1) # never select the first exon cuz then you can't pad, sort of hacky :(
        exon_subset_indices.sort() # make sure the order is sorted

        # initialize exon_boundaries
        exon_boundaries = [] # list of tuples

        # store this in our folder
        


        # create the genome.fa
        file_path = os.path.join(loc, 'genome.fa')
        chromosome_id = transcript_group[transcript]["seqname"][0]
        sequence = sequences[chromosome_id]
        len_of_sequence = len(sequence)


        if not sequence:
            print("SEQUENCE IS NONE")
            continue
        
        length_of_ref_genome = 0
        with open(file_path, 'w') as file:
            file.write('>' + 'cDNA_' + transcript + '\n')
            # at 2 nt extension
            reference_genome = sequence.seq[ transcript_group[transcript]['start'][0] - 2 : transcript_group[transcript]['end'][-1] + 1 ]
            file.write(str(reference_genome))
            length_of_ref_genome = len(reference_genome)
            assert(transcript_group[transcript]['end'][-1] - (transcript_group[transcript]['start'][0] - 1) <= len(reference_genome))

        file_path = os.path.join(loc, 'ans.txt')
        with open(file_path, 'w') as file:  # Using 'with' ensures the file is closed properly after its suite finishes
            for exon_idx in exon_subset_indices: 
                start, end = transcript_group[transcript]['start'][exon_idx], transcript_group[transcript]['end'][exon_idx]

                assert(type(start) == type(end) == int)
                offset = transcript_group[transcript]['start'][0] ## chr 1000, transcript: 10, transcript_start = 990
                exon_boundaries.append((start, end))
                file.write(str(start - offset + 2) + '/' + str(end - offset + 2) + '\n')
                assert(start - offset <= length_of_ref_genome and end - offset <= length_of_ref_genome)
            file.write('+')

        # create file to know what sequence this transcript comes from (to know what to align it to)
        file_path = os.path.join(loc, 'sequence_ref.txt')
        with open(file_path, 'w') as file:
            file.write(transcript_group[transcript]['seqname'][0])
        
        # load in the transcripts
        file_path = os.path.join(loc, 'cdna.fa')
        with open(file_path, 'w') as file:  # Using 'with' ensures the file is closed properly after its suite finishes
            file.write('>' + 'cDNA_' + transcript + '\n')

            if sequence: 

                # ensuring the data we are adding is actually being added
                current_length = 0
                expected_length = 0


                for start, end in exon_boundaries: 
                    assert(type(start) == type(end) == int)

                    expected_length += end - start + 1

                    if (start >= len(sequence)): 
                        print("starting index is greater than length of sequence")
                        continue
                    
                    if (end >= len(sequence)):
                        print("ending index is greater than the length of sequence")
                        continue
                    
                    # 2nt padding 
                    exon_seq = sequence[start - 1: end]

                    current_length += len(exon_seq) 
                    
                    if mutation_rate:
                        nucleotide_sequence = mutate_dna(str(exon_seq.seq), mutation_rate)
                        file.write(nucleotide_sequence)
                    else:
                        file.write(str(exon_seq.seq))
            else: 
                print("Sequence type is NONE")
                continue
            
            if (expected_length != current_length):
                print("Expected Length: " + str(expected_length), "Current Length: " + str(current_length))

            file.write('\n')

    print("Successfully Created Artificial cDNA strands")
############################################################


###################
def extractSequences():
    print("Attempting to Extract Sequences from Genome")

    # load the sequence
    sequences = load_fasta()

    # make a folder to store the sequences
    loc = 'nucleotideSequences'
    os.makedirs(loc, exist_ok=True)

    # for every sequence, store it as a file
    for sequence in list(sequences.keys()):

        file_path = os.path.join(loc, sequence + '_genome.fa')
        with open(file_path, 'w') as file:  # Using 'with' ensures the file is closed properly after its suite finishes
            file.write('>' + sequence + '\n')
            file.write(str(sequences[sequence]))


    print("Successfully Extracted Sequences from Genome")
############################################################


transcript_group = loadingGenome()


# analyzeGenome(transcript_group)
createTranscripts(transcript_group, 'mutatedRNAs', mutation_rate = 0.01)
createTranscripts(transcript_group, 'mutatedRNAs', mutation_rate = 0.05)
createTranscripts(transcript_group, 'mutatedRNAs', mutation_rate = 0.1)
createTranscripts(transcript_group, 'mutatedRNAs', mutation_rate = 0.25)
# createMutatedTranscripts(transcript_group)
# extractSequences()
    

