from multiprocessing import Pool
import subprocess
import os

# referenceGenomePath = "nucleotideSequences"

def run_pairagon(sub_dir: str, transcript_id: str):

    print("running", transcript_id)

    print(f"Pairagon.pl -o ./{sub_dir}/{transcript_id}/ -t ./{sub_dir}/{transcript_id}/genome.fa -q ./{sub_dir}/{transcript_id}/cdna.fa -gtf -p pairagon_seeded -alignment -a GMap")
    subprocess.run(
        f"Pairagon.pl -o ./{sub_dir}/{transcript_id}/ -t ./{sub_dir}/{transcript_id}/genome.fa -q ./{sub_dir}/{transcript_id}/cdna.fa -gtf -p pairagon_seeded -alignment -a GMap", 
        check=True,
        shell=True,
        stdout=subprocess.DEVNULL,
    )
    print("done!", transcript_id)

def run_gmap(sub_dir: str, transcript_id:str):
    print("running", transcript_id)



    # get the reference sequence
    reference_seq = ''
    with open(f'./{sub_dir}/{transcript_id}/sequence_ref.txt', 'r') as file:
        reference_seq = file.read().strip()

    subprocess.run(f"gmap ./{sub_dir}/{transcript_id}/cdna.fa -g ./{sub_dir}/{transcript_id}/genome.fa -E genomic -f 3 > {sub_dir}/{transcript_id}/gmap.gff",check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    print("done!", transcript_id)



pairagon_output_file = "pairagon_seeded.gtf"
gmap_output_file = "gmap.gff"

def test_subdir_pairagon(sub_dir: str):
    ids = [id for x in os.walk(f"./{sub_dir}") if (id := x[0].split("/")[-1]) != sub_dir and pairagon_output_file not in x[2]] # we haven't processed it yet (pairagon_seeded.gtf is Pairagon's output file)
    with Pool(35) as p:
       p.starmap(run_pairagon, [(sub_dir, id) for id in ids])

def test_subdir_gmap(sub_dir: str):
    ids = [id for x in os.walk(f"./{sub_dir}") if (id := x[0].split("/")[-1]) != sub_dir and gmap_output_file not in x[2]] # we haven't processed it yet (pairagon_seeded.gtf is Pairagon's output file)
    with Pool(35) as p:
       p.starmap(run_gmap, [(sub_dir, id) for id in ids])


if __name__ == "__main__":

    # run simulated mRNAs
    # test_subdir_pairagon("simulatedmRNAs")
    # test_subdir_gmap("mutatedRNAs-0.1")
    # test_subdir_gmap("mutatedRNAs-0.01")
    # test_subdir_gmap("mutatedRNAs-0.05")
    # test_subdir_gmap("mutatedRNAs-0.25")

    test_subdir_pairagon("mutatedRNAs-0.1")
    test_subdir_pairagon("mutatedRNAs-0.01")
    test_subdir_pairagon("mutatedRNAs-0.05")
    test_subdir_pairagon("mutatedRNAs-0.25")

    # # # run mutated simulated mRNAs
    # test_subdir_pairagon("mutatedSimulatedmRNAs")
    # test_subdir_gmap("mutatedSimulatedmRNAs")




    
