from multiprocessing import Pool
import subprocess
import os

sub_dir = "test_stricter_human_genome_with_1_nt_ext"
pairagon = True
output_file = "pairagon_seeded_actual.gtf" if pairagon else "gmap.gff"
def run_pairagon(transcript_id: str):
    print("running", transcript_id)
    subprocess.run(
        f"Pairagon.pl -o ./{sub_dir}/{transcript_id}/ -t ./{sub_dir}/{transcript_id}/genome.fa -q ./{sub_dir}/{transcript_id}/cdna.fa -gtf -p {output_file.split('.')[0]} -alignment -a GMap -cross", # prefix is pairagon and also output alignment in case the gff is cancer
        check=True,
        shell=True,
        stdout=subprocess.DEVNULL,
    )
    print("done!", transcript_id)

def run_gmap(transcript_id:str):
    print("running", transcript_id)
    subprocess.run(f"gmap ./{sub_dir}/{transcript_id}/cdna.fa -g ./{sub_dir}/{transcript_id}/genome.fa -E genomic -f 3 > {sub_dir}/{transcript_id}/{output_file}",check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    print("done!", transcript_id)
if __name__ == "__main__":
    ids = [id for x in os.walk(f"./{sub_dir}") if (id := x[0].split("/")[-1]) != sub_dir and output_file not in x[2]] # we haven't processed it yet (pairagon_seeded.gtf is Pairagon's output file)
    print(len([x for x in os.walk(f"./{sub_dir}")])-len(ids)-1,"alignments already done. Skipping...")
    with Pool(50) as p:
       p.map(run_pairagon if pairagon else run_gmap, ids)
