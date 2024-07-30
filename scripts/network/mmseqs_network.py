from scripts.cluster.mmseqs_cluster import mmseqs_makedb
from scripts.mmseqs_utils import mmseqs_search
from scripts.utils import running_message
from Bio import SeqIO # type: ignore
from tqdm import tqdm
import threading
import logging
import os

def break_fasta(fasta, outdir, threads):
    total_size = os.path.getsize(fasta)
    os.makedirs(outdir, exist_ok=True)

    files_list= [f"{outdir}/part_{i}.fasta" for i in range(threads)]
    verify_run = True

    for file in files_list:
        if os.path.exists(file):
            pass
        else:
            verify_run = False
            break
    
    if verify_run:
        logging.info("prior split files found, using existing files")
        return files_list
    
    broken_files = [open(f"{outdir}/part_{i}.fasta", 'w') for i in range(threads)]

    with open(fasta) as f, tqdm(total = total_size,
                                  desc= "Reading Fasta",
                                  unit = 'B',
                                  unit_scale = True,
                                  unit_divisor = 1024) as pbar:
        records_written = 0
        for record in SeqIO.parse(f, 'fasta'):
            sublist_index = records_written % threads
            SeqIO.write(record, broken_files[sublist_index], 'fasta')
            records_written += 1
            pbar.update(len(record.format('fasta')))

    for broken_file in broken_files:
        broken_file.close()
    
    return [f"{outdir}/part_{i}.fasta" for i in range(threads)]
            

@running_message
def mmseqs_network(query, outdir, reference_db, threads):
    ref_db_path = f"{outdir}/reference_db"
    query_db_path = f"{outdir}/query_db"
    tmp_path = f"/media/biolinux/nvme0n1p1/sulman_mmseqs_tem"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(ref_db_path, exist_ok=True)
    os.makedirs(query_db_path, exist_ok=True)

    
    ref_files_list = break_fasta(reference_db, ref_db_path, threads)
    
    ref_db_path_list=[]
    for file in ref_files_list:
        ref_path = f"{'.'.join(file.split('.')[0:-1])}"
        os.makedirs(ref_path, exist_ok=True)
        ref_db_path_list.append(ref_path)
    
    thread_list = range(threads)
    threading_process = []

    for i in thread_list:
        x = threading.Thread(target=mmseqs_makedb, args=(ref_files_list[i], ref_db_path_list[i]))
        threading_process.append(x)
        x.start()
    
    for index, threading_process in enumerate(threading_process):
        threading_process.join()
        print(f"Thread {index} finished.")
        
    mmseqs_makedb(query, query_db_path)

    for ref_db in tqdm(ref_db_path_list, desc="searching homologs in reference database", unit="Part"):
        tmp_path = f"{tmp_path}/{ref_db}"
        out_path = f"{outdir}/{ref_db}"
        os.makedirs(out_path, exist_ok=True)
        os.makedirs(tmp_path, exist_ok=True)
        
        mmseqs_search(query_db_path, ref_db, out_path, tmp_path)