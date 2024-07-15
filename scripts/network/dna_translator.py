from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from tqdm import tqdm
from scripts.utils import running_message


def translate_dna(dna_sequence):
    frames = []
    seq_len = len(dna_sequence)
    
    # Pad the sequence with 'N' if its length is not a multiple of three
    if seq_len % 3 != 0:
        dna_sequence += 'N' * (3 - seq_len % 3)
        seq_len = len(dna_sequence)
    
    for i in range(3):
        protein = str(Seq(dna_sequence[i:]).translate())
        frames.append(protein[:len(protein) // 3 * 3])

    reverse_complement = str(Seq(dna_sequence).reverse_complement())
    for i in range(3):
        protein = str(Seq(reverse_complement[i:]).translate())
        frames.append(protein[:len(protein) // 3 * 3])

    return frames

def process_record(record):
    dna_sequence = str(record.seq)
    protein_frames = translate_dna(dna_sequence)
    return [f'>{record.id}_frame{i+1}\n{protein}' for i, protein in enumerate(protein_frames)]

@running_message
def translate(input_file, output_file, threads):
    chunk_size = 1000  # Process records in chunks of 1000
    
    with open(output_file, 'w') as output_handle:
        with Pool(processes=threads) as pool:
            # Count total records for progress bar
            total_records = sum(1 for _ in SeqIO.parse(input_file, 'fasta'))
            progress_bar = tqdm(total=total_records, unit='records')

            def process_chunk(records):
                results = pool.map(process_record, records)
                for chunk_results in results:
                    output_handle.write('\n'.join(chunk_results) + '\n')
                progress_bar.update(len(records))

            records = []
            for record in SeqIO.parse(input_file, 'fasta'):
                records.append(record)
                if len(records) == chunk_size:
                    process_chunk(records)
                    records = []

            # Process any remaining records
            if records:
                process_chunk(records)

            progress_bar.close()