#python split_fasta.py input.fasta 5000 output_

import argparse

def split_fasta(input_file, batch_size, output_prefix):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    batch_count = 1
    current_batch = []
    
    for line in lines:
        if line.startswith(">"):
            if len(current_batch) >= batch_size:
                write_output(output_prefix + f'{batch_count:02d}.fasta', current_batch)
                current_batch = []
                batch_count += 1
        current_batch.append(line)

    if current_batch:
        write_output(output_prefix + f'{batch_count:02d}.fasta', current_batch)

    total_reads = len(lines) // 2
    last_batch_reads = len(current_batch) // 2
    remaining_reads = total_reads - (batch_count - 1) * (batch_size // 2) - last_batch_reads

    with open(output_prefix + 'summary.txt', 'w') as summary_file:
        summary_file.write(f"Total Reads: {total_reads}\n")
        summary_file.write(f"Number of Batches: {batch_count - 1}\n")
        summary_file.write(f"Last Batch Reads: {last_batch_reads}\n")
        summary_file.write(f"Remaining Reads: {remaining_reads}\n")

def write_output(output_file, data):
    with open(output_file, 'w') as f:
        f.writelines(data)

def main():
    parser = argparse.ArgumentParser(description='Split Fasta file into batches.')
    parser.add_argument('input_file', help='Path to the input Fasta file')
    parser.add_argument('batch_size', type=int, help='Number of reads per batch')
    parser.add_argument('output_prefix', help='Prefix for the output files')
    args = parser.parse_args()

    split_fasta(args.input_file, args.batch_size, args.output_prefix)

if __name__ == '__main__':
    main()
