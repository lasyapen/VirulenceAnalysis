
import os
import csv
import sys
from collections import defaultdict
import re


BASE_DIR = "/scratch/lp27068/Olivia1_Data_Reads/virulence_analysis"
MAG_LENGTHS_FILE = os.path.join(BASE_DIR, "mag_lengths.txt")
PARSED_OUTPUT_DIR = os.path.join(BASE_DIR, "parsed_output")

CLASSIFICATION_FILE = os.path.join(BASE_DIR, "annotated_blast_results2.txt")  
OUTPUT_CSV = os.path.join(BASE_DIR, "virulence_abundance.csv")


def parse_mag_lengths(mag_lengths_file):
    mag_lengths = {}
    with open(mag_lengths_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                mag_path, length = line.split('\t')
                mag_path = os.path.abspath(mag_path)
                mag_lengths[mag_path] = int(length)
            except ValueError:
           
    return mag_lengths

def clean_taxonomy(taxonomy_str):
        if '__' in taxonomy_str:
        return taxonomy_str.split('__', 1)[1]
    return taxonomy_str

def parse_classification(classification_file):
        mag_classification = defaultdict(list)
    with open(classification_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                                columns = line.split('\t')
                if len(columns) < 1:
                    continue
                first_field = columns[0
                parts = first_field.split('>')
                if len(parts) != 2:
                    continue
                mag_path_tax = parts[1]

                if ';' not in mag_path_tax:
                    
                    continue
                mag_path, taxonomy_info = mag_path_tax.split(';', 1)

                
                taxonomy_info = taxonomy_info.rstrip('$')

                
                taxonomy_parts = taxonomy_info.split(';')
                raw_organism = taxonomy_parts[-1] if taxonomy_parts else "Unknown"
                organism = clean_taxonomy(raw_organism)

                                mag_path = os.path.abspath(mag_path)

                                mag_classification[mag_path].append(organism)

            except Exception as e:

    final_mag_classification = {}
    for mag_path, organisms in mag_classification.items():
        
        organism_counts = defaultdict(int)
        for org in organisms:
            organism_counts[org] += 1
        sorted_organisms = sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)
        final_organism = sorted_organisms[0][0] if sorted_organisms else "Unknown"
        final_mag_classification[mag_path] = final_organism

        
    for mag_path, organism in list(final_mag_classification.items())[:10]:          print(f"{mag_path} => {organism}", file=sys.stderr)

    return final_mag_classification

def get_all_mag_files(parsed_output_dir):
    mag_files = []
    for location in os.listdir(parsed_output_dir):
        location_path = os.path.join(parsed_output_dir, location)
        if not os.path.isdir(location_path):
            continue
        for distance in os.listdir(location_path):
            distance_path = os.path.join(location_path, distance)
            if not os.path.isdir(distance_path):
                continue
            for sample in os.listdir(distance_path):
                sample_path = os.path.join(distance_path, sample)
                if not os.path.isdir(sample_path):
                    continue
                for mag in os.listdir(sample_path):
                    mag_path_dir = os.path.join(sample_path, mag)
                    if not os.path.isdir(mag_path_dir):
                        continue
                    mag_name = mag
                    blast_file = f"blast_results_{mag_name}.txt"
                    blast_file_path = os.path.join(mag_path_dir, blast_file)
                    if os.path.isfile(blast_file_path):
                                                reconstructed_mag_path = reconstruct_mag_path_from_directory(location, distance, sample, mag_name)
                        mag_files.append( (location, distance, sample, reconstructed_mag_path, blast_file_path) )
                    else:
                        
    return mag_files

def reconstruct_mag_path_from_directory(location, distance, sample, mag_name):
            sample_number = ''.join(filter(str.isdigit, sample))
    sample_name = f"Sample{sample_number}"

        mag_name = re.sub(r'_fa$', '.fa', mag_name)

    mag_name = re.sub(r'_(\d+)\.fa$', r'.\1.fa', mag_name)

    reconstructed_mag_path = os.path.join(
        BASE_DIR,
        f"{location}_AirSamples_Ground_{distance}_{sample_name}_{mag_name}"
    )



    return os.path.abspath(reconstructed_mag_path)

def calculate_hit_length(blast_file_path):
        total_hit_length = 0
    with open(blast_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 8:
                
                continue
            try:
                col7 = int(cols[6])
                col8 = int(cols[7])
                hit_length = col8 - col7
                total_hit_length += hit_length
            except ValueError:
                
    return total_hit_length

def main():
        if not os.path.isfile(MAG_LENGTHS_FILE):
       
        sys.exit(1)

    mag_lengths = parse_mag_lengths(MAG_LENGTHS_FILE)
    if not mag_lengths:
              sys.exit(1)

    
    if not os.path.isfile(CLASSIFICATION_FILE):
        
        sys.exit(1)

    mag_classifications = parse_classification(CLASSIFICATION_FILE)

        if not os.path.isdir(PARSED_OUTPUT_DIR):
        
        sys.exit(1)

    mag_files = get_all_mag_files(PARSED_OUTPUT_DIR)
    if not mag_files:
                sys.exit(1)

    sample_to_mags = defaultdict(list)
    for location, distance, sample, reconstructed_mag_path, blast_file in mag_files:
        sample_key = (location, distance, sample)
        sample_to_mags[sample_key].append( (reconstructed_mag_path, blast_file) )

        with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
                csvwriter.writerow(["Organism", "Location", "Distance_Category", "Sample", "MAG_Name", "Abundance"])

        for sample_key, mags in sample_to_mags.items():
            location, distance, sample = sample_key
            total_sample_length = 0
            mag_lengths_in_sample = {}
            mag_paths_in_sample = {}

            for mag_path, blast_file in mags:
                if mag_path in mag_lengths:
                    mag_length = mag_lengths[mag_path]
                    mag_lengths_in_sample[mag_path] = mag_length
                    mag_paths_in_sample[mag_path] = mag_path
                    total_sample_length += mag_length
                else:
                    print(f"Length for MAG not found: {mag_path}", file=sys.stderr)

            if total_sample_length == 0:
                print(f"Total length for Sample {sample} is zero. Skipping.", file=sys.stderr)
                continue

            for mag_path, blast_file in mags:
                if mag_path not in mag_lengths_in_sample:
                    print(f"Skipping MAG {mag_path} in Sample {sample} due to missing length.", file=sys.stderr)
                    continue
                mag_length = mag_lengths_in_sample[mag_path]
                hit_length_sum = calculate_hit_length(blast_file)
                abundance = hit_length_sum / total_sample_length
                                organism = mag_classifications.get(mag_path, "Unknown")
                original_mag_name = os.path.basename(mag_path)
                csvwriter.writerow([organism, location, distance, sample, original_mag_name, f"{abundance:.6f}"])


if __name__ == "__main__":
    main()
