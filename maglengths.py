cat *.fa > concatenated_query.fasta
fasta_file_path = "concatenated_query.fasta" 
output_file_path = "mag_lengths.txt"

mag_lengths = {}
current_mag = None

with open(fasta_file_path, "r") as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
          if "k141_" not in line:
            current_mag = line[1:]
            if current_mag not in mag_lengths:
              mag_lengths[current_mag] = 0
         else:
            mag_lengths[current_mag] += len(line)

with open(output_file_path, "w") as output_file:
    for mag, length in mag_lengths.items():
        output_file.write(f"{mag}\t{length}\n")
