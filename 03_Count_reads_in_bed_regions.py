import argparse
import pysam

# Parse command line arguments
parser = argparse.ArgumentParser(description="Count reads in each region of a BED file from a SAM file.")
parser.add_argument("samfile", help="Input SAM file")
parser.add_argument("bedfile", help="Input BED file")
parser.add_argument("output_samfile", help="Output SAM file")
parser.add_argument("output_counts_file", help="Output file for read counts")
parser.add_argument("-q", "--mapq", type=int, default=0, help="Minimum mapping quality threshold (default: 0)")
args = parser.parse_args()

# Load the BED regions into memory
bed_regions = []
with open(args.bedfile) as bedfile:
    for line in bedfile:
        fields = line.strip().split()
        region = (fields[0], int(fields[1]), int(fields[2]))
        bed_regions.append(region)

# Initialize a dictionary to store the read counts for each region
read_counts = {region: 0 for region in bed_regions}

# Open the input and output SAM files
input_samfile = pysam.AlignmentFile(args.samfile, "r")
output_samfile = pysam.AlignmentFile(args.output_samfile, "wh", template=input_samfile)

# Loop through each read in the input SAM file
for read in input_samfile:
    # Check if the read meets the minimum mapping quality threshold
    if read.mapping_quality < args.mapq:
        continue
    # Check if the read overlaps any of the regions in the BED file
    for region in bed_regions:
        if read.reference_name == region[0] and read.reference_start is not None and read.reference_end is not None and read.reference_start < region[2] and read.reference_end > region[1]:
            # Increment the read count for this region
            read_counts[region] += 1
            # Write the read to the output SAM file
            output_samfile.write(read)
            break

# Close the input and output SAM files
input_samfile.close()
output_samfile.close()

# Write the read counts to a file
with open(args.output_counts_file, "w") as outfile:
    for region in bed_regions:
        outfile.write(f"{region[0]}\t{region[1]}\t{region[2]}\t{read_counts[region]}\n")
