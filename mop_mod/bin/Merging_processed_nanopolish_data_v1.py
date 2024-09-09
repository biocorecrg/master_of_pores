#!/usr/bin/env python
desc= "Concatenate file(s) containing median current intensity values per position values."

# Import required libraries:
import argparse
import csv
import gzip
import statistics


def process_dicts(files):
	intensity = dict()
	coverage = dict()
	
	for file in files:
		#Open and read input file:
		with gzip.open(file, mode='rt') as csv_file:
			next(csv_file)
			#csv_reader = csv.reader(csv_file, delimiter='\t')
			for line in csv_file:
				line = line.rstrip()
				row = line.split("\t")
				if (row[1] != "position"):
					#Update both dictionaries:
					#for row in csv_reader:
					key = ','.join(row[0:4])

					#Update intensity dict:
					if key in intensity:
						intensity[key].append(float(row[4]))
					else:
						intensity[key] = [float(row[4])]

					#Update coverage dict:
					if key in coverage:
						coverage[key].append(int(row[5]))
					else:
						coverage[key] = [int(row[5])]

		#Close input file:
		#csv_file.close()
	
	#Data processing of data stored in both dictionaries:
	for key in intensity:
		intensity[key] = statistics.median(intensity[key])
		coverage[key] = sum(coverage[key])
	
	return intensity,coverage

def generate_output(output_file, intensity, coverage):
	#Create output file:
	out_file = output_file+'.tsv'
	f = open(out_file, 'w')
	
	#Print header:
	header=['contig','position','reference_kmer','read_name','median','coverage']
	print('\t'.join(header),file=f)
	
	#Print data stored in both dictionaries:
	for key in intensity:
		splitted_key = key.split(",")
		print('\t'.join(splitted_key), "{:.3f}".format(intensity.get(key)), coverage.get(key), file=f, sep='\t')

	f.close()

def main():
	parser  = argparse.ArgumentParser(description=desc)

	parser.add_argument('-i', '--input', nargs='+' ,help='Input file(s) to process.')
	parser.add_argument('-o', '--output', help='Output filename')

	a = parser.parse_args()

	#Read, parse and merge data from individual eventalign files:
	intensity, coverage = process_dicts(a.input)
	
	#Generate output file with processed data:
	generate_output(a.output, intensity, coverage)

if __name__=='__main__': 
	main()
