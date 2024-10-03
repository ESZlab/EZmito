#!/usr/bin/env python

import argparse
import time
import sys
import os
import pandas as pd
import shutil
import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Nexus import Nexus
import warnings
import numpy as np
import statistics
import pyfiglet

# Remove warnings
warnings.filterwarnings("ignore")


# Start tracking time at the beginning of the script
start_time = time.time()

# Genetic code descriptions
code_help = (
	"Specify the genetic code to use:\n"
	"  2  - Vertebrate Mitochondrial Code\n"
	"  3  - Yeast Mitochondrial Code\n"
	"  4  - Mold, Protozoan, and Coelenterate Mitochondrial Code\n"
	"  5  - Invertebrate Mitochondrial Code\n"
	"  9  - Echinoderm and Flatworm Mitochondrial Code\n"
	"  13 - Ascidian Mitochondrial Code\n"
	"  14 - Alternative Flatworm Mitochondrial Code\n"
	"  16 - Chlorophycean Mitochondrial Code\n"
	"  21 - Trematode Mitochondrial Code\n"
	"  22 - Scenedesmus obliquus Mitochondrial Code\n"
	"  23 - Thraustochytrium Mitochondrial Code\n"
	"  33 - Cephalodiscidae Mitochondrial UAA-Tyr Code"
)


## General functions

# Concatenation
def concatenate(path, strand):
	file_l=[] #empty list to be filled
	for file in os.listdir(path):
		file_l.append(os.path.join(path,file))


	nexi = [(fname, Nexus.Nexus(fname)) for fname in file_l] #function to prepare nexus file to be concatenated

	combined = Nexus.combine(nexi) #combination of nexus files
	output_file = path+'/'+strand+'.nexus'
	combined.write_nexus_data(filename=open(output_file, "w")) #saving the concatenation
	
	return output_file
	
	
# NEXUS to FASTA
def nexus2fasta(nexus):
	output_fasta = nexus.replace('.nexus', '.fasta')
	
	# Create a list to store all the records
	records = []
	
	with open(nexus, "r") as nexus_file:
		# Parse the fasta file and annotate each record with molecule type 'DNA'
		for record in SeqIO.parse(nexus_file, "nexus"):
			record.annotations['molecule_type'] = 'DNA'
			record.seq = record.seq.replace('-','')
			records.append(record)
	
	# Write all the records to the Nexus file
	with open(output_fasta, "w") as fasta_file:
		SeqIO.write(records, fasta_file, "fasta")
#	os.remove(nexus)
	return output_fasta

# FASTA to NEXUS
def fasta2nexus(fasta):
	output_nexus = fasta.replace('.fasta', '.nexus')
	
	# Create a list to store all the records
	records = []
	
	with open(fasta, "r") as fasta_file:
		# Parse the fasta file and annotate each record with molecule type 'DNA'
		for record in SeqIO.parse(fasta_file, "fasta"):
			record.annotations['molecule_type'] = 'DNA'
			records.append(record)
	
	# Write all the records to the Nexus file
	with open(output_nexus, "w") as nexus_file:
		SeqIO.write(records, nexus_file, "nexus")
	os.remove(fasta)
	return output_nexus


# Pad the sequences
def pad_sequence(fasta):
	# Read all sequences from the given file
	sequences = list(SeqIO.parse(fasta, 'fasta'))  # Read all sequences at once
	
	# Determine the maximum sequence length
	max_length = max(len(seq.seq) for seq in sequences)
	
	padded_sequences = []
	for sequence in sequences:
		seq_length = len(sequence.seq)
		if seq_length < max_length:
			# Pad the sequence with gaps ('-') to match the maximum length
			sequence.seq = Seq(str(sequence.seq).ljust(max_length, '-'))
		padded_sequences.append(sequence)
	
	# Write all padded sequences to a new file
	processed_filename = fasta[:fasta.rfind('_')]
	processed_filename = processed_filename+'_aligned.fasta'
	with open(processed_filename, 'w') as outfile:
		SeqIO.write(padded_sequences, outfile, 'fasta')
	return processed_filename

# Remove stop codons
def isStopCodon(filename, geneticCode):

	processed_filename, filename_extension = os.path.splitext(filename)
	processed_filename = processed_filename[:processed_filename.rfind('_')]
	processed_filename = processed_filename+'_checked.fasta'
	
	with open(processed_filename, 'w') as outfile:
		for sequence in SeqIO.parse(filename, 'fasta'):
			translated=sequence.seq.translate(table=geneticCode, to_stop=False) #translating the nucleotide sequence to aa

			if len(sequence.seq) % 3 == 0: #check if the sequence is made of codon or there are more or less nucleotide
				if str(translated).count('*') == 1 and str(translated).find('*')+1 == len(translated) : #check if stop codon is at the end
					sequence.seq=sequence.seq[:-3]					
					print('Warning: Stop codon found at the end of the following sequence: ' + str(filename[filename.rfind('/')+1:])+' of '+ str(sequence.id) +'. It will be automatically removed\n')
				elif str(translated).count('*') > 1:
					print('***CODE ERROR: Multiple stop codons found in the following sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. Check again the data set and re-submit it to the Web Server***\n')
					sys.exit()

				else:
					sequence.seq=sequence.seq


			elif len(sequence.seq) % 3 == 2:
				if str(translated).count('*') >= 1:
					print('***CODE ERROR: Stop codons found in the following truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. Check again the data set and re-submit it to the Web Server***\n')
					sys.exit()

				else:
					print('Warning: Found a truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. It will be automatically adjusted and the analysis is continuing...\n')
					sequence.seq=sequence.seq[:-2]

			elif len(sequence.seq) % 3 == 1:
				if str(translated).count('*') >= 1:
					print('***CODE ERROR: Stop codons found in the following truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. Check again the data set and re-submit it to the Web Server***\n')
					sys.exit()

				else:
					print('Warning: Found a truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. It will be automatically adjusted and the analysis is continuing...\n')
					sequence.seq=sequence.seq[:-1]
		
					
			SeqIO.write(sequence, outfile, "fasta")
	
		return processed_filename


# Replacing _verify alphabet function of Bio.Alphabet
def _verify_alphabet(sequence, alphabet):
	alphabet = set(alphabet) 
	return all(letter in alphabet for letter in sequence)

# Check IUPAC ambiguities
def isIUPAC(filename):
	for sequence in SeqIO.parse(filename, 'fasta'):
		seq = Seq(str(sequence.seq).upper()) #forcing to read the sequence using only the ambiguos DNA alphabet
		test = _verify_alphabet(seq, 'ATCGNRYBDKMHVSW') #testing it
		if test == False: #crashes if it founds a non recognized nucleotide
			print('***CODE ERROR: Found a not recognized nucleotide in the following file ' + str(filename) + ' of the following id' + str(sequence.id)+'. Check your matrix and re-submit the compressed folder to the web server.***\n')
			sys.exit()

# Check the lengths
def check_length(filename):
	seq_length=[]
	for seqrecord in SeqIO.parse(filename ,'fasta'):
		seq_length.append(len(seqrecord.seq))
	avg=statistics.mean(seq_length)
	stdev=statistics.stdev(seq_length)
	treshold = avg+(stdev*2)
	for seqrecord in SeqIO.parse(filename ,'fasta'):
		if len(seqrecord.seq) > treshold or len(seqrecord.seq) < treshold:
			print('Warning: the length of ' + str(seqrecord.id)+ ' in ' + 
			os.path.basename(filename).replace('_noempty.fasta','') + ' differs more than three standard deviation')

# Check the FASTA file
def is_fasta(fasta):
	with open(fasta, "r") as handle:
		fasta_records = list(SeqIO.parse(handle, "fasta"))
		if not fasta_records:
			print("***CODE ERROR. Invalid or empty fasta file detected. Please re-submit a valid file.***")
			sys.exit("Error: No valid FASTA records found.")
		else:
			return fasta_records

# Check the GFF3 file
def is_gff3 (gff_file):
	from BCBio.GFF import GFFExaminer
	examiner = GFFExaminer()
	in_handle = open(gff_file)
	first = in_handle.readline()
	if '##gff-version 3' in first:
		return True
	else:
		print("***CODE ERROR. Invalid GFF3 file detected. Please re-submit a valid file. Check the GFF3 format at: http://www.ensembl.org/info/website/upload/gff3.html***")
		sys.exit("Error: No valid GFF3 file.")
	in_handle.close()

# Check if there are any repeated IDs in the fasta file
def is_duplicated(filename):
	seen_ids = set()

	for record in SeqIO.parse(filename, 'fasta'):
		if record.id in seen_ids:
			print(f"***CODE ERROR. Duplicated ID: {record.id} found in {filename}. Check your matrix and re-submit.***")
			sys.exit(f"Error: Duplicated ID '{record.id}' found.")

		seen_ids.add(record.id)

#	print("Quality check | No duplicated IDs found. Continuing the analysis...")

# Check if the input file is a valid fasta and process the sequences
def check_fasta(filename, outdir):

	fasta_counter = 0
	
	processed_filename, filename_extension = os.path.splitext(filename)
	processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_noempty.fasta')

	fasta_records = is_fasta(filename)
	if fasta_records:
		# Writing the cleaned file
		with open(processed_filename, 'w') as cleaned_file:
			for record in fasta_records:
				record.description = record.description.replace(' ', '_')
				SeqIO.write(record, cleaned_file, 'fasta')
				fasta_counter += 1

#	if fasta_counter:
#		print("Quality check | Fasta file processed without errors. Continuing...")
	#look for duplicated ids
	is_duplicated(processed_filename)
	return processed_filename

# Remove gaps ('-') from sequences
def remove_gaps(filename):
	processed_filename, filename_extension = os.path.splitext(filename)
	processed_filename = processed_filename[:processed_filename.rfind('_')]
	processed_filename = processed_filename+'_degapped.fasta'

	with open(processed_filename, 'w') as outfile:
		for record in SeqIO.parse(filename, 'fasta'):
			# Check for gaps and remove them
			if '-' in str(record.seq):
				print(f"Warning: Found gap (-) in {record.id}. Removing gap...")
				record.seq = Seq(str(record.seq).replace("-", ""))
			# Write each record to the file
			SeqIO.write(record, outfile, 'fasta')

#	print("Quality check | Gaps removed from sequences. File is ready.")
	
	return processed_filename
	
def replace_dir(directory):
	# If the directory exists, remove it
	if os.path.exists(directory):
		shutil.rmtree(directory)
	
	# Create the new directory
	os.makedirs(directory, exist_ok=True)

# Placeholder function for each subcommand
def ez_circular_subcommand(args):
	from Bio.SeqRecord import SeqRecord
	import pybedtools
	import re

	fasta_file = args.input
	bed_file = args.bed
	output_fasta = os.path.join(args.outdir, 'output.fasta')
	output_bed = os.path.join(args.outdir, 'output.bed')
	gene_name = args.start
	linear = '' if args.feature == 'circular' else 'linear'
	
	
	print(pyfiglet.figlet_format("EZcircular")) 
	
	print(f"Running EZcircular with the following parameters:\ninput: {fasta_file}\nbed: {bed_file}\noutdir: {args.outdir}\nfeature: {linear}\nstarting gene: {gene_name}\n\n")

	# Remove parentheses from gene name if present
	if '(' in gene_name or ')' in gene_name:
		gene_name = re.sub(r'\([^)]*\)', '', gene_name)
	
	def make_circular(fasta_file, bed_file):
		# Modify the FASTA file
		with open(fasta_file + '_linear.fa', 'w') as circularized:
			for record in SeqIO.parse(fasta_file, 'fasta'):
				seq_len = len(record.seq)
				record.seq = record.seq + 'N' * 100
				SeqIO.write(record, circularized, 'fasta')
				break

		# Modify the BED file
		columns = ["chrom", "start", "end", "name", "score", "strand"]
		bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=columns)
		bed_df.loc[len(bed_df.index)] = [bed_df['chrom'][0], seq_len, seq_len + 100, 
			'### Artifact of circularization of EZcircularize ###', '', '']
		bed_df.to_csv(bed_file + '_linear.bed', index=False, header=None, sep='\t')
		return fasta_file + '_linear.fa', bed_file + '_linear.bed'
	
	def parse_bed(bed_file, gene_name):
		try:
			bedtool = pybedtools.BedTool(bed_file)
			interval = bedtool.filter(lambda x: gene_name.lower() in x.name.lower())[0]
			return int(interval.start), int(interval.end)
		except Exception as e:
			print(f"Error parsing BED file: {e}")
			return None

	def reorder_and_write_fasta(gene_coordinates, fasta_file, output_fasta, gene_name):
		# Ensure output directory exists
		replace_dir(os.path.dirname(output_fasta))
		
		if gene_coordinates is None:
			print(f"Error: Unable to extract coordinates for gene {gene_name} from BED.")
		else:
			gene_start, gene_end = gene_coordinates
			gene_record = SeqRecord(Seq(""), id="", description="")
			for record in SeqIO.parse(fasta_file, "fasta"):
				gene_record.seq += record.seq[gene_start:]  # Adjust coordinates to 0-based index
				gene_record.seq += record.seq[:gene_start]
				gene_record.id = record.id + '_arranged_from_' + gene_name
			with open(output_fasta, "w") as output_handle:
				SeqIO.write(gene_record, output_handle, "fasta")
		return gene_record.id


	
	def write_bed(bed_file, output_bed, gene_coordinates, gene_record_id):
		columns = ["chrom", "start", "end", "name", "score", "strand"]
		bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=columns)
		length_seq = len(next(SeqIO.parse(fasta_file, "fasta")).seq)
		bed_df_first = bed_df[bed_df["start"] >= gene_coordinates[0]].copy()
		bed_df_first["start"] -= gene_coordinates[0]
		bed_df_first["end"] -= gene_coordinates[0]
		bed_df_second = bed_df[bed_df["start"] < gene_coordinates[0]].copy()
		to_add = length_seq - int(bed_df['end'].iloc[-1]) + bed_df_first['end'].iloc[-1]
		bed_df_second["start"] += to_add
		bed_df_second["end"] += to_add
		bed_df_ordered = pd.concat([bed_df_first, bed_df_second])
		bed_df_ordered["chrom"] = gene_record_id
		bed_df_ordered.to_csv(output_bed, sep='\t', index=False, header=None)

	# Apply circularization if needed
	if linear == 'linear':
		fasta_file, bed_file = make_circular(fasta_file, bed_file)

	# Parse the BED file and reorder sequences
	gene_coordinates = parse_bed(bed_file, gene_name)
	gene_record_id = reorder_and_write_fasta(gene_coordinates, fasta_file, output_fasta, gene_name)
	write_bed(bed_file, output_bed, gene_coordinates, gene_record_id)

	# Cleanup temporary files if created
	if linear == 'linear':
		os.remove(fasta_file)
		os.remove(bed_file)

	warnings.resetwarnings()

	

def ez_codon_subcommand(args):

	from cai2.CAI import RSCU
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.patches as patches


	def ezcodon_main(path, genetic_code, strand, outdir):
		for filename in os.listdir(path):
			filename = os.path.join(path,filename)
			processed_file = check_fasta(filename, outdir)   # Check if the fasta file is valid and free of duplicates
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)  # Check and remove gaps
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			padded_file = pad_sequence(clean_file)
			os.remove(clean_file)
			nexus_file = fasta2nexus(padded_file)
		concatenated = concatenate(outdir, strand)
		return concatenated
		
	def CodonUsage(fasta, folder, geneticCode, strand):
		df_codon=pd.DataFrame(columns=['Species', 'Codon', 'RSCU']) #create and empty df
		df_aa=pd.DataFrame(columns=['Species', 'AA', 'Freq'])
	
		for sequence in SeqIO.parse(fasta, 'fasta'): #for each entry in the fasta combined sequences
			#RSCU
			whole=str('')
			d={}
			codons = [sequence.seq[i:i+3] for i in range(0, len(sequence.seq), 3)] #divide the sequence in triplets
			for codon in codons: #loop to parse single codons
				seq= Seq(str(codon).upper()).translate(table=geneticCode, to_stop=False)
				test = _verify_alphabet(seq, 'ACDEFGHIKLMNPQRSTVWY')
				if test == True:
					string=str(codon)
					whole += string
					d[str(codon).upper()]=sequence.id
				else:
					pass
					print(f'Warning: the codon ({codon}) in {sequence.id} produces an ambiguos aminoacid and it will be excluded from the analysis')
			seqlist=[]
			seqlist.append(whole)
			rscu=RSCU(seqlist, genetic_code=int(geneticCode))
			rscu=pd.DataFrame(rscu.items(), columns=['Codon', 'RSCU']) #calculates the RSCU for each seq
			df= pd.DataFrame(d.items(),columns=['Codon','Species']) #dictionary converted to a dataframe
			df= df.merge(rscu, how='left', on='Codon')
			df_codon = pd.concat([df_codon, df], ignore_index=True, sort=False) #add the current df to the empty one
			#AAfreq
			d_aa={'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0 ,'V':0, 'W':0, 'Y':0}
			translated=sequence.seq.translate(table=geneticCode, to_stop=False)
			for aa in translated:
				for k,v in d_aa.items():
					if aa == k:
						d_aa[k]=v+1
			for k,v in d_aa.items():
				d_aa[k]=v/len(translated)*100
			df= pd.DataFrame(d_aa.items(),columns=['AA','Freq'])
			df['Species']=str(sequence.id)
			df_aa = pd.concat([df_aa, df], ignore_index=True, sort=False)

		#RSCU

		df_codon=df_codon.reset_index() #reset index
		df_codon=df_codon[['Species','Codon','RSCU']] #discard the 1st column (it was the previous index)
		

		df_codon['AA']='' #add a new empty AA column
		for i in range(0,len(df_codon)): #for each triplet
			triplet=str(Seq(df_codon.iloc[i][1]).translate(table=geneticCode, to_stop=False)) #translate it to AA
			df_codon.iloc[i,3]=triplet #Add this result to the new column o
		
		output_RSCU = os.path.join(folder,strand+'_RSCU.csv')
		df_codon.to_csv(output_RSCU)


		#AAfreq
		df_aa=df_aa.reset_index()
		df_aa=df_aa[['Species','AA','Freq']]
		output_AAfreq =  os.path.join(folder,strand+'_AAfreq.csv')
		df_aa.to_csv(output_AAfreq)
		
		return output_RSCU, output_AAfreq
		
	def AA_freq_plot(table, strand, outdir):
	# Load data
		if os.path.exists(table):
			df = pd.read_csv(table).drop(columns=['Unnamed: 0'])
			nsp = df['Species'].nunique()

			# For amino acid frequency plot
			output_fig = os.path.join(outdir, f'{strand}_Aminoacid_frequency.pdf')
			
			if nsp <= 20:
				# Line plot for species with <= 20 unique species
				plt.figure(figsize=(20, 15))
				for species in df['Species'].unique():
					species_data = df[df['Species'] == species]
					plt.plot(species_data['AA'], species_data['Freq'], marker='o', label=species)
					
				plt.title(f'{strand} strand Amino Acid Frequency', fontsize=15)
				plt.xlabel('Amino Acids', fontsize=12)
				plt.ylabel('Frequency', fontsize=12)
				plt.legend(title="Species")
				plt.tight_layout()
				
				# Enable the grid
				plt.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5, color='#cccccc')
				plt.gca().set_axisbelow(True)
				plt.savefig(output_fig, dpi=300)
				plt.close()  # Close the plot to avoid display

			else:
				# Box plot for species with > 20 unique species
				plt.figure(figsize=(20, 15))
				df.boxplot(column='Freq', by='AA', grid=False, showfliers=False)
				plt.title(f'{strand} strand Amino Acid Frequency', fontsize=15)
				plt.xlabel('Amino Acids', fontsize=12)
				plt.ylabel('Frequency (%)', fontsize=12)
				plt.suptitle('')  # Suppress the automatic 'by' title
				plt.tight_layout()
				
				# Enable the grid
				plt.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5, color='#cccccc')
				plt.gca().set_axisbelow(True)
				plt.savefig(output_fig, dpi=300)
				plt.close()  # Close the plot to avoid display
				
	def RSCU_plot(table, strand, outdir):
		if os.path.exists(table):
			codon = pd.read_csv(table).drop(columns=['Unnamed: 0'])
			species_groups = codon.groupby('Species')

			output_fig = os.path.join(outdir, f'{strand}_strand_RSCU.pdf')

			# Create a PDF file to store all plots
			with PdfPages(output_fig) as pdf:

				for name, df in species_groups:
					df['Col'] = df.groupby('AA').cumcount()

					# Sort the dataframe by AA alphabetically
					df = df.sort_values(by='AA')

					# Create a figure with two subplots: one for bar plot and one for tile plot
					fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [8, 2]}, figsize=(15, 10))

					# Get unique AA and assign positions after sorting
					unique_AA = df['AA'].unique()
					x_positions = np.arange(len(unique_AA))

					# Create a dictionary to map AA to x-axis positions
					aa_to_position = {aa: i for i, aa in enumerate(unique_AA)}

					# Initialize bottom to 0 for stacking bars
					bottoms = {aa: 0 for aa in unique_AA}
					
					# Define the colormap with a consistent number of unique colors
					colors = plt.get_cmap('Dark2', len(df['Col'].unique()))

					# Plot the stacked bars for codons on ax1
					for i, row in df.iterrows():
						x = aa_to_position[row['AA']]  # Get x position from the AA
						ax1.bar(x, row['RSCU'], bottom=bottoms[row['AA']], 
								color=colors(row['Col']), edgecolor='black', width=0.8)
						
						# Update the bottom position for stacking
						bottoms[row['AA']] += row['RSCU']

					# Set x-axis labels as amino acids
					ax1.set_xticks(x_positions)
					ax1.set_xticklabels(unique_AA, position=(1,-.05))  # Assign amino acids (AA) as labels on the x-axis

					# Add titles and labels
					ax1.set_title(f"{name} {strand} strand RSCU")
					ax1.set_ylabel('RSCU')

					# Enable the grid on ax1
					ax1.grid(True, which='both', axis='y', linestyle='-', linewidth=0.5, color='#cccccc')
					ax1.grid(True, which='both', axis='x', linestyle='-', linewidth=0.5, color='#cccccc')

					# Ensure the grid is in the background
					ax1.set_axisbelow(True)

					# Adjust layout
					ax1.set_xlim(-0.5, len(df['AA'].unique()) - 0.5)

					# Tile plot (Codons) on ax2
					amino_acids = df['AA'].unique()

					# Reverse the order for y-axis and iteration for proper stacking
					col_values_reversed = df['Col'].max() - df['Col']

					for idx, row in df.iterrows():
						aa = row['AA']
						codon = row['Codon']
						col_value = col_values_reversed[idx]  # Use the reversed value
						aa_idx = list(amino_acids).index(aa)
						codon_idx = col_value  # Use the reversed index
						
						# Draw the rectangle tile for the codon
						ax2.add_patch(patches.Rectangle((aa_idx, codon_idx), 1, 1, 
														fc=colors(row['Col']), ec='black'))

						# Add codon label
						ax2.text(aa_idx + 0.5, codon_idx + 0.5, codon, ha='center', va='center', 
								color='white', fontweight='bold', size=7)

					# Set the axis limits and labels
					ax2.axes.get_yaxis().set_visible(False)
					ax2.axes.get_xaxis().set_visible(False)

					# Set axis limits
					ax2.set_xlim(0, len(amino_acids))
					ax2.set_ylim(0, df['Col'].max() + 1)

					# Hide axis spines and ticks for ax2
					ax2.spines[:].set_visible(False)
					ax2.tick_params(left=False, bottom=False)

					# Remove grid for ax2
					ax2.grid(False)

					# Save the figure to the PDF
					pdf.savefig(fig)
					plt.close()


	genetic_code = args.code
	outdir = args.outdir
	
	replace_dir(outdir)
	
	tmp = os.path.join(outdir, 'tmp')
	replace_dir(tmp)
	
	plots = os.path.join(outdir, 'plots')
	replace_dir(plots)
	
	tables = os.path.join(outdir, 'tables')
	replace_dir(tables)
	
	print(pyfiglet.figlet_format("EZcodon")) 
	
	if args.heavy and not args.light:
		J_path = args.heavy
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		print(f"Running EZcodon with the following parameters:\nheavy chain fasta files: {J_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and not args.heavy:
		N_path = args.light
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		print(f"Running EZcodon with the following parameters:\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and args.heavy:
		J_path = args.heavy
		N_path = args.light
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		JN_tmp = os.path.join(tmp, 'JN')
		replace_dir(JN_tmp)
		print(f"Running EZcodon with the following parameters:\nheavy chain fasta files: {J_path}\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	# Ensure at least one of --heavy or --light is provided
	if not args.heavy and not args.light:
		parser.error("At least one of --heavy (-J) or --light (-N) must be provided.")
	
	
	
	
	if args.heavy and not args.light:
		concatenated = ezcodon_main(J_path, genetic_code, 'J', J_tmp)
		final_file = nexus2fasta(concatenated)
		output_RSCU, output_AAfreq = CodonUsage(final_file, tables, genetic_code, 'J')
		AA_freq_plot(output_AAfreq, 'J', plots)
		RSCU_plot(output_RSCU, 'J', plots)
		
	if args.light and not args.heavy:
		concatenated = ezcodon_main(N_path, genetic_code, 'N', N_tmp)
		final_file = nexus2fasta(concatenated)
		output_RSCU, output_AAfreq = CodonUsage(final_file, tables, genetic_code, 'N')
		AA_freq_plot(output_AAfreq, 'N', plots)
		RSCU_plot(output_RSCU, 'N', plots)
		
	if args.light and args.heavy:
	
		concatenatedJ = ezcodon_main(J_path, genetic_code, 'J', J_tmp)
		concatenatedN = ezcodon_main(N_path, genetic_code, 'N', N_tmp)
		
		# copy the files to the JN folder
		src_files = os.listdir(J_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(J_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'J.nexus':
				shutil.copy(full_file_name, JN_tmp)
		src_files = os.listdir(N_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(N_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'N.nexus':
				shutil.copy(full_file_name, JN_tmp)
		concatenated = concatenate(JN_tmp, 'JN')
		final_file = nexus2fasta(concatenated)
		output_RSCU, output_AAfreq = CodonUsage(final_file, tables, genetic_code, 'JN')
		AA_freq_plot(output_AAfreq, 'JN', plots)
		RSCU_plot(output_RSCU, 'JN', plots)
		
	shutil.rmtree(tmp)

def ez_map_subcommand(args):
	from pycirclize.parser import Gff
	
	gff_file = args.gff
	colorJ = args.colorJ
	colorN = args.colorN
	feature = args.feature
	outdir = args.outdir
	
	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZmap")) 
	
	print(f"Running EZmap with the following parameters:\nmt_feature: {feature}\nGFF file: {gff_file}\ncolors: {colorJ}, {colorN}\noutdir: {outdir}\n\n")
	
		
	def plot_linear(gff_file, colorJ, colorN, outdir):
		gff = Gff(gff_file)

		f_cds_feats = gff.extract_features(["CDS","tRNA","rRNA"], target_strand=1)
		r_cds_feats = gff.extract_features(["CDS","tRNA","rRNA"], target_strand=-1)

		AT_feat = gff.extract_features("sequence_feature")

		plt.figure(figsize=(30, 5))
		plt.ylim(-3,6)
		plt.xlim(0,gff.range_size
		+1)

		plt.arrow(0, 0, gff.range_size+1, 0, head_width= 0, head_length= 0, width= 0.01, fc='black', ec='black', lw = 1 )

		cnt = 1
		for feat in f_cds_feats:
			start, end = int(str(feat.location.start)), int(str(feat.location.end))
			label = feat.qualifiers.get("product", [""])[0]
			plt.arrow(start, 0, end - start, 0, head_width= 0.4, head_length= (end - start) *0.25, width= 0.25, length_includes_head= True, fc=colorJ, ec='black', lw = 1 )
			if len(label) > 10:
				plt.text((start + end) / 2, .5, label, ha='center',  va = 'baseline', fontsize=10, rotation=90)
			else:
				if cnt % 2 == 0:
					plt.text((start + end) / 2, .5, label, ha='center',  va = 'baseline', fontsize=10, rotation=90)
				else:
					plt.text((start + end) / 2, -2, label, ha='center', fontsize=10, rotation=90)
				cnt +=1

		for feat in r_cds_feats:
			end, start = int(str(feat.location.start)), int(str(feat.location.end))
			label = feat.qualifiers.get("product", [""])[0]
			plt.arrow(start, 0, end - start, 0, head_width= 0.4, head_length= (start - end) *0.25, width= 0.25, length_includes_head= True, fc=colorN, ec='black', lw = 1)
			if len(label) > 10:
				plt.text((start + end) / 2, .5, label, ha='center',  va = 'baseline', fontsize=10, rotation=90)
			else:
				if cnt % 2 == 0:
					plt.text((start + end) / 2, .5, label, ha='center',  va = 'baseline', fontsize=10, rotation=90)
				else:
					plt.text((start + end) / 2, -2, label, ha='center', fontsize=10, rotation=90)
				cnt +=1

		for feat in AT_feat:
			start, end = int(str(feat.location.start)), int(str(feat.location.end))
			label = 'A+T rich region / Control region'
			plt.arrow(start, 0, end - start, 0,head_width= 0.4, head_length= (end - start) *0.25, width= 0.25, length_includes_head= True, fc=colorJ, ec='black', lw = 1)
			if len(label) > 10:
				plt.text((start + end) / 2, .5, label, ha='center',  va = 'baseline', fontsize=10, rotation=90)
			else:
				if cnt % 2 == 0:
					plt.text((start + end) / 2, .5, label, ha='center',  va = 'baseline', fontsize=10, rotation=90)
				else:
					plt.text((start + end) / 2, -2, label, ha='center', fontsize=10, rotation=90)
				cnt +=1


		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['left'].set_visible(False)
		plt.gca().yaxis.set_visible(False)
		plt.xticks(np.arange(0, gff.range_size+1, 500),rotation=90)
		plt.xlabel('Genomic position (bp)')
		plt.savefig(os.path.join(outdir, 'mt_linear_output.pdf'), bbox_inches='tight')


	def plot_circular(gff_file, colorJ, colorN, outdir):
		from pycirclize import Circos
		from pycirclize.parser import Gff

		# Load GFF file
		gff = Gff(gff_file)
		circos = Circos(sectors={gff.name: gff.range_size})
		circos.text(gff.target_seqid, size=15)

		sector = circos.sectors[0]
		cds_track = sector.add_track((80, 100))
		cds_track.axis(fc="white", ec="none")




		# Plot forward & reverse CDS
		region = gff.extract_features("region")

		f_cds_feats = gff.extract_features(["CDS","tRNA","rRNA"], target_strand=1)
		r_cds_feats = gff.extract_features(["CDS","tRNA","rRNA"], target_strand=-1)

		AT_feat = gff.extract_features("sequence_feature")


		# Plot forward & reverse CDS with increased arrow sizes

		cds_track.genomic_features(
			region,
			plotstyle="box",
			r_lim=(90, 91),  # Keep within valid range, larger radial range for bigger arrows
			fc="black",
			lw=.5  # Increase line width for thicker arrows
		)

		cds_track.genomic_features(
			f_cds_feats,
			plotstyle="arrow",
			r_lim=(80, 100),  # Keep within valid range, larger radial range for bigger arrows
			fc=colorJ,
			lw=1.0  # Increase line width for thicker arrows
		)

		cds_track.genomic_features(
			r_cds_feats,
			plotstyle="arrow",
			r_lim=(80, 100),  # Keep within valid range for reverse strand
			fc=colorN,
			lw=1.0  # Increase line width for thicker arrows
		)

		cds_track.genomic_features(
			AT_feat,
			plotstyle="arrow",
			r_lim=(80, 100),  # Keep within valid range, larger radial range for bigger arrows
			fc=colorJ,
			lw=1.0  # Increase line width for thicker arrows
		)



		# Extract CDS product labels
		pos_list, labels = [], []
		for feat in gff.extract_features(["CDS", 'tRNA', 'rRNA','sequence_feature']):
			start, end = int(str(feat.location.end)), int(str(feat.location.start))
			pos = (start + end) / 2
			if  'A+T' in feat.type or 'AT region' in feat.type:
				label = 'A+T rich region / Control region'
			else:
				label = feat.qualifiers.get("product", [""])[0]
			if label == "" or label.startswith("hypothetical"):
				continue
			pos_list.append(pos)
			labels.append(label)


		# Plot CDS product labels on outer position
		cds_track.xticks(
			pos_list,
			labels,
			label_orientation="vertical",
			show_bottom_line=True,
			label_size=6,
			line_kws=dict(ec="white"),
		)
		# Plot xticks & intervals on inner position
		cds_track.xticks_by_interval(
			interval=500,
			outer=False,
			show_bottom_line=True,
			label_formatter= lambda v: f"{int(v)} bp", #lambda v: f"{v/ 1000:.1f} Kb",
			label_orientation="vertical",
			line_kws=dict(ec="black"),
		)

		fig = circos.plotfig()
		fig.savefig(f'{outdir}/circular_plot.pdf', bbox_inches='tight')  # Save the plot as a PNG file
		plt.close(fig)
		
	is_gff3(gff_file)
	
	if feature == 'linear':
		plot_linear(gff_file, colorJ, colorN, outdir)
	
	else:
		plot_circular(gff_file, colorJ, colorN, outdir)
	

def ez_mix_subcommand(args):
	import matplotlib.cm as cm
	
	fasta_file = args.input
	length = args.length
	identity = args.identity*100
	outdir = args.outdir
	blastn = os.path.join(args.blastn, 'blastn')
	
	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZmix")) 
	
	print(f"Running EZmix with the following parameters:\ninput: {fasta_file}\noutdir: {outdir}\nidentity: {identity}\nlength: {length}\npath to blastn: {blastn}\n\n")
	
	def run_blast(query, subject, outdir):
			output_file = os.path.join(outdir, "blout")
			command = f"{blastn} -query {query} -subject {subject} -outfmt '6 qseqid qstart qend sseqid sstart send pident length evalue bitscore slen qlen' -evalue 0.5 > {output_file}"
			subprocess.run(command, shell=True)
			return output_file if os.path.exists(output_file) and os.path.getsize(output_file) > 0 else None

	# Function to parse BLAST output
	def parse_blast_output(blast_file, min_length, min_similarity):
		# Define column names based on the BLAST output format
		columns = ['qseqid', 'qstart', 'qend', 'sseqid', 'sstart', 'send', 'pc', 'length', 'evalue', 'bitscore', 'slen', 'qlen']
		
		# Read the BLAST output into a pandas DataFrame
		df = pd.read_csv(blast_file, sep="\t", header=None, names=columns)
		
		# Filter rows based on minimum similarity (pc) and length
		filtered_df = df[(df['pc'] >= min_similarity) & (df['length'] >= min_length)]
		
		# Select only the relevant columns (you can modify this based on what you need)
		hits_df = filtered_df[['qseqid','qstart', 'qend', 'sseqid' , 'sstart', 'send', 'pc']]

		# Convert the filtered DataFrame to a list of dictionaries
		hits = hits_df.to_dict(orient='records')
		
		return hits
		
	# Function to create a PDF plot
	def create_plot(names, df, max_length, output_pdf, min_similarity, min_length):
		plt.figure(figsize=(10,7))  # A4 size in inches
		plt.title('EZmix output (min length: '+ str(min_length) +  'bp ; min percent identity: ' +  str(min_similarity) + '%)')
		plt.xlabel('Assembly, bp')
		plt.xlim(-100, max_length+100)
		plt.ylim(-1, len(names))

		
		n = len(df.qseqid.unique())
		cmap = cm.get_cmap('Accent')
		colors = [cmap(i) for i in range(n)]

		my_colors = {key:value for key, value in zip(df.qseqid.unique(), colors)}

		for index, row in df.iterrows():
   		# Plot horizontal bars
			plt.plot([0,row['length']],  [row['qseqid'],row['qseqid']], color='black')
			if row['qstart'] != '':
				plt.plot([row['qstart'],row['qend']],  [row['qseqid'],row['qseqid']], color='black', lw=5)
				plt.plot([row['sstart'],row['send']],  [row['sseqid'],row['sseqid']], color='black', lw=5)
				plt.plot([(row['qstart'] + row['qend']) / 2, (row['sstart'] + row['send']) / 2], 
						 [row['qseqid'], row['sseqid']], lw=2, color=my_colors[row['qseqid']])

		plt.xticks(np.arange(0, max_length+1, 500), rotation=75)
		plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5, color='lightgray')
		plt.savefig(output_pdf)
		plt.close()		
		

	processed_file = check_fasta(fasta_file, outdir)   # Check if the fasta file is valid and free of duplicates
	checked_file = remove_gaps(processed_file)  # Check and remove gaps
	
	os.remove(processed_file)
	
	
	sequences = list(SeqIO.parse(checked_file, 'fasta'))
	seq_names = [record.id for record in sequences]  # Truncate names to 10 characters
	seq_lengths = [len(record.seq) for record in sequences]
	lengths_df = pd.DataFrame({name:length for name, length in zip(seq_names, seq_lengths)}.items(), columns=['qseqid','length'])
	lengths_df['qstart'] = ''
	max_length = max(seq_lengths)
	
	hits = []
	for i in range(len(sequences) - 1):
		for j in range(i + 1, len(sequences)):
			query_file = os.path.join(outdir, f"primo_{i}.fasta")
			subject_file = os.path.join(outdir, f"secondo_{j}.fasta")
			
			SeqIO.write(sequences[i], query_file, 'fasta')
			SeqIO.write(sequences[j], subject_file, 'fasta')
			
			blast_output = run_blast(query_file, subject_file, outdir)
			if blast_output is None:
				continue
			else:
				hits += parse_blast_output(blast_output, length, identity)
			
			os.remove(query_file)
			os.remove(subject_file)
			if os.path.exists(blast_output):
				os.remove(blast_output)
				
	#### Attenzione, non crea un df correttamente. ci sono sseqid che mancano nei qseqid. quindi aggiungili a mano se non li trovi (missing obj). poi c'è da capire come non far cadere il primo qseqid nello 0 dell'y axis 			
				
	hits_df = pd.DataFrame(hits)
	
	hits_df = pd.concat([hits_df, lengths_df], ignore_index=True)
	# Create the plot
	output_pdf = os.path.join(outdir, f"{os.path.basename(fasta_file)}_output.pdf")
	create_plot(seq_names, hits_df, max_length, output_pdf, identity, length)
	print(f"Plot saved to {output_pdf}")
	

def ez_pipe_subcommand(args):

	from Bio.SeqRecord import SeqRecord
	from itaxotools.pygblocks import compute_mask, trim_sequence, Options
	from Bio import  AlignIO
	import re
	
	# Function to write the fixed initial part of the config file
	def write_initial_cfg(path):
		with open(f'{path}/partition_finder.cfg', 'w') as cfg_file:
				cfg_file.write('## ALIGNMENT FILE ##\n')
				cfg_file.write('alignment = infile.phy;\n\n')
				cfg_file.write('## BRANCHLENGTHS: linked | unlinked ##\n')
				cfg_file.write('branchlengths = linked;\n\n')
				cfg_file.write('## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##\n')
				cfg_file.write('models = mrbayes;\n\n')
				cfg_file.write('# MODEL SELECCTION: AIC | AICc | BIC #\n')
				cfg_file.write('model_selection = aicc;\n\n')
				cfg_file.write('## DATA BLOCKS: see manual for how to define ##\n')
				cfg_file.write('[data_blocks]\n')

# Function to write the final part of the config file
	def write_final_cfg(path):
		with open(f'{path}/partition_finder.cfg', 'a') as cfg_file:
				cfg_file.write('\n## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##\n')
				cfg_file.write('[schemes]\n')
				cfg_file.write('search = greedy;\n')
	
	
	def write_charsets(nexus):

		charset_file, filename_extension = os.path.splitext(nexus)
		charset_file = charset_file + '.charsets'
		new_content = ""

		# Open the input file for reading
		with open(nexus, 'r') as f:
				# Read the content and split to find the relevant portion after 'begin sets;'
				charsets = f.read().split('begin sets;')[1]
				
				# Split the content by semicolon to process each part separately
				path = charsets.split(';')
				
				for p in path:
						if 'charpartition combined' not in p:
								# Extract the gene name from the path
								gene_name = p.split('/')[-1].split('.')[0]
								
								gene_name = gene_name.split('_checked')[0]
								
								# Define the regex pattern
								pattern = r"(charset\s+)([\S]+\.nexus)(\s+=\s+\d+-\d+)"
								
								# Replace the matched path while keeping the rest intact and prepend with a tab
								new_text = re.sub(pattern, fr'\tcharset {gene_name}.nexus\3', p)
								
								# Append the modified string to new_content
								new_content += new_text + ';'
						else:
								break

		# Write the new content to the output file
		with open(charset_file, 'w') as f_out:
				f_out.write('BEGIN ASSUMPTIONS;')  # Write the 'begin sets;' line back
				f_out.write(new_content)
				f_out.write('\nEND;')  # Write the 'end sets;' line back 
				
		return charset_file	
				
	# NEXUS to PHYLIP
	def nexus2phylip(nexus):
		output_phylip = nexus.replace('.nexus', '.phy')
		
		# Create a list to store all the records
		records = []
		
		with open(nexus, "r") as nexus_file:
			# Parse the fasta file and annotate each record with molecule type 'DNA'
			for record in SeqIO.parse(nexus_file, "nexus"):
				record.annotations['molecule_type'] = 'DNA'
				records.append(record)
		
		# Write all the records to the Nexus file
		with open(output_phylip, "w") as phylip_file:
			SeqIO.write(records, phylip_file, "phylip-relaxed")
		os.remove(nexus)
		return output_phylip

	def remove_third_codon_position(filename, outdir):
		
		processed_filename, filename_extension = os.path.splitext(filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_twopos.fasta')
	
		with open(processed_filename, "w") as output_handle:
			for record in SeqIO.parse(filename, "fasta"):
				sequence = str(record.seq)
				trimmed_seq = ''.join([sequence[i:i+2] for i in range(0, len(sequence), 3)])
				trimmed_record = SeqRecord(Seq(trimmed_seq), id=record.id)
				SeqIO.write(trimmed_record, output_handle, "fasta")
				
		return processed_filename
			

	def Gblocks(filename_aa, filename_nt, outdir):
		
		
#		command = f'Gblocks {filename} -t=c -b3=4 -e=gblocks'
#		
#		subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#		
#		output_file = filename+'-gb'
#		htm_file = filename+'-gb.htm'
#		
#		renamed_file, filename_extension = os.path.splitext(filename)
#		renamed_file = os.path.join(outdir,os.path.basename(renamed_file)+'_trimmed.fasta')
#		
#		shutil.copy(output_file, renamed_file)
#		os.remove(htm_file)
#		os.remove(output_file)

		trimmed_file, filename_extension = os.path.splitext(filename_nt)
		trimmed_file = os.path.join(outdir,os.path.basename(trimmed_file)+'_trimmed.fasta')
		
		sequences = AlignIO.read(filename_aa, 'fasta')
		
		options = Options(
			    IS = (len(sequences)*0.5)+1,
			    FS = len(sequences)*0.85,
			    CP = 4,
			    BL1 = 10,
			    BL2 = 10,
			    GT = 0,
			    GC = '-'
			)
			
		mask = compute_mask(sequences, options, log=False)
		mask = ''.join([char * 3 for char in mask])
		
		with open(trimmed_file, 'w') as outfile:
			for record in SeqIO.parse(filename_nt, 'fasta'):
				record.seq = trim_sequence(record.seq, mask)
				SeqIO.write(record, outfile, 'fasta')
		
		return trimmed_file
		

	def alignment(aa_filename, nt_filename, outdir): 
	
		tmp_filename, filename_extension = os.path.splitext(aa_filename)
		tmp_filename = os.path.join(outdir,os.path.basename(tmp_filename)+'_aligned.fasta')
		
		processed_filename, filename_extension = os.path.splitext(aa_filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_ntaligned.fasta')
		
		command = f'mafft --quiet {aa_filename} > {tmp_filename}'
		subprocess.run(command, shell=True)
	
		with open(processed_filename,'w') as nt_alignment:
			seq_length=[]
			for al_sequence in SeqIO.parse(tmp_filename,'fasta'):
				seq_length.append(len(al_sequence.seq))
				for nt_sequence in SeqIO.parse(nt_filename, 'fasta'):
					if al_sequence.id == nt_sequence.id:
						seq='' #it is 'refilled' with adjacent codon
						c=0
						for x in al_sequence.seq:
								
							if x != '-':
								x = 3*int(c+1)-3
								codon=nt_sequence.seq[x:x+3]
								seq=seq+str(codon)
								c+=1
							else:
								gap='---'
								seq=seq+gap
									
						nt_alignment.write('>'+nt_sequence.id+'\n'+seq+'\n')
						break

		return processed_filename, tmp_filename

	def translate(filename, genetic_code, outdir):
		processed_filename, filename_extension = os.path.splitext(filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_translated.fasta')
		
		with open(processed_filename, 'w') as outfile:
			for record in SeqIO.parse(filename, 'fasta'):
				record.seq = record.seq.translate(table=genetic_code, to_stop=False)
				SeqIO.write(record, outfile, 'fasta')
				
		return processed_filename
	
	def ezpipe_main(path, genetic_code, outdir, positions):
		for filename in os.listdir(path):
			filename = os.path.join(path,filename)
			processed_file = check_fasta(filename, outdir)   # Check if the fasta file is valid and free of duplicates
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)  # Check and remove gaps
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			translated_file = translate(clean_file, genetic_code, outdir)
			alignment_nt_file, alignment_aa_file = alignment(translated_file, clean_file, outdir)
			gblocked = Gblocks(alignment_aa_file, alignment_nt_file, outdir)
			os.remove(alignment_aa_file)
			os.remove(alignment_nt_file)
			os.remove(translated_file)
			os.remove(clean_file)
			if positions == 2:
				tmp = remove_third_codon_position(gblocked, outdir)
				os.remove(gblocked)
				gblocked = tmp
			nexus_file = fasta2nexus(gblocked)
			
		concatenated = concatenate(outdir, 'infile')
		charsets = write_charsets(concatenated)
		nexus2phylip(concatenated)
				# Write the initial part of the config file
		write_initial_cfg(outdir)
		
		# Open the partitions file for reading and the config file for appending
		with open(charsets, 'r') as pf, open(f'{outdir}/partition_finder.cfg', 'a') as cfg_file:
				for line in pf:
					match = re.search(r"charset\s+(.*)\s+=\s+([0-9]+)-([0-9]+)", line)
					if match:
						part = match.group(1).split('.')[0]  # Partition name without extension
						from_pos = int(match.group(2))  # First position in the partition
						to_pos = int(match.group(3))		# Last position in the partition

						# Option for three codon positions
						if positions == 3:
							for pos in range(3):
								cfg_file.write(f"{part}_p{pos+1} = {from_pos+pos}-{to_pos}\\3;\n")

						# Option for two codon positions
						elif positions == 2:
							for pos in range(2):
								cfg_file.write(f"{part}_p{pos+1} = {from_pos+pos}-{to_pos}\\2;\n")

		# Write the final part of the config file
		write_final_cfg(outdir)
			

			
	genetic_code = args.code
	outdir = args.outdir
	positions = args.positions
	
	replace_dir(outdir)
	
	tmp = os.path.join(outdir, 'tmp')
	replace_dir(tmp)
	
	genes = args.input
	
	print(pyfiglet.figlet_format("EZpipe")) 
	
	print(f"Running EZpipe with the following parameters:\ngenes path {genes}\ngenetic code: {genetic_code}\npositions: {positions}\noutdir: {outdir}\n\n")
	
	ezpipe_main(genes, genetic_code, tmp, positions)
	shutil.copy(os.path.join(tmp,'infile.phy'), os.path.join(outdir,'infile.phy'))
	shutil.copy(os.path.join(tmp,'partition_finder.cfg'), os.path.join(outdir,'partition_finder.cfg'))
	shutil.rmtree(tmp)
	
	
	

def ez_skew_subcommand(args):
	
	from collections import Counter
	from Bio.SeqUtils import GC123

	def ezskew_main(path, genetic_code, strand, outdir):
		for filename in os.listdir(path):
			filename = os.path.join(path,filename)
			processed_file = check_fasta(filename, outdir)   # Check if the fasta file is valid and free of duplicates
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)  # Check and remove gaps
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			padded_file = pad_sequence(clean_file)
			os.remove(clean_file)
			nexus_file = fasta2nexus(padded_file)
		concatenated = concatenate(outdir, strand)
		return concatenated
		
	def first_skew(fasta):
		ATdic = {}
		for record in SeqIO.parse(fasta, 'fasta'):
			gc_tot, gc_1, gc_2, gc_3 = GC123(record.seq)
			at_tot = 100 - gc_tot
			ATdic[record.id] = at_tot
		first_bias =  pd.DataFrame(ATdic.items(), columns=['Species','AT%'])
		return first_bias
	
	def second_skew(fasta, strand):
		ACdic = {}
		GTdic = {}

		for record in SeqIO.parse(fasta, 'fasta'):
			if strand == 'J':
				dic = dict(Counter(record.seq))
				AC = int(dic.get('A')) + int(dic.get('C'))
				ACperc = AC / len(record.seq) * 100
				ACdic[record.id] = ACperc

			
			if strand == 'N':
				dic = dict(Counter(record.seq))
				GT = int(dic.get('G')) + int(dic.get('T'))
				GTperc = GT / len(record.seq) * 100
				GTdic[record.id] = GTperc

			
			if strand == 'JN':
				dic = dict(Counter(record.seq))
				AC = int(dic.get('A')) + int(dic.get('C'))
				ACperc = AC / len(record.seq) * 100
				ACdic[record.id] = ACperc
				dic = dict(Counter(record.seq))
				GT = int(dic.get('G')) + int(dic.get('T'))
				GTperc = GT / len(record.seq) * 100
				GTdic[record.id] = GTperc
				

		Jbias =  pd.DataFrame(ACdic.items(), columns=['Species','AC%'])
		Nbias =  pd.DataFrame(GTdic.items(), columns=['Species','GT%'])  

		return Jbias, Nbias
		
	def skew_main(string, x, y): #calculates at skew
		i=0 # i is A or G
		j=0 # j is T or C
		for nt in string:
			if nt.upper() == x.upper():
				i+=1
			elif nt.upper() == y.upper():
				j+=1
		if i+j == 0:
			tot_skew = 0 #just in case a+t is zero 
		else:
			tot_skew=(i-j)/(i+j) #at skew
		return tot_skew
		
	def third_skew (fasta, strand): #strand = J or N
		diz={} #empty dic
		for element in SeqIO.parse(fasta, 'fasta'): #element in this case correspond to the taxa present in the dataset
			N=[] #list for the first nucleotide position
			NN=[] #list for the second nucleotide position
			NNN=[] #list for the third nucleotide position
			#get the codons
			if len(element.seq) % 3 == 0:
				codons = [element.seq[i:i+3] for i in range(0, len(element.seq), 3)]
				for codon in codons: #loop to parse single codons
					first=N.append(codon[0])
					second=NN.append(codon[1])
					third=NNN.append(codon[2])
					
							
			
			diz[element.id]=[N,NN,NNN] #this dictionary is completed by the taxa id and the corresponig nucleotides 1st, 2nd and 3rd position
		
		skew = pd.DataFrame.from_dict(diz) #dictionary converted to a dataframe
		skew=skew.transpose()
		skew=skew.reset_index()
		skew.columns=['Species','N','NN','NNN'] #new column names

		#new columns (empty for the moment)
		skew['FIRST_AT_SKEW_'+strand]=''
		skew['FIRST_CG_SKEW_'+strand]=''
		skew['SECOND_AT_SKEW_'+strand]=''
		skew['SECOND_CG_SKEW_'+strand]=''
		skew['THIRD_AT_SKEW_'+strand]=''
		skew['THIRD_CG_SKEW_'+strand]=''

	 	#remove extra characters from the df columns
		skew['NNN'] = skew['NNN'].astype(str).str.replace('[','',regex=False).str.replace(']','',regex=False).str.replace("'",'',regex=False).str.replace(', ','',regex=False)
		skew['N'] = skew['N'].astype(str).str.replace('[','',regex=False).str.replace(']','',regex=False).str.replace("'",'',regex=False).str.replace(', ','',regex=False)
		skew['NN'] = skew['NN'].astype(str).str.replace('[','',regex=False).str.replace(']','',regex=False).str.replace("'",'',regex=False).str.replace(', ','',regex=False)

		#AT-CG skew calculations
		length = len(skew)
		for index in range(0,int(length)): #iterate over index
			row=skew.iloc[index]
			first=row[1] #N col
			second=row[2] #NN col
			third=row[3] #NNN col

			#at and cg skew for each N NN or NNN col and associate this value to the corresponding col (eg third_at_skew col)
			at=skew_main(first, 'A', 'T')
			cg=skew_main(first, 'C', 'G')
			skew.at[index,'FIRST_AT_SKEW_'+strand]=at
			skew.at[index,'FIRST_CG_SKEW_'+strand]=cg
			
			at=skew_main(second, 'A', 'T')
			cg=skew_main(second, 'C', 'G')
			skew.at[index,'SECOND_AT_SKEW_'+strand]=at
			skew.at[index,'SECOND_CG_SKEW_'+strand]=cg
			
			at=skew_main(third, 'A', 'T')
			cg=skew_main(third, 'C', 'G')
			skew.at[index,'THIRD_AT_SKEW_'+strand]=at
			skew.at[index,'THIRD_CG_SKEW_'+strand]=cg
			
			
		skew=skew[['Species', 'FIRST_AT_SKEW_'+strand, 'FIRST_CG_SKEW_'+strand, 'SECOND_AT_SKEW_'+strand, 'SECOND_CG_SKEW_'+strand, 'THIRD_AT_SKEW_'+strand, 'THIRD_CG_SKEW_'+strand]]
		return skew


	def plot_codon_skew(df, x_col, y_col, strand_col, title, ax, legend):
		markers = {"J": "o", "N": "^"}  # Markers for different strands
		species_unique = df['Species'].unique()
		
		# Create a color map for species
		species_colors = plt.get_cmap('tab20', len(species_unique))
		species_color_map = {species: species_colors(i) for i, species in enumerate(species_unique)}
		
		# Plotting each species with its respective color and strand with its marker
		for species in species_unique:
			for strand in markers.keys():
				subset = df[(df['Species'] == species) & (df[strand_col] == strand)]
				ax.scatter(subset[x_col], subset[y_col], 
						   label=species, 
						   marker=markers.get(strand), 
						   color=species_color_map[species])


		# Add horizontal and vertical lines at zero
		ax.axhline(0, color='gray', linestyle='--')
		ax.axvline(0, color='gray', linestyle='--')
		
		# Set axis limits based on the max skew values
		ax.set_xlim(-max_first_AT, max_first_AT)
		ax.set_ylim(-max_first_CG, max_first_CG)
		
		# Axis labels and title
		ax.set_xlabel('AT skew')
		ax.set_ylabel('CG skew')
		ax.set_title(title)
		
		if legend:
			# Create a single legend with species names and strand shape as markers
			handles, labels = ax.get_legend_handles_labels()
			by_label = dict(zip(labels, handles))
			ax.legend(by_label.values(), by_label.keys(), title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')

			
		
	genetic_code = args.code
	outdir = args.outdir
	
	replace_dir(outdir)
	
	tmp = os.path.join(outdir, 'tmp')
	replace_dir(tmp)
	
	plots = os.path.join(outdir, 'plots')
	replace_dir(plots)
	
	tables = os.path.join(outdir, 'tables')
	replace_dir(tables)
	
	print(pyfiglet.figlet_format("EZskew")) 
	
	if args.heavy and not args.light:
		J_path = args.heavy
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		print(f"Running EZskew with the following parameters:\nheavy chain fasta files: {J_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and not args.heavy:
		N_path = args.light
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		print(f"Running EZskew with the following parameters:\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and args.heavy:
		J_path = args.heavy
		N_path = args.light
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		JN_tmp = os.path.join(tmp, 'JN')
		replace_dir(JN_tmp)
		print(f"Running EZskew with the following parameters:\nheavy chain fasta files: {J_path}\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	# Ensure at least one of --heavy or --light is provided
	if not args.heavy and not args.light:
		parser.error("At least one of --heavy (-J) or --light (-N) must be provided.")
		
	
	
	if args.heavy and not args.light:
		concatenated = ezskew_main(J_path, genetic_code, 'J', J_tmp)
		final_file = nexus2fasta(concatenated)
		FSK = first_skew(final_file)
		SSK = second_skew(final_file,'J')[0]
		TSK = third_skew(final_file,'J')


	if args.light and not args.heavy:
		concatenated = ezskew_main(N_path, genetic_code, 'N', N_tmp)
		final_file = nexus2fasta(concatenated)
		FSK = first_skew(final_file)
		SSK = second_skew(final_file,'N')[1]
		TSK = third_skew(final_file,'N')

	if args.light and args.heavy:
		
		concatenatedJ = ezskew_main(J_path, genetic_code, 'J', J_tmp)
		concatenatedN = ezskew_main(N_path, genetic_code, 'N', N_tmp)
		
		# copy the files to the JN folder
		src_files = os.listdir(J_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(J_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'J.nexus':
				shutil.copy(full_file_name, JN_tmp)
		src_files = os.listdir(N_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(N_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'N.nexus':
				shutil.copy(full_file_name, JN_tmp)
				
				
		concatenated = concatenate(JN_tmp, 'JN')	
		
		final_file = nexus2fasta(concatenated)		
		final_fileJ = nexus2fasta(concatenatedJ)
		final_fileN = nexus2fasta(concatenatedN)
		
		FSK = first_skew(final_file)
		
		SSKJ = second_skew(final_fileJ,'J')[0]
		SSKN = second_skew(final_fileN,'N')[1]
		
		TSKJ = third_skew(final_fileJ,'J')
		TSKN = third_skew(final_fileN,'N')
		
		
		SSK = SSKJ.merge(SSKN, on = 'Species', how='left')
		TSK = TSKJ.merge(TSKN, on = 'Species', how='left')
		
		
	df_tmp = FSK.merge(SSK, on = 'Species', how='left')
	df_final = df_tmp.merge(TSK, on = 'Species', how='left')
	
	output_table = os.path.join(tables,'Final_table.csv')
	df_final.to_csv(output_table)
	
	# Taking the relevant columns based on gene
	if args.light and args.heavy:
		AT = df_final[['Species', 'AT%']]
		AC = df_final[['Species', 'AC%']]
		GT = df_final[['Species', 'GT%']]
		AT.columns = ['Species', 'value']
		AC.columns = ['Species', 'value']
		GT.columns = ['Species', 'value']
		AT.loc[:, 'strand'] = 'AT%'
		AC.loc[:, 'strand'] = 'AC% (J/heavy strand)'
		GT.loc[:, 'strand'] = 'GT% (N/light strand)'
		df_final2 = pd.concat([AT, AC, GT])
	if args.heavy and not args.light:
		AT = df_final[['Species', 'AT%']]
		AC = df_final[['Species', 'AC%']]
		AT.columns = ['Species', 'value']
		AC.columns = ['Species', 'value']
		AT.loc[:, 'strand'] = 'AT%'
		AC.loc[:, 'strand'] = 'AC% (J/heavy strand)'
		df_final2 = pd.concat([AT, AC])
	if args.light and not args.heavy:
		AT = df_final[['Species', 'AT%']]
		GT = df_final[['Species', 'GT%']]
		AT.columns = ['Species', 'value']
		GT.columns = ['Species', 'value']
		AT.loc[:, 'strand'] = 'AT%'
		GT.loc[:, 'strand'] = 'GT% (N/light strand)'
		df_final2 = pd.concat([AT, GT])
		
	

	# Plotting freqpoly plot
	if len(df_final) > 30:
		fig, ax = plt.subplots(figsize=(10, 6))

		strand_colors = {'AT%': '#e41a1c', 'AC% (J/heavy strand)': '#377eb8', 'GT (N/light strand)': '#4daf4a'}

		for s in df_final2['strand'].unique():
			subset = df_final2[df_final2['strand'] == s]
			ax.hist(subset['value'], histtype='step', 
					stacked=True, fill=True, bins=10, label=s, color=strand_colors.get(s))

		ax.set_xlabel('%')
		ax.set_ylabel('Frequency')
		ax.legend(title='Bias percentages')

		if args.light and args.heavy == 'JN':
			ax.set_title("First and second bias frequencies: AT%, AC% and GT%")
		elif args.light == 'N':
			ax.set_title("First and second bias frequencies: AT% and GT%")
		elif args.heavy == 'J':
			ax.set_title("First and second bias frequencies: AT% and AC%")

		output_file = os.path.join(plots, 'First_and_second_bias_frequency.pdf')
		plt.savefig(output_file, dpi=500, bbox_inches='tight')
		 
		 
	#plotting barplots	 
	fig, ax = plt.subplots(figsize=(14, 8))

	# Get unique species and strands
	species = df_final2['Species'].unique()
	strands = df_final2['strand'].unique()
	n_strands = len(strands)

	# Set the width of each bar and create offsets
	bar_width = 0.2
	x = np.arange(len(species))

	strand_colors = {'AT%': '#e41a1c', 'AC% (J/heavy strand)': '#377eb8', 'GT (N/light strand)': '#4daf4a'}


	# Plot bars for each strand with an offset to dodge them
	for i, s in enumerate(strands):
		subset = df_final2[df_final2['strand'] == s]
		ax.bar(x + i * bar_width, subset['value'], width=bar_width, label=s, color=strand_colors.get(s))

	# Adjust the x-ticks and labels
	ax.set_xticks(x + bar_width * (n_strands - 1) / 2)
	ax.set_xticklabels(species, rotation=90)

	# Set labels and title based on the strand
	if 'JN' in strands:
		ax.set_title("First and second bias: AT%, AC% and GT%")
	elif 'N' in strands:
		ax.set_title("First and second bias: AT% and GT%")
	elif 'J' in strands:
		ax.set_title("First and second bias: AT% and AC%")

	# Set axis labels
	ax.set_ylabel('%')
	ax.set_xlabel('Species')

	# Add legend
	ax.legend(labels=strands)
	
	output_file = os.path.join(plots, 'First_and_second_bias.pdf')
	plt.savefig(output_file, dpi=500, bbox_inches='tight')
	
	# Plotting the skews
	if args.light and args.heavy:
		J = df_final[['Species','FIRST_AT_SKEW_J', 'FIRST_CG_SKEW_J', 'SECOND_AT_SKEW_J', 'SECOND_CG_SKEW_J', 'THIRD_AT_SKEW_J', 'THIRD_CG_SKEW_J']]
		N = df_final[['Species','FIRST_AT_SKEW_N', 'FIRST_CG_SKEW_N', 'SECOND_AT_SKEW_N', 'SECOND_CG_SKEW_N', 'THIRD_AT_SKEW_N', 'THIRD_CG_SKEW_N']]
		J.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		N.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		J['strand'] = 'J'
		N['strand'] = 'N'
		JN = pd.concat([J, N])
	if args.heavy and not args.light:
		JN = df_final[['Species','FIRST_AT_SKEW_J', 'FIRST_CG_SKEW_J', 'SECOND_AT_SKEW_J', 'SECOND_CG_SKEW_J', 'THIRD_AT_SKEW_J', 'THIRD_CG_SKEW_J']]
		JN.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		JN['strand'] = 'J'
	if args.light and not args.heavy:
		JN = df_final[['Species','FIRST_AT_SKEW_N', 'FIRST_CG_SKEW_N', 'SECOND_AT_SKEW_N', 'SECOND_CG_SKEW_N', 'THIRD_AT_SKEW_N', 'THIRD_CG_SKEW_N']]
		JN.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		JN['strand'] = 'N'

	# Dividing data for codon positions
	first = JN[['Species', 'FAT', 'FCG', 'strand']]
	second = JN[['Species', 'SAT', 'SCG', 'strand']]
	third = JN[['Species', 'TAT', 'TCG', 'strand']]

	# Calculate max values for axis limits
	max_first_AT = max(abs(first['FAT']).max(), abs(second['SAT']).max(), abs(third['TAT']).max())
	max_first_CG = max(abs(first['FCG']).max(), abs(second['SCG']).max(), abs(third['TCG']).max())

	# Unique species list for color assignment
	species_unique = first['Species'].unique()
	species_colors = plt.get_cmap('tab20', len(species_unique))

	# Map species to specific colors
	species_color_map = {species: species_colors(i) for i, species in enumerate(species_unique)}
	
	fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
		
	plot_codon_skew(first, 'FAT', 'FCG', 'strand', 'First codon position', axes[0], legend = False)
	plot_codon_skew(second, 'SAT', 'SCG', 'strand', 'Second codon position', axes[1], legend = False)
	plot_codon_skew(third, 'TAT', 'TCG', 'strand', 'Third codon position', axes[2], legend = True)

	# Show plots
	plt.tight_layout()

	output_file = os.path.join(plots, 'Third_bias.pdf')
	plt.savefig(output_file, dpi=500)
			
	shutil.rmtree(tmp)

def ez_split_subcommand(args):
	
	from BCBio.GFF import GFFExaminer
	from BCBio import GFF

	fasta_file = args.input
	gff_file = args.gff
	outdir = args.outdir

	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZsplit")) 
	
	print(f"Running EZsplit with the following parameters:\nfasta file: {fasta_file}\nGFF file: {gff_file}\noutdir: {outdir}\n\n")
	
	is_gff3(gff_file)
	is_fasta(fasta_file)
	
	examiner = GFFExaminer()
	in_handle = open(gff_file)

	ids = [list(i)[0] for i in examiner.available_limits(in_handle).get('gff_id').keys()]

	names={	'ATP6': ['ATP6', 'A6','MT-ATP6', 'ATPASE6'] ,
		'ATP8': ['ATP8','A8', 'MT-ATP8', 'ATPASE8'] ,
		'COX1': ['COX1','MT-CO1', 'CO1',  'COXI', 'MTCO1', 'MTCOX1','MTCOXI'] ,	   
		'COX2': ['COX2','MT-CO2', 'CO2',  'COXII', 'MTCO2', 'MTCOX2','MTCOXII'] ,
		'COX3': ['COX3','MT-CO3', 'CO3',  'COXIII', 'MTCO3', 'MTCOX3','MTCOXIII'] ,
		'COB': ['COB','MT-CYB', 'CYTB', 'MTCYB'] ,
		'ND1': ['NAD1', 'MT-ND1', 'MTND1', 'NADH1', 'ND1'] ,
		'ND2': ['NAD2', 'MT-ND2', 'MTND2', 'NADH2', 'ND2' ] ,
		'ND3': ['NAD3', 'MT-ND3', 'MTND3', 'NADH3', 'ND3'] ,
		'ND4': ['NAD4', 'MT-ND4', 'MTND4', 'NADH4', 'ND4'] ,
		'NDL': ['NADL','MT-ND4L', 'MTND4L', 'NADH4L', 'ND4L', 'NDL', 'NAD4L', 'NADHL'] ,
		'ND5': ['NAD5', 'MT-ND5', 'MTND5', 'NADH5', 'ND5'] ,
		'ND6': ['NAD6', 'MT-ND6', 'MTND6', 'NADH6', 'ND6']
		}
		
	list_of_genes = list(names.keys())
	d_found_genes = {}
	for rec in GFF.parse(gff_file):
		rec_id = rec.id
		genes_found = []
		for feature in rec.features:
			try:
				if feature.qualifiers.get('gene_biotype')[0] == 'protein_coding':
					start, end, strand = int(str(feature.location.start)), int(str(feature.location.end)), str(feature.location.strand)
					gene = feature.qualifiers.get('Name')[0].upper()
					for k, v in names.items():
						if gene in v:
							gene_syn = k
						else:
							pass
					for record in SeqIO.parse(fasta_file, 'fasta'):
						if rec_id == record.id:
							record.description = record.description
							record.seq = record.seq[start:end]
							if strand == '-1':
								record.seq = record.seq.reverse_complement()

							genes_found.append(gene_syn)
							with open(f'{outdir}/{gene_syn.lower()}.fasta', 'a') as outfile:
								outfile.write('>'+str(record.description)+'\n'+str(record.seq)+'\n')

			except:
				pass
			d_found_genes[rec_id] = genes_found

	with open(f'{outdir}/missing_genes.txt', 'w') as outfile:
		for k,v in d_found_genes.items():
			missing_genes = list(set(list_of_genes) - set(v))
			if len(missing_genes) > 0:
				outfile.write(f'{k}\t{missing_genes}\n')
	



def main():
	parser = argparse.ArgumentParser(
		description="Choose the EZmito program you want to use.",
		formatter_class=argparse.RawTextHelpFormatter
	)
	
	subparsers = parser.add_subparsers(
		title="Available commands",
#		description="Subcommands for different analyses.",
		help="Use one of these subcommands for the corresponding analysis.",
		dest="command"
	)

	# EZcircular subcommand
	parser_ez_circular = subparsers.add_parser(
		'ezcircular', 
		help='Circularize a sequence from another starting gene',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_circular.add_argument("-f", "--feature", help="Genome feature. [circular or linear]", default='circular', type=str)
	parser_ez_circular.add_argument("-s", "--start", help="Starting gene", default='cox1', type=str)
	parser_ez_circular.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_circular = parser_ez_circular.add_argument_group('required named arguments')
	required_ez_circular.add_argument("-i", "--input", required=True, help="FASTA input file")
	required_ez_circular.add_argument("-b", "--bed", required=True, help="BED input file")
	parser_ez_circular.set_defaults(func=ez_circular_subcommand)

	# EZcodon subcommand
	parser_ez_codon = subparsers.add_parser(
		'ezcodon', 
		help='Analyze codon usage of J (heavy) and N (light) chains',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_codon.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_codon = parser_ez_codon.add_argument_group('required named arguments')
	required_ez_codon.add_argument("-J", "--heavy", help="Path to J chain FASTA files")
	required_ez_codon.add_argument("-N", "--light", help="Path to N chain FASTA files")
	required_ez_codon.add_argument("-c", "--code", required=True, help="Genetic code as an integer", type=int)

	# Set the function to be called by the ezcodon subcommand
	parser_ez_codon.set_defaults(func=ez_codon_subcommand)


	
	# EZmap subcommand
	parser_ez_map = subparsers.add_parser(
		'ezmap', 
		help='Create a custom map plot of your mitogenome',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_map.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_map.add_argument("-f", "--feature", help="Genome feature. [circular or linear]", default='circular', type=str)
	parser_ez_map.add_argument("-colJ", "--colorJ", help="J (heavy) strand color", default='#add8e6', type=str)
	parser_ez_map.add_argument("-colN", "--colorN", help="N (light) strand color", default='#B22222', type=str)
	required_ez_map = parser_ez_map.add_argument_group('required named arguments')
	required_ez_map.add_argument("-g", "--gff", required=True, help="GFF3 input file")
	parser_ez_map.set_defaults(func=ez_map_subcommand)
	
	# EZmix subcommand
	parser_ez_mix = subparsers.add_parser(
		'ezmix', 
		help='Check for possible chimeras occurred during your mitogenome assemblies',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_mix.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_mix.add_argument("-bn", "--blastn", help="BLASTn path", default='', type=str)
	parser_ez_mix.add_argument("-id", "--identity", help="Identity threshold [0.5-1]", default=0.95, type=float)
	parser_ez_mix.add_argument("-len", "--length", help="Length threshold (bp)", default=200, type=int)
	required_ez_mix = parser_ez_mix.add_argument_group('required named arguments')
	required_ez_mix.add_argument("-i", "--input", required=True, help="MultiFASTA input file")
	parser_ez_mix.set_defaults(func=ez_mix_subcommand)

	# EZpipe subcommand
	parser_ez_pipe = subparsers.add_parser(
		'ezpipe', 
		help='Prepare J and N (heavy and light) FASTA files for phylogenetic analysis',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_pipe.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_pipe.add_argument("-p", "--positions", help="Number of codon positions to analyze (should the program remove 3rd codon positions?) [2,3]", default=3, type=int)
	required_ez_pipe = parser_ez_pipe.add_argument_group('required named arguments')
	required_ez_pipe.add_argument("-i", "--input", required = True, help="Path to FASTA files")
	required_ez_pipe.add_argument("-c", "--code", required=True, help=code_help, type=int)
	parser_ez_pipe.set_defaults(func=ez_pipe_subcommand)


	# EZskew subcommand
	parser_ez_skew = subparsers.add_parser(
		'ezskew', 
		help='Analyze codon skew biases of J (heavy) and N (light) chains',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_skew.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_skew = parser_ez_skew.add_argument_group('required named arguments')
	required_ez_skew.add_argument("-J", "--heavy", help="Path to J chain FASTA files")
	required_ez_skew.add_argument("-N", "--light", help="Path to N chain FASTA files")
	required_ez_skew.add_argument("-c", "--code", required=True, help="Genetic code as an integer", type=int)

	# Set the function to be called by the ezskew subcommand
	parser_ez_skew.set_defaults(func=ez_skew_subcommand)
	
	

	# EZsplit subcommand
	parser_ez_split = subparsers.add_parser(
		'ezsplit', 
		help='Split multiple Genbank downloaded multifasta to a set of PCGs',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_split.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_split = parser_ez_split.add_argument_group('required named arguments')
	required_ez_split.add_argument("-g", "--gff", required=True, help="GFF3 input file")
	required_ez_split.add_argument("-i", "--input", required=True, help="MultiFASTA input file of complete mitogenomes")
	parser_ez_split.set_defaults(func=ez_split_subcommand)

	# Parse the arguments and call the appropriate subcommand function
	args = parser.parse_args()
	
	if args.command is None:
		parser.print_help()
	else:
		args.func(args)

if __name__ == "__main__":
	main()
	# End of the script: measure runtime
	end_time = time.time()
	runtime = round(end_time - start_time, 2)  # Runtime in seconds

	# Print runtime and citation information to the screen
	print(f"\n\n-------------------The process was correctly completed in {runtime} seconds-------------------")
	print("Thank you for using this code. If it helped you, please cite:")
	print("Cucini C., Leo C., Iannotti N., Boschi S., Brunetti C., Pons J., Fanciulli P. P., Frati F., Carapelli A., & Nardi F. (2021)")
	print("EZmito: a simple and fast tool for multiple mitogenome analyses, Mitochondrial DNA Part B, 6(3), 1101-1109.")
	print("Doi: 10.1080/23802359.2021.1899865")

