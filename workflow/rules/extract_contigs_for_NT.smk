#----- Extract the Probable Phage Contig Sequences to a New File to Query against NT with MMSeqs2 -----#

# extract_contigs_for_NT.smk

rule extract_contigs_for_NT:
	input:
		identifiers = os.path.join("results", "phage_contigs_for_NT", "phage_contigs_for_NT.txt"),
		sequences = os.path.join("results", "contig_dictionary_filtered", "assembly_1kb.fasta")
	output:
		os.path.join("results", "phage_contigs_for_NT", "phage_contigs_for_NT.fasta")
	run:
		# open output file to hold sequences, and the input sequence file to read from
		output_file = open(output[0], "a")
		contig_fasta_file = open(input[1], "r")

		# initiate empty contig dictionary to hold the contig IDs
		contig_dictionary = {}

		# open the contig ID file
		contig_ID_file = open(input[0], "r")

		# Put contig IDs into the contig dictionary
		for line in contig_ID_file:
			contig_dictionary[line.rstrip("\r\n")] = True
			copyNextLine = False

		# Process the contig fasta file line-by-line to copy the appropriate contig sequences
		for line in contig_fasta_file:
			if ">" in line:
				copyNextLine = False
				columns = line.split()
				contig_ID = columns[0].replace('>', '').rstrip("\n")

				if contig_ID in contig_dictionary.keys():
					output_file.write(line)
					copyNextLine = True
			
			else:
				if copyNextLine == True:
					output_file.write(line)

		output_file.close()
		contig_fasta_file.close()
		contig_ID_file.close()

		
