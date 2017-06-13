#!/usr/bin/env python2.7

Usage= """
dna2aa.py - version 1
This script finds the longest open reading frame and saves a file with the result.

Usage: dna2aa.py FILE"""

import os
import sys
import operator


# Usage and input check
def init():
	if len(sys.argv) < 2:
		print Usage
		return(0)
	else: 
		ImportFilePath = sys.argv[1]

		if os.path.isfile(ImportFilePath) == 0:
			sys.stderr.write("File does not exist!\n")
			return(0)
		else:
			return(1)

# Function for producing the compliment strand of the DNA sequence
def compliment_string(String):
	# Define disctionary of complimentary bases
	ComplimentBase = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
	
	# Loop sequence backwards and build string with bases from dictionary
	ReverseComplimentString = ''
	for i in range(len(String)):
		Complimentary = ComplimentBase.get(String[-(i+1)],'X')
		ReverseComplimentString += Complimentary
	
	return(ReverseComplimentString)

# Function returning a aminoacid sequence from a DNA sequence(first parameter)
# with reading frame set to: 1st, 2nd or 3rd base(2nd parameter)
def translation_dna(String, Start=1):
	CodonDict = {'ATT':'I', 'ATC':'I', 'ATA':'I',
	'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 'TTG':'L',
	'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
	'TTT':'F', 'TTC':'F',
	'ATG':'M',
	'TGT':'C', 'TGC':'C',
	'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
	'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
	'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
	'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
	'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
	'TAT':'Y', 'TAC':'Y',
	'TGG':'W',
	'CAA':'Q', 'CAG':'Q',
	'AAT':'N','AAC':'N',
	'CAT':'H', 'CAC':'H',
	'GAA':'E', 'GAG':'E',
	'GAT':'D', 'GAC':'D',
	'AAA':'K', 'AAG':'K',
	'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
	'TAA':'*', 'TAG':'*', 'TGA':'*'}
	
	# Shift reading frame to 2nd or 3rd base
	if Start == 2:
		String = String[1:]
	if Start == 3:
		String = String[2:]
	
	Peptide = ''
	
	# Loop through the string 3 bases at the time
	for Codon in range(0,len(String),3): 
		Peptide += CodonDict.get(String[Codon:(Codon+3)],'X')
	return Peptide

# Function returning length of open reading frame from start acid(2nd param.) to stop codon.
def orf_length(AminoAcids, Start):
	Length_ORF = 0
	for i in range(Start,len(AminoAcids)):
		if AminoAcids[i] == '*' or AminoAcids[i] == 'X':
			return Length_ORF
		else:
			Length_ORF += 1
		return Length_ORF

# Function returning the longest open reading frame in a given aminoacid sequence
def longest_orf(AminoAcids):
	AminoAcids = AminoAcids.replace('X','*')
	ORF = AminoAcids.split('*')
	Items = []
	for Seq in ORF:
		Items.append([len(Seq),Seq])
	LongestORF = max(Items, key=operator.itemgetter(1))
	return LongestORF[1]
	
# Function returning open reading frame of a sequence based on 
# reading frame and number of first and last amino acid
def open_reading_frame(String,RF,Start,Stop):
	# Cut bases depending on reading frame
	Sequence = String[(RF-1):]
	# Cut bases until start codon
	Sequence = Sequence[(Start*3):]
	# Keep bases until stop
	Sequence = Sequence[:(Stop*3)]
	return Sequence

# Function reading file and parsing FASTA sequences into a list [name][seq]
def read_sequence():
	ImportFileName = sys.argv[1]
	ImportFile = open(ImportFileName, 'r')
	
	Sequences = []
	Counter = -1
	
	for Line in ImportFile:
		Line = Line.strip()
		if Line[0] == '>':
			Name = Line[1:]
			Sequences.append([Name,''])
			Counter += 1
		else:
			if Counter > -1:
				Sequences[Counter][1] += Line.upper()
	
	return(Sequences)
	
# Function for printing a FASTA formated file
def write_fasta_file(SequenceList):
	OutputFileName = "%s_longest_orf.fasta" % sys.argv[1]
	OutputFile = open(OutputFileName, 'w')
	for Seq in SequenceList:
		# Split long sequences over several lines
		Lines = [Seq[1][i:70+i] for i in xrange(0, len(Seq[1]), 70)] # Bases per line set to 70
		
		# Write header
		OutputFile.write('>'+Seq[0]+'\n')
		# Loop through sequences lines
		for Line in Lines:
			OutputFile.write(Line+'\n')
		
	OutputFile.close()

## Main
if init() == 0:
	quit()

# Read DNA sequences
DNASequences = read_sequence()

# Define list for resulting longest open reading frame (ORF)
FinalORF = [] # [Sequence name, sequence]

# Step through imported sequences
for Seq in DNASequences:
	OriginalStrand = Seq[1]
	ComplimentaryStrand = compliment_string(OriginalStrand)
	
	# Loop through possible reading frames
	IntraStrandVariants = [] # [AA Sequence, reading frame, forward/reverse]
	for i in range(1,4):
		IntraStrandVariants.append([translation_dna(OriginalStrand, i),i,'forward'])
		IntraStrandVariants.append([translation_dna(ComplimentaryStrand, i),i,'reverse'])
	
	# Loop through all aminoacid sequences
	IntraStrandORF = [] # [AA Sequence, reading frame, forward/reverse, StartAminoAcid, Length of ORF]
	for Strand in IntraStrandVariants:
		#Loop through the whole length of the sequence
		for Start in range(len(Strand[0])):
			LenORF = orf_length(Strand[0], Start)
			IntraStrandORF.append([Strand[0],Strand[1],Strand[2], Start, LenORF])
	
	# Determine longest ORF
	MaxORF = max(IntraStrandORF, key=operator.itemgetter(4))
 	
 	# Extract DNA sequence of longest ORF
 	if MaxORF[2] == 'forward':
 		FinalORF.append([Seq[0],open_reading_frame(OriginalStrand, MaxORF[1],MaxORF[3],MaxORF[4])])
 	else: # Eg. the reverse strand
 		FinalORF.append([Seq[0],open_reading_frame(ComplimentaryStrand, MaxORF[1],MaxORF[3],MaxORF[4])])

# Write output to FASTA formated file

write_fasta_file(FinalORF)

print "Aminoacid sequence(s):"
for ORF in FinalORF:
	print '>'+ORF[0]
	print translation_dna(ORF[1], 1)
