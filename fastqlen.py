#!/usr/bin/env python2.7

Usage= """
fastqlen.py - version 1
This script reads FastQ files and reports the number of sequences and their identifiers.

Usage: fastaqlen.py FILE"""

import os
import sys

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
		
if init() == 0:
	quit()
	
ImportFile = open(sys.argv[1], 'r')

LineCounter = 1
Sequences = []
InSeq = 0

for Line in ImportFile:
	if Line[0] == '@' and InSeq == 0:	# Line starts with @ and not in lines following previous sequence
		# Strip name from unnecessary signs
		Name = Line[1:]
		Name = Name.strip('\n')
		Sequences.append(Name)
		InSeq = 1		# Set in sequence marker in case quality line commence with @
		LineCounter = 1 # Reset counter for new sequence
	else:
		if InSeq == 1:
			LineCounter = LineCounter + 1 # Count lines
	# Leave sequence at 4th line
	if LineCounter == 4: 
		InSeq = 0

ImportFile.close()

print len(Sequences)
for i in range(len(Sequences)):
	print Sequences[i]
