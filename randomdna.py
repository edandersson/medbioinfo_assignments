#!/usr/bin/env python2.7
## This script generates a random DNA sequence of length specified by the user

import random

# User input of sequence length to generate
SeqLength = int(raw_input("Length: "))

# Define list with elements to build the sequence and their proportions in the generated sequence 
Elements = ['A','T','C','G']
Proportions = [1,1,1,1]
Sample = ''
# Generate string of elements with desired distribution of elements
x = 0
for i in Elements:
	y = Proportions[x]
	Sample= Sample + (i * y)
	x = x+1

# Convert the string into a list to be used by the choice function
Sample = list(Sample)


# Build the randomly pick elements from the generated distribution
RandomSequence = ''

for i in range(SeqLength):
	Random = random.choice(Sample)
	RandomSequence = RandomSequence + Random

# Print sequence in FASTA format
print ">Random.sequence.%d\n%s" % (SeqLength,RandomSequence)
