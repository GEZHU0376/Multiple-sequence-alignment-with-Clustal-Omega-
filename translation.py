#!/usr/bin/env python2
# script.py
#This program translate nucleotides to amino acid
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def CodonsProtein( input_value, gencode):
    protein = ''
    
    for i in range(0, len(input_value), 3):
        codon = str(input_value)[i:i+3]
        if (codon in gencode.keys()):
            protein += gencode[codon]
    return protein

input_file = 'APOE_refseq_transcript.fasta'
ouput_file = 'APOE_refseq_transcript.output.fasta'

# Bio.Seq convert fasta file into a dictionary
seq_dict = SeqIO.to_dict(SeqIO.parse("APOE_refseq_transcript.fasta","fasta" ))

translated_ouput = open(ouput_file,"w")

for element in seq_dict.keys():
    codon = seq_dict[element]
    protein = CodonsProtein(codon.seq, gencode)
    header = ">" + element + "\n"
    translated_ouput.write(header)

    protein_sequence = protein + "\n"
    translated_ouput.write(protein_sequence)

translated_ouput.close()
