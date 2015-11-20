# Gregory McInnes
# gmcinnes@stanford.edu
# November 20, 2015

# compound-hets.py
# Identify compound heterozygous variants from a vcf-like file.  This script will look over each gene
# in the list and determine if there are any pairs of variants for each sample that were not present in either
# the mother or father.

# Questions
# Is one frequency required for all datasets or one for each?
# Decide how to handle variants where there is no reported frequency for a provided metric
# Some variants have multiple reported consequences.  How to handle these cases?
# Should frequency cutoffs be minimums or maximums?
# What is the ped file for?

import argparse
from collections import defaultdict

# File index
CHROMOSOME = 0
POSITION = 1
REF = 2
ALT = 3
INFO = 4
DENOVO = 5
DENOVO_COUNT_AFFECTED = 6
DENOVO_IDS_AFFECTED = 7
DENOVO_COUNT_UNAFFECTED = 8
DENOVO_IDS_UNAFFECTED = 9
NEWLY_HOMOZYGOUS = 10
NEWLY_HOMOZYGOUS_COUNT_AFFECTED = 11
NEWLY_HOMOZYGOUS_IDS_AFFECTED = 12
NEWLY_HOMOZYGOUS_COUNT_UNAFFECTED = 13
NEWLY_HOMOZYGOUS_IDS_UNAFFECTED = 14
HEMIZYGOUS = 15
HEMIZYGOUS_COUNT_AFFECTED = 16
HEMIZYGOUS_IDS_AFFECTED = 17
HEMIZYGOUS_COUNT_UNAFFECTED = 18
HEMIZYGOUS_IDS_UNAFFECTED = 19
FROM_MOTHER = 20
FROM_MOTHER_COUNT_AFFECTED = 21
FROM_MOTHER_IDS_AFFECTED = 22
FROM_MOTHER_COUNT_UNAFFECTED = 23
FROM_MOTHER_IDS_UNAFFECTED = 24
FROM_FATHER = 25
FROM_FATHER_COUNT_AFFECTED = 26
FROM_FATHER_IDS_AFFECTED = 27
FROM_FATHER_COUNT_UNAFFECTED = 28
FROM_FATHER_IDS_UNAFFECTED = 29
UNKNOWN_PHASE = 30
UNKNOWN_PHASE_COUNT_AFFECTED = 31
UNKNOWN_PHASE_IDS_AFFECTED = 32
UNKNOWN_PHASE_COUNT_UNAFFECTED = 33
UNKNOWN_PHASE_IDS_UNAFFECTED = 34
MISSING = 35
MISSING_COUNT_AFFECTED = 36
MISSING_IDS_AFFECTED = 37
MISSING_COUNT_UNAFFECTED = 38
MISSING_IDS_UNAFFECTED = 39
UNCERTAIN = 40
UNCERTAIN_COUNT_AFFECTED = 41
UNCERTAIN_IDS_AFFECTED = 42
UNCERTAIN_COUNT_UNAFFECTED = 43
UNCERTAIN_IDS_UNAFFECTED = 44

# Class to parse file and identify compound hets for each sample.
# Input:
#   input: input file
#   ped: ped file
#   thousand_freq: maximum thousand genomes frequency allowable
#   exac_freq: maximum ExAC frequency allowable
#   esp_freq: maximum esp frequency allowable
#   complete_freq: maximum complete genomics frequency allowable
#   consequence: required predicted variant consequence
#   output: output file prefix
class CompoundHets(object):
    def __init__(self, input, ped=None, thousand_freq=None, exac_freq=None, esp_freq=None, complete_freq=None,
                 consequence=None, output=None):
        self.thousand_cutoff = thousand_freq
        self.exac_cutoff = exac_freq
        self.esp_cutoff = esp_freq
        self.complete_cutoff = complete_freq
        self.consequence = consequence
        genes = self.read_input_file(input)
        samples, positions = self.check_genes(genes)
        self.print_positions(positions, output)
        self.print_samples(samples, output)

    # Print each compound het pair for each sample to a sample file
    # Input:
    #   samples: dictionary of lists of tuples of Variants
    #   output: output prefix
    def print_samples(self, samples, output):
        if not samples:
            return
        for s in samples:
            filename = "%s.compound_het_pairs.txt" % s
            if output:
                filename = "%s.%s.compound_het_pairs.txt" % (output, s)
            f = open(filename, 'w')
            f.write("#PAIR\tCHR\tPOS\tGENE\n")
            count = 1
            for pair in samples[s]:
                v1 = pair[0]
                v2 = pair[1]
                f.write(v1.for_sample_file(count) + "\n")
                f.write(v2.for_sample_file(count) + "\n")
                count += 1
            f.close()

    # Print all pairs of compound hets to file
    # Input:
    #   positions: A list of tuples of variants
    #   output: output file prefix
    def print_positions(self, positions, output):
        if not positions:
            return
        filename = "all_compound_het_pairs.txt"
        if output:
            filename = "%s.all_compound_het_pairs.txt" % output
        f = open(filename, 'w')
        f.write("#PAIR\tCHR\tPOS\tGENE\tFROM_MOTHER\tFROM_FATHER\n")
        count = 1
        for p in positions:
            v1 = p[0]
            v2 = p[1]
            f.write(v1.for_variant_file(count) + "\n")
            f.write(v2.for_variant_file(count) + "\n")
            count += 1
        f.close()

    # For each gene in the gene list, check whether each sample has a compound het pair within the gene
    # Input:
    #   genes: dictionary of lists of Variants
    # Return:
    #   samples: dictionary of lists of tuples of Variants
    #   positions: A list of tuples of variants
    def check_genes(self, genes):
        samples = defaultdict(list)
        positions = []
        for g in genes:
            if len(genes[g]) == 1:
                continue
            for v1 in genes[g]:
                for paternal_id in v1.from_father:
                    for v2 in genes[g]:
                        if v2.has_maternal(paternal_id):
                            if (v1, v2) not in positions and (v2, v1) not in positions:
                                positions.append((v1, v2))
                            samples[paternal_id].append((v1, v2))
        return samples, positions

    # Parse each line in the input file, convert it to a Variant object, and assign it to a gene within the gene
    #  dictionary
    # Input:
    #   input: input filename
    # Return:
    #   genes: dictionary of lists of Variants
    def read_input_file(self, input):
        genes = defaultdict(list)
        with open(input) as f:
            for line in f:
                variant = self.parse_line(line)
                if variant is not None:
                    genes[variant.gene].append(variant)
        return genes

    # Parse a line to determine chromosome, position, and other relevant information
    # Input:
    #   line: A line directly from the input file
    # Return:
    #   variant: Variant object. None if it does not meet the input criteria
    def parse_line(self, line):
        fields = line.split("\t")
        chr = fields[CHROMOSOME]
        pos = fields[POSITION]
        info = fields[INFO]
        from_mother = fields[FROM_MOTHER_IDS_AFFECTED]
        from_father = fields[FROM_FATHER_IDS_AFFECTED]
        variant = Variant(chr, pos, info, from_mother, from_father)
        if variant.pass_1kg(self.thousand_cutoff) and variant.pass_complete(self.complete_cutoff) and \
            variant.pass_esp(self.esp_cutoff) and variant.pass_exac(self.exac_cutoff) and \
            variant.pass_consequence(self.consequence) and variant.has_ids():
            return variant
        return None

# Variant Class
#  This class creates and object that stores information for each variant.  Upon initialization, the object
#  fills in information regarding which gene it is found within, various frequencies if availalble, and variant
#  consequence.
# Input:
#   chr: Chromosome
#   pos: Position
#   info: Info field from VCF
#   from_mother: IDs of the samples that inherited this variant from the mother
#   from_father: IDs of the samples that inherited this variant from the father
#   debug: Print debugging information
class Variant(object):
    def __init__(self, chr, pos, info, from_mother, from_father, debug=False):
        self.chr = chr
        self.pos = pos
        self.info = info
        self.from_mother = self.parse_ids(from_mother)
        self.from_father = self.parse_ids(from_father)
        self.debug = debug
        fields = self.info.split(";")
        self.k1g_freq = self.get_1kg_freq(fields)
        self.esp_freq = self.get_esp_freq(fields)
        self.exac_freq = self.get_exac_freq(fields)
        self.complete_freq = self.get_complete_freq(fields)
        self.gene, self.consequence = self.get_gene_and_consequence(fields)

    # Print the variant to the console.  Debugging method.
    def print_variant(self):
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chr, self.pos, self.gene, ",".join(self.consequence),
                                                          self.k1g_freq, self.esp_freq,
                                                          self.exac_freq, self.complete_freq,
                                                          ",".join(self.from_mother), ",".join(self.from_father))

    # Prepare a string for printing to the sample file
    # Input:
    #   pair: The pair id for the variant
    # Return:
    #   string: a string of information for printing
    def for_sample_file(self, pair):
        string = "%s\t%s\t%s\t%s" % (pair, self.chr, self.pos, self.gene)
        return string

    # Prepare a string for printing to the variant file
    # Input:
    #   pair: The pair id for the variant
    # Return:
    #   string: a string of information for printing
    def for_variant_file(self, pair):
        string = "%s\t%s\t%s\t%s\t%s\t%s" % (pair, self.chr, self.pos, self.gene, ",".join(self.from_mother),
                                          ",".join(self.from_father))
        return string

    # Get the gene symbol and consequence from the info field
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   gene: gene symbol
    #   consequence: variant consequence
    def get_gene_and_consequence(self, fields):
        for f in fields:
            if f.startswith("CSQ"):
                info = f.split("|")
                consequence = info[1].split('&')
                gene = info[3]
                return gene, consequence
        if self.debug:
            print "Gene and consequence could not be identified!"

    # Parse the ids from either the mother or the father
    # Input:
    #   id_string: string of ids separated by commas
    # Return:
    #   ids: list of ids
    def parse_ids(self, id_string):
        if id_string == "NA":
            return []
        ids = id_string.split(",")
        return ids

    # Get the thousand genome allele frequency from the VCF info
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   1000 genomes allele frequency
    def get_1kg_freq(self, fields):
        for f in fields:
            if f.startswith("1000g2014oct_all"):
                return f.split("=")[1]
        if self.debug:
            print "1,000 Genomes allele frequency not found!"

    # Get the esp allele frequency from the VCF info
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   esp allele frequency
    def get_esp_freq(self, fields):
        for f in fields:
            if f.startswith("esp6500siv2_all"):
                return f.split("=")[1]
        if self.debug:
            print "esp allele frequency not found!"

    # Get the ExAC allele frequency from the VCF info
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   ExAC allele frequency
    def get_exac_freq(self, fields):
        for f in fields:
            if f.startswith("ExAC_ALL"):
                return f.split("=")[1]
        if self.debug:
            print "ExAC allele frequency not found!"

    # Get the Complete Genomics allele frequency from the VCF info
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   Complete Genomics allele frequency
    def get_complete_freq(self, fields):
        for f in fields:
            if f.startswith("cg46"):
                return f.split("=")[1]
        if self.debug:
            print "Complete Genomics allele frequency not found!"

    # Check if the thousand genomes allele frequency is less than a provided cutoff
    # Input:
    #   cutoff: the maximum allele frequency allowed
    # Return:
    #   True or False based on meeting criteria
    def pass_1kg(self, cutoff):
        if cutoff is None or self.k1g_freq < cutoff and self.k1g_freq != '.':
            return True
        return False

    # Check if the esp allele frequency is less than a provided cutoff
    # Input:
    #   cutoff: the maximum allele frequency allowed
    # Return:
    #   True or False based on meeting criteria
    def pass_esp(self, cutoff):
        if cutoff is None or self.esp_freq < cutoff and self.k1g_freq != '.':
            return True
        return False

    # Check if the ExAC allele frequency is less than a provided cutoff
    # Input:
    #   cutoff: the maximum allele frequency allowed
    # Return:
    #   True or False based on meeting criteria
    def pass_exac(self, cutoff):
        if cutoff is None or self.exac_freq < cutoff and self.k1g_freq != '.':
            return True
        return False

    # Check if the Complete Genomics allele frequency is less than a provided cutoff
    # Input:
    #   cutoff: the maximum allele frequency allowed
    # Return:
    #   True or False based on meeting criteria
    def pass_complete(self, cutoff):
        if cutoff is None or self.complete_freq < cutoff and self.k1g_freq != '.':
            return True
        return False

    # Check if the consequence matches a provided consequence
    # Input;
    #   consequence: desired consequence
    # Return:
    #   True or False based on meeting criteria
    def pass_consequence(self, consequence):
        if consequence is None or consequence in self.consequence:
            return True
        return False

    # Check if there are any samples that inherited this variant
    # Return:
    #   True or False based on meeting criteria
    def has_ids(self):
        if len(self.from_father) > 0 or len(self.from_mother) > 0:
            return True
        return False

    # Check if a provided id inherited this variant from the father
    # Input:
    #   id: Sample ID
    # Return:
    #   True or False based on meeting criteria
    def has_paternal(self, id):
        if id in self.from_father:
            return True
        return False

    # Check if a provided id inherited this variant from the mother
    # Input:
    #   id: Sample ID
    # Return:
    #   True or False based on meeting criteria
    def has_maternal(self, id):
        if id in self.from_mother:
            return True
        return False

# Parse the command line
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script identifies compound heterozygous variants from a vcf-like file.')
    parser.add_argument("-i", "--input",
                                help="Input file.  VCF-like file containing information about variant inheritance.")
    parser.add_argument("-p", "--ped", default=None, help="Ped file.")
    parser.add_argument("--thousand_freq", help="1000 Genomes maximum frequency")
    parser.add_argument("--exac_freq", help="ExAC maximum frequency")
    parser.add_argument("--esp_freq", help="esp maximum frequency")
    parser.add_argument("--complete_freq", help="Complete Genomics maximum frequency")
    parser.add_argument("--consequence", help="Variant consequence")
    parser.add_argument("-o", "--output", help="File output prefix.  Not required, but recommended to prevent "
                                               "overwriting of old files")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    if options.input is None:
        print "Input required.  --input \n\n"
        parser.print_help()
        exit(0)
    return options

# Main
if __name__ == "__main__":
    options = parse_command_line()
    CompoundHets(options.input, options.ped, options.thousand_freq, options.exac_freq, options.esp_freq,
                 options.complete_freq, options.consequence, options.output)