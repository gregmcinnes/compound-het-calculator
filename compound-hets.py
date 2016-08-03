# Gregory McInnes
# gmcinnes@stanford.edu
# November 20, 2015

# compound-hets.py
# Identify compound heterozygous variants from a vcf-like file.  This script will look over each gene
# in the list and determine if there are any pairs of variants for each sample that were not present in either
# the mother or father.

# Questions
# Is one frequency required for all datasets or one for each? - yes
# Decide how to handle variants where there is no reported frequency for a provided metric - as long as one of them has one
# Some variants have multiple reported consequences.  How to handle these cases? - include as long as one of them matches
# Should frequency cutoffs be minimums or maximums? -
# What is the ped file for?

# Round 2
# Are there grandparents?
# Functionality check: I only output compound hets where there is at least one affected child.  Is this the desired result?


from __future__ import division
from __future__ import print_function
import argparse
import os
from collections import defaultdict
from utilities import Utilities
from classes import CompoundHet
import shutil


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
#   gene_lookback: How many genes to save in memory.  Default is 5
class CompoundHets(object):
    def __init__(self, input, ped, frequency_cutoff=None, consequence=None, output=None, prefix=None, transcript=False, debug=False,
                 gene_lookback=5, force=False):
        self.frequency_cutoff = frequency_cutoff
        self.consequence = consequence
        self.transcript = transcript
        self.output = output
        self.prefix = prefix
        self.force = force
        self.aggregate_file = os.path.join(output, 'compound_hets.aggregate.tsv')
        if self.prefix:
            self.aggregate_file = os.path.join(output, self.prefix + '.compound_hets.aggregate.tsv')
        self.debug = debug
        self.Utils = Utilities(self.frequency_cutoff, self.consequence, self.transcript, self.debug)
        self.families = self.Utils.read_ped_file(ped)
        self.lookback = gene_lookback
        self.setup_output(output, self.families)
        self.process_file(input)


    # Print each compound het pair for each sample to a sample file
    # Input:
    #   samples: dictionary of lists of tuples of Variants
    #   output: output prefix
    def print_samples(self, samples, output):
        if not samples:
            return
        for s in samples:
            filename = "%s.compound_het_pairs.txt" % s
            if self.prefix:
                filename = self.prefix + "." + filename
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

    # Parse each line in the input file, convert it to a Variant object, and assign it to a gene within the gene
    #  dictionary
    # Input:
    #   input: input filename
    # Return:
    #   genes: dictionary of lists of Variants
    def process_file(self, input):
        current_genes = defaultdict(list)
        gene_list = []
        checked_genes = []
        with open(input) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                variant = self.Utils.parse_line(line)
                if variant is None:
                    continue
                if variant.gene not in current_genes:
                    if variant.gene in checked_genes:
                        print("%s ALREADY CHECKED! File may be out of order!" % variant.gene)
                    else:
                        checked_genes.append(variant.gene)
                    gene_list.append(variant.gene)
                    # Keep a lookback of 5 genes to account for gene overlaps
                    if len(current_genes) > self.lookback:
                        to_process = gene_list.pop(0)
                        self.process_gene(current_genes[to_process])
                        del current_genes[to_process]
                current_genes[variant.gene].append(variant)
        for gene in current_genes:
            self.process_gene(current_genes[gene])

    # Get any het calls for a gene and print the results
    # Input:
    #   gene: A list of variants (Variant objects) within a gene
    # Return:
    #   Nothing
    def process_gene(self, gene):
        if self.debug:
            print("Processing gene")
        hets = self.check_gene(gene)
        if hets:
            self.print_hets(hets)

    # For each gene in the gene list, check whether each sample has a compound het pair within the gene
    # Input:
    #   genes: dictionary of lists of Variants
    # Return:
    #   samples: dictionary of lists of tuples of Variants
    #   positions: A list of tuples of variants
    def check_gene(self, gene):
        hets = []
        identified_pairs = []
        for v1 in gene:
            for patient_id in v1.from_father_affected + v1.from_father_unaffected:
                family_id = self.Utils.patient_id_to_family(patient_id)
                if not self.families[family_id].members[patient_id].has_disease():
                    continue
                for v2 in gene:
                    if v1.pos == v2.pos:
                        continue
                    if v2.has_maternal(patient_id):
                        key = "-".join(sorted([v1.pos, v2.pos, self.Utils.patient_id_to_family(patient_id)]))
                        het = CompoundHet(patient_id, v2, v1, self.families)
                        if not key in identified_pairs:
                            if self.transcript:
                                if v1.transcript == v2.transcript:
                                    hets.append(het)
                                    identified_pairs.append(key)
                            else:
                                hets.append(het)
                                identified_pairs.append(key)
        return hets

    # Print out the hets - debugging function
    # input:
    #   hets: a CompoundHet object
    def print_hets(self, hets):
        for h in hets:
            h.print_aggregate(self.aggregate_file)
            h.print_family(self.output, self.prefix)

    # Create output folder and initialize all output files with headers
    # input:
    #   output: the output directory
    #   families: A dictionary of Family objects.  Family IDs are the keys of the dictionary.
    def setup_output(self, output, families):
        # Due to popular demand this has been removed.
        #if os.path.exists(output) and not self.force:
        #    print("Output folder found!  Please choose a different output path.")
        #    exit(0)
        #elif os.path.exists(output) and self.force:
        #    print("Existing output folder being moved to %s" % output + ".old")
        #    if os.path.exists(output + ".old"):
        #        print("Removing existing old output: %s" % output + ".old")
        #        #os.remove(output + ".old")
        #        shutil.rmtree(output + ".old")
        #    os.rename(output, output + ".old")
        if not os.path.exists(output):
            os.mkdir(output)

        file = open(self.aggregate_file,'w')
        print('Gene	Maternal_VarID	Maternal_Var_CSQ	Paternal_VarID	Paternal_Var_CSQ	'
              'Maternal_Var_CSQ-Paternal_Var_CSQ	family	fam_n_children	fam_n_aff	fam_n_unaff	fam_n_missing_aff'
              '	fam_n_missing_unaff	n_uncertain_aff	n_uncertain_unaff	n_aff_carriers	n_unaff_carriers	'
              'n_noncarrier_aff	n_noncarrier_unaff	frac_of_aff frac_of_unaff	frac_of_aff_missing_adjusted	'
              'frac_of_unaff_missing_adjusted	frac_of_aff_uncertain_adjusted	frac_of_unaff_uncertain_adjusted	'
              'frac_of_aff_carriers_missing_uncertain_adj	frac_of_unaff_carriers_missing_uncertain_adj	'
              'Info_Variant1	Info_Variant2', file=file)
        file.close()

        for f in families:
            filename = f + '.family_file.tsv'
            if self.prefix:
                filename = self.prefix + "." + filename
            family_file = os.path.join(output, filename)
            file = open(family_file,'w')
            print('Gene	Maternal_Var_CSQ-Paternal_Var_CSQ	family	child_id	is_aff	Inheritance_Variant1	Chr	'
                  'Position	Ref	Alt	VariantID	esp6500siv2_all	ExAC_ALL	ThousandGenomes_2014oct_all	cg46	CADD'
                  '	CADD_Phred	Polyphen2_HDIV_score	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_pred'
                  '	consequence	Inheritance_Variant2	Chr	Position	Ref	Alt	VariantID	esp6500siv2_all	ExAC_ALL'
                  '	ThousandGenomes_2014oct_all	cg46	CADD	CADD_Phred	Polyphen2_HDIV_score	Polyphen2_HDIV_pred'
                  '	Polyphen2_HVAR_score	Polyphen2_HVAR_pred	consequence	Info_Variant1', file=file)
            file.close()

# Parse the command line
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script identifies and reports compound heterozygous variants.')
    parser.add_argument("-i", "--input", help="An input file listing genomic positions and the patients that have variants at each positions.  Each row represents a single genomic position, with the chromosome, locus, reference, and alternate alleles listed (just like a VCF).  Also like a VCF, each row must have an INFO column with information about which gene the position is in, allele frequency in popular genomics resources, and predicted allele effect.")
    parser.add_argument("-p", "--ped", help="A PED file detailing family and individual identification numbers, paternal IDs, sex, and phenotype for each patient.  For phenotype identification, 1 is considered unaffected, 2 is affected by disease.")
    parser.add_argument("-t", "--transcript", action='store_true', default=False, help="Optional argument.  If set, only variants occuring within the same transcript will be output.")
    parser.add_argument("--frequency", help="Optional argument.  The maximum allele frequency cutoff allowed.  1000 Genomes, ExAC, and Complete Genomics allele frequencies are considered.  If any of the three resources reports an allele frequency greater than the set cutoff the variant will not be considered. If not specified all variants will all allele frequencies will be considered.")
    parser.add_argument("--consequence", nargs='*', help="Optional argument. The desired variant consequence (e.g. missense_variant, nonsense_variant, etc).  Only variants with the specified consequence will be considered.  Multiple consequences may be specified.  If not specified all variants will be considered.")
    parser.add_argument("-o", "--output", help="The output directory.")
    parser.add_argument("--prefix", help="Prefix for all output files.  By default no prefix is set.")
    parser.add_argument("--lookback", default=5, help="Number of genes to save in local memory.  Default=5")
    parser.add_argument("--force", action='store_true', help="Overwrite existing output directory")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    if options.input is None:
        print("Input required.  --input \n\n")
        parser.print_help()
        exit(0)

    if options.ped is None:
        print("Ped file required.  --ped \n\n")
        parser.print_help()
        exit(0)

    if options.output is None:
        print("Output directory required.  --output \n\n")
        parser.print_help()
        exit(0)



    return options

# Main
if __name__ == "__main__":
    options = parse_command_line()
    CompoundHets(options.input, options.ped, options.frequency,
                 options.consequence, options.output, options.prefix, options.transcript, options.debug,
                 options.lookback, options.force)