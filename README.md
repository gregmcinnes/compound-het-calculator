# compound-het-calculator

Use compount-hets.py to identify samples with compount heterozygous variants.  Provide a VCF-like file with phasing information for which parent a variant was inherted from and the script will output files specifying all compound hets found for all samples, and a file for each sample containing only the compound hets for that sample.

## Installation

To install the software clone the git repository to you local machine or cluster.  Once you have cloned the repository it should be ready to use.

```
git clone https://github.com/gregmcinnes/compound-het-calculator 
```

## Options

```
$ python compound-hets.py --help
usage: compound-hets.py [-h] [-i INPUT] [-p PED] [-t TRANSCRIPT]
                        [--frequency FREQUENCY] [--consequence CONSEQUENCE]
                        [-o OUTPUT] [--lookback LOOKBACK] [-d]

This script identifies and reports compound heterozygous variants.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        An input file listing genomic positions and the
                        patients that have variants at each positions. Each
                        row represents a single genomic position, with the
                        chromosome, locus, reference, and alternate alleles
                        listed (just like a VCF). Also like a VCF, each row
                        must have an INFO column with information about which
                        gene the position is in, allele frequency in popular
                        genomics resources, and predicted allele effect.
  -p PED, --ped PED     A PED file detailing family and individual
                        identification numbers, paternal IDs, sex, and
                        phenotype for each patient. For phenotype
                        identification, 1 is considered unaffected, 2 is
                        affected by disease.
  -t TRANSCRIPT, --transcript TRANSCRIPT
                        Optional argument. The name of the gene transcript to
                        match. Only variants from the specified transcript
                        will be returned. Only one transcript allowed. If not
                        specified all transcripts will be considered.
  --frequency FREQUENCY
                        Optional argument. The maximum allele frequency cutoff
                        allowed. 1000 Genomes, ExAC, and Complete Genomics
                        allele frequencies are considered. If any of the three
                        resources reports an allele frequency greater than the
                        set cutoff the variant will not be considered. If not
                        specified all variants will all allele frequencies
                        will be considered.
  --consequence CONSEQUENCE
                        Optional argument. The desired variant consequence
                        (e.g. missense_variant, nonsense_variant, etc). Only
                        variants with the specified consequence will be
                        considered. Multiple consequences may be specified. If
                        not specified all variants will be considered.
  -o OUTPUT, --output OUTPUT
                        The output directory. The specified output directory
                        must not already exist. The created directory will
                        contain all relevant output files.
  --lookback LOOKBACK   Number of genes to save in local memory. Default=5
  -d, --debug           Output debugging messages. May be very verbose.
```

## Output

All output files are directed to the directory specified with --output.  Two types of files are output as variants are processed: an aggregate file containing all identified compound hets (compound_hets.aggregate.tsv), and family files for each family listing information about each family member and any compound hets they carry (FAMILY_ID.family_file.tsv).
