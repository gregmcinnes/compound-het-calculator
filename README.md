# compound-het-calculator

Use compount-hets.py to identify samples with compount heterozygous variants.  Provide a VCF-like file with phasing information for which parent a variant was inherted from and the script will output files specifying all compound hets found for all samples, and a file for each sample containing only the compound hets for that sample.

## Example

#### Help statement
```
python compound-hets.py --help
usage: compound-hets.py [-h] [-i INPUT] [-p PED]
                        [--thousand_freq THOUSAND_FREQ]
                        [--exac_freq EXAC_FREQ] [--esp_freq ESP_FREQ]
                        [--complete_freq COMPLETE_FREQ]
                        [--consequence CONSEQUENCE] [-o OUTPUT] [-d]

This script identifies compound heterozygous variants from a vcf-like file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file. VCF-like file containing information about
                        variant inheritance.
  -p PED, --ped PED     Ped file.
  -t TRANSCRIPT, --transcript TRANSCRIPT
                        Transcript to match
  --thousand_freq THOUSAND_FREQ
                        1000 Genomes maximum frequency
  --exac_freq EXAC_FREQ
                        ExAC maximum frequency
  --esp_freq ESP_FREQ   esp maximum frequency
  --complete_freq COMPLETE_FREQ
                        Complete Genomics maximum frequency
  --consequence CONSEQUENCE
                        Variant consequence
  -o OUTPUT, --output OUTPUT
                        File output prefix. Not required, but recommended to
                        prevent overwriting of old files
  -d, --debug           Output debugging messages. May be very verbose.
```

#### Specifying options

```
python compount-hets.py --thousand_freq 0.01 -i chr22.tsv -o chr22
```

Output files

chr22.AU032504.compound_het_pairs.txt
chr22.AU032505.compound_het_pairs.txt
chr22.AU035208.compound_het_pairs.txt
chr22.AU037103.compound_het_pairs.txt
chr22.AU037104.compound_het_pairs.txt
chr22.AU066004.compound_het_pairs.txt
chr22.AU1000303.compound_het_pairs.txt
chr22.AU1000304.compound_het_pairs.txt
chr22.AU1072311.compound_het_pairs.txt
chr22.AU1164301.compound_het_pairs.txt
chr22.AU1164303.compound_het_pairs.txt
chr22.AU1174303.compound_het_pairs.txt
chr22.AU1274301.compound_het_pairs.txt
chr22.AU1274304.compound_het_pairs.txt
chr22.AU1393301.compound_het_pairs.txt
chr22.AU1393303.compound_het_pairs.txt
chr22.AU1397301.compound_het_pairs.txt
chr22.AU1397302.compound_het_pairs.txt
chr22.AU1397303.compound_het_pairs.txt
chr22.AU1399302.compound_het_pairs.txt
chr22.AU1438303.compound_het_pairs.txt
chr22.AU1438304.compound_het_pairs.txt
chr22.AU1608303.compound_het_pairs.txt
chr22.AU1796304.compound_het_pairs.txt
chr22.AU1802302.compound_het_pairs.txt
chr22.AU1886301.compound_het_pairs.txt
chr22.AU1886302.compound_het_pairs.txt
chr22.AU1886303.compound_het_pairs.txt
chr22.AU1921304.compound_het_pairs.txt
chr22.AU1921305.compound_het_pairs.txt
chr22.AU1953301.compound_het_pairs.txt
chr22.AU2263301.compound_het_pairs.txt
chr22.AU2263302.compound_het_pairs.txt
chr22.all_compound_het_pairs.txt

``` 
head chr22.AU032504.compound_het_pairs.txt

#PAIR	CHR	POS	GENE
1	22	38165269	TRIOBP
1	22	38119197	TRIOBP
```

```
head chr22.all_compound_het_pairs.txt

#PAIR	CHR	POS	GENE	FROM_MOTHER	FROM_FATHER
1	22	50687800	HDAC10	AU022705,AU022704	AU2492301,AU1335301,AU1399302,AU1511303,AU1399301,AU1312301,AU1335302,AU1312302,AU1433303,AU1312303
1	22	50685332	HDAC10	AU1399303,AU1796303,AU1796304,AU1399302	AU1674301,AU1674302,AU0331302,AU1772302,AU019706,AU1772303,AU0331301,AU1772301,AU1820303
2	22	50216725	BRD1	AU1159302,AU1592303,AU1159301,AU1592302	AU1886301,AU1886302,AU020703,AU1886303,AU016904
2	22	50217649	BRD1	AU1886303,AU1886301,AU1886302
3	22	46653285	PKDREJ		AU1393303,AU1393301
3	22	46652736	PKDREJ	AU1393306,AU1393303	AU1072302
4	22	45218339	ARHGAP8	AU1137201,AU1137202	AU1867303,AU0616302,AU0616301,AU1867302,AU035208
4	22	45255682	ARHGAP8	AU035208	AU037104,AU037103,AU1942303,AU1942302
5	22	40037051	CACNA1I	AU021005,AU021003	AU032504,AU032505,AU1802302,AU0958304,AU1802301

```

