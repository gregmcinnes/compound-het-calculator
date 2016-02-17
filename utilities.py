from classes import Variant
from classes import Family
from classes import Patient


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

# Utilities class
#   Functions required for processing files
class Utilities(object):
    def __init__(self, frequency_cutoff, consequence, transcript, debug=False):
        self.frequency_cutoff = frequency_cutoff
        self.consequence = consequence
        self.transcript = transcript
        self.debug = debug

    # Return the family ID from the patient ID
    def patient_id_to_family(self, id):
        return id[0:6]

    # Read in the ped file.  Create a dictionary of Family objects for each family.
    def read_ped_file(self, file):
        # Format
        # ID	IID	Paternal	Maternal	Sex	Phenotype
        ID = 0
        IID = 1
        PATERNAL = 2
        MATERNAL = 3
        SEX = 4
        PHENOTYPE = 5
        families = {}
        with open(file) as f:
            for line in f:
                if line.startswith('ID'):
                    continue
                line = line.rstrip()
                fields = line.split("\t")
                family_id = fields[ID]
                id = fields[IID]
                patient = Patient(family_id, id, fields[PATERNAL],
                                  fields[MATERNAL], int(fields[SEX]), int(fields[PHENOTYPE]))
                if family_id in families:
                    families[family_id].add_member(id, patient)
                else:
                    family = Family(family_id)
                    families[family_id] = family
                    families[family_id].add_member(id, patient)
        if self.debug:
            for f in families:
                families[f].print_family()
        return families



    # Parse a line to determine chromosome, position, and other relevant information
    # Input:
    #   line: A line directly from the input file
    # Return:
    #   variant: Variant object. None if it does not meet the input criteria
    def parse_line(self, line):
        #print line
        line = line.rstrip()
        fields = line.split("\t")
        chr = fields[CHROMOSOME]
        pos = fields[POSITION]
        ref = fields[REF]
        alt = fields[ALT]
        info = fields[INFO]
        from_mother_affected = fields[FROM_MOTHER_IDS_AFFECTED]
        from_mother_unaffected = fields[FROM_MOTHER_IDS_UNAFFECTED]
        from_father_affected = fields[FROM_FATHER_IDS_AFFECTED]
        from_father_unaffected = fields[FROM_FATHER_IDS_UNAFFECTED]
        unknown_affected = fields[UNKNOWN_PHASE_IDS_AFFECTED]
        unknown_unaffected = fields[UNKNOWN_PHASE_IDS_UNAFFECTED]
        missing_affected = fields[MISSING_IDS_AFFECTED]
        missing_unaffected = fields[MISSING_IDS_UNAFFECTED]
        uncertain_affected = fields[UNCERTAIN_IDS_AFFECTED]
        uncertain_unaffected = fields[UNCERTAIN_IDS_UNAFFECTED]

        variant = Variant(chr, pos, ref, alt, info, from_mother_affected, from_mother_unaffected, from_father_affected,
                          from_father_unaffected, unknown_affected, unknown_unaffected, missing_affected, missing_unaffected,
                          uncertain_affected, uncertain_unaffected)
        if variant.meets_requirements(self.frequency_cutoff, self.consequence, self.transcript):
            return variant
        return None