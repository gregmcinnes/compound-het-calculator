from __future__ import division
from __future__ import print_function
from collections import defaultdict

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
    def __init__(self, chr, pos, ref, alt, info, from_mother_affected, from_mother_unaffected,
                 from_father_affected, from_father_unaffected, unknown_affected, unknown_unaffected,
                 missing_affected, missing_unaffected, uncertain_affected, uncertain_unaffected,
                 debug=False):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info
        self.from_mother_affected = self.parse_ids(from_mother_affected)
        self.from_mother_unaffected = self.parse_ids(from_mother_unaffected)
        self.from_father_affected = self.parse_ids(from_father_affected)
        self.from_father_unaffected = self.parse_ids(from_father_unaffected)
        self.unknown_affected = self.parse_ids(unknown_affected)
        self.unknown_unaffected = self.parse_ids(unknown_unaffected)
        self.missing_affected = self.parse_ids(missing_affected)
        self.missing_unaffected = self.parse_ids(missing_unaffected)
        self.uncertain_affected = self.parse_ids(uncertain_affected)
        self.uncertain_unaffected = self.parse_ids(uncertain_unaffected)
        self.debug = debug
        fields = self.info.split(";")
        self.k1g_freq = self.get_1kg_freq(fields)
        self.esp_freq = self.get_esp_freq(fields)
        self.exac_freq = self.get_exac_freq(fields)
        self.complete_freq = self.get_complete_freq(fields)
        self.gene, self.consequence, self.transcript = self.get_gene_consequence_and_transcript(fields)
        self.id = "%s_%s_%s_%s" % (self.chr, self.pos, self.ref, self.alt)
        self.polyphen_hdiv_pred = self.get_polyphen_hdiv_pred(fields)
        self.polyphen_hdiv_score = self.get_polyphen_hdiv_score(fields)
        self.polyphen_hvar_pred = self.get_polyphen_hvar_pred(fields)
        self.polyphen_hvar_score = self.get_polyphen_hvar_score(fields)
        self.cadd = self.get_cadd(fields)
        self.cadd_phred = self.get_cadd_phred(fields)

    def print_short(self):
        print("%s\t%s\t%s\t%s" % (self.chr, self.pos, self.gene, self.consequence))

    # Print the variant to the console.  Debugging method.
    def print_variant(self):
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chr, self.pos, self.gene, ",".join(self.consequence),
                                                          self.k1g_freq, self.esp_freq,
                                                          self.exac_freq, self.complete_freq,
                                                          ",".join(self.from_mother_affected),
                                                          ",".join(self.from_father_affected)))

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
        string = "%s\t%s\t%s\t%s\t%s\t%s" % (pair, self.chr, self.pos, self.gene, ",".join(self.from_mother_affected),
                                          ",".join(self.from_father_affected))
        return string

    # Get the gene symbol, consequence, and transcript from the info field
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   gene: gene symbol
    #   consequence: variant consequence
    #   transcript: transcript ID
    def get_gene_consequence_and_transcript(self, fields):
        for f in fields:
            if f.startswith("CSQ"):
                info = f.split("|")
                consequence = info[1].split('&')
                transcript = info[6]
                gene = info[3]
                #print(transcript)
                return gene, consequence, transcript
        if self.debug:
            print("Gene and consequence could not be identified!")

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
            print("1,000 Genomes allele frequency not found!")

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
            print("esp allele frequency not found!")

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
                print("ExAC allele frequency not found!")

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
            print("Complete Genomics allele frequency not found!")

    # Get the Polyphen HDIV prediction
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   Polyphen HDIV prediction
    def get_polyphen_hdiv_pred(self, fields):
        for f in fields:
            if f.startswith("Polyphen2_HDIV_pred"):
                return f.split("=")[1]
        if self.debug:
            print("Polyphen2 HDIV prediction not found!")

    # Get the Polyphen HDIV score
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   Polyphen HDIV score
    def get_polyphen_hdiv_score(self, fields):
        for f in fields:
            if f.startswith("Polyphen2_HDIV_score"):
                return f.split("=")[1]
        if self.debug:
            print("Polyphen2 HDIV score not found!")

    # Get the Polyphen HVAR prediction
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   Polyphen HVAR prediction
    def get_polyphen_hvar_pred(self, fields):
        for f in fields:
            if f.startswith("Polyphen2_HVAR_pred"):
                return f.split("=")[1]
        if self.debug:
            print("Polyphen2 HVAR prediction not found!")

    # Get the Polyphen HVAR score
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   Polyphen HVAR score
    def get_polyphen_hvar_score(self, fields):
        for f in fields:
            if f.startswith("Polyphen2_HVAR_score"):
                return f.split("=")[1]
        if self.debug:
            print("Polyphen2 HVAR score not found!")

    # Get the CADD score
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   CADD score
    def get_cadd(self, fields):
        for f in fields:
            if f.startswith("CADD="):
                return f.split("=")[1]
        if self.debug:
            print("CADD not found!")

    # Get the CADD Phred score
    # Input:
    #   fields: A list of strings from the VCF info
    # Return:
    #   CADD Phred score
    def get_cadd_phred(self, fields):
        for f in fields:
            if f.startswith("CADD_Phred"):
                return f.split("=")[1]
        if self.debug:
            print("CADD Phred not found!")


    # Check if the variant meets all the requirements: frequency cutoffs, consequence, and transcript
    # Input:
    #   frequency_cutoff: A float that serves as the maximum value for allele frequency
    #   consequence: A string that the variant consequence must match
    #   transcript: The transcript that the variant must match
    # Return:
    #   CADD score
    def meets_requirements(self, frequency_cutoff, consequence, transcript):
        if self.pass_allele_frequency(frequency_cutoff) and \
                self.pass_consequence(consequence) and self.has_ids():
                #and self.match_transcript(transcript):
            return True
        return False

    # Check if the all allele frequencies are less than a provided cutoff
    # Input:
    #   cutoff: the maximum allele frequency allowed
    # Return:
    #   True or False based on meeting criteria
    def pass_allele_frequency(self, cutoff):
        if cutoff is None:
            return True
        if self.pass_complete(cutoff) and self.pass_exac(cutoff) and self.pass_esp(cutoff) and self.pass_1kg(cutoff):
            return True
        return False

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
    # Input:
    #   consequence: desired consequence
    # Return:
    #   True or False based on meeting criteria
    def pass_consequence(self, consequence):
        if consequence is None:
            return True
        for c in consequence:
            if c in self.consequence:
                return True
        return False

    # Check if the transcript matches a provided transcript
    # Input:
    #   transcript: desired transcript
    # Return:
    #   True or False based on meeting criteria
    def match_transcript(self, transcript):
        if transcript is None or transcript == self.transcript:
            return True
        return False

    # Check if there are any samples that inherited this variant
    # Return:
    #   True or False based on meeting criteria
    def has_ids(self):
        if len(self.from_father_affected) > 0 or len(self.from_mother_affected) > 0:
            return True
        return False

    # Check if a provided id inherited this variant from the father
    # Input:
    #   id: Sample ID
    # Return:
    #   True or False based on meeting criteria
    def has_paternal(self, id):
        if id in self.from_father_affected:
            return True
        return False

    # Check if a provided id inherited this variant from the mother
    # Input:
    #   id: Sample ID
    # Return:
    #   True or False based on meeting criteria
    def has_maternal(self, id):
        if id in self.from_mother_affected + self.from_mother_unaffected:
            return True
        return False

    # Get the status of a patient id from the variant.
    # Input:
    #   id: Patient ID
    # Returns:
    #   A string for the origin of the variant
    def get_status(self, id):
        if id in self.from_mother_affected:
            return "MOTHER_AFFECTED"
        if id in self.from_mother_unaffected:
            return "MOTHER_UNAFFECTED"
        if id in self.from_father_affected:
            return "FATHER_AFFECTED"
        if id in self.from_father_unaffected:
            return "FATHER_AFFECTED"
        if id in self.unknown_affected:
            return "UNKNOWN_AFFECTED"
        if id in self.unknown_unaffected:
            return "UNKNOWN_UNAFFECTED"
        if id in self.missing_affected:
            return "MISSING_AFFECTED"
        if id in self.missing_unaffected:
            return "MISSING_UNAFFECTED"
        if id in self.uncertain_affected:
            return "UNCERTAIN_AFFECTED"
        if id in self.uncertain_unaffected:
            return "UNCERTAIN_UNAFFECTED"
        return None

# Patient class
# This class stores information about individual patients.  Each patient has a family id, patient id, id of the father,
#   id of the mother, sex, and disease status
# Input:
#   family_id: family ID
#   id: patient ID
#   paternal: paternal ID
#   maternal: maternal ID
#   sex: string representing gender
#   phenotype: 1 for non-disease, 2 for disease
class Patient(object):
    def __init__(self, family_id, id, paternal, maternal, sex, phenotype):
        self.family_id = family_id
        self.id = id
        self.paternal = paternal
        self.maternal = maternal
        self.sex = sex
        self.phenotype = phenotype

    # Check if the patient is a child.  If it is missing ids for both parents, we assume it is a child.
    def is_child(self):
        if self.maternal == '0' and self.paternal == '0':
            return False
        return True

    # Print out information about a patient - debugging function
    def print_patient(self):
        print("%s\t%s\t%s\t%s\t%s\t%s" % (self.family_id, self.id, self.paternal, self.maternal, self.sex, self.phenotype))

    # Return true if the patient has disease
    def has_disease(self):
        if self.phenotype == 2:
            return True
        return False

# Family class
# Store information about each family
# Class can be loaded with 2 variants in order to determine which members of the family are carriers of the compound het
# Input:
#   family_id: Family ID
class Family(object):
    def __init__(self, family_id):
        self.family_id = family_id
        self.members = {}
        self.statuses = {}
        self.initialized = False
        self.affected_carriers_count = None
        self.affected_family_count = None
        self.unaffected_family_count = None
        self.unaffected_carriers_count = None
        self.affected_missing_count = None
        self.unaffected_missing_count = None
        self.affected_uncertain_count = None
        self.unaffected_uncertain_count = None
        self.member_statuses = None

    # Add a family member
    # Input:
    #   id: Patient ID
    #   patient: Patient object
    def add_member(self, id, patient):
        self.members[id] = patient

    # Return the patient object given a patient ID
    # Input:
    #   id: Patient id
    # Return:
    #   Patient object
    def get_member(self, id):
        return self.members[id]

    # Print out information about the family - debugging function
    def print_family(self):
        print("FAMILY %s" % self.family_id)
        for m in self.members:
            self.members[m].print_patient()
        print()

    # Count the number of children for the family
    def child_count(self):
        count = 0
        for m in self.members:
            if self.members[m].is_child():
                count += 1
        return count

    # Count the number of family members that are affected with disease
    def affected_count(self):
        count = 0
        for m in self.members:
            if self.members[m].has_disease():
                count += 1
        self.affected_family_count = count
        return count

    # Count the number of family members that do not have disease
    def unaffected_count(self):
        count = 0
        for m in self.members:
            if not self.members[m].has_disease():
                count += 1
        self.unaffected_family_count = count
        return count

    # Given two variants, initialize compound het carrier statuses
    def initialize_variants(self, variant_1, variant_2):
        if self.initialized:
            self.reset_variant_counts()
        self.get_status_counts(variant_1)
        self.get_status_counts(variant_2)
        self.classify_family_members()
        self.get_number_missing_affected()
        self.get_number_missing_unaffected()
        self.get_number_uncertain_affected()
        self.get_number_uncertain_unaffected()
        self.get_number_affected_carriers()
        self.get_number_unaffected_carriers()
        self.get_number_noncarrier_affected()
        self.get_number_noncarrier_unaffected()
        self.initialized = True

    # Determine the carrier status of each family member for a given variant.  Store the statuses in a dictionary
    def get_status_counts(self, variant):
        statuses = defaultdict(list)
        for m in self.members:
            status = variant.get_status(m)
            statuses[status].append(m)
        self.statuses[variant.id] = statuses

    # Determine the number of family members that are MISSING and affected by disease for the compound het
    def get_number_missing_affected(self):
        # Rule
        # A missing call for either variant
        count = 0
        for m in self.members:
            patient = self.members[m]
            if patient.has_disease() and self.member_statuses[m] == 'MISSING':
                count += 1
        self.affected_missing_count = count
        return count

    # Determine the number of family members that are MISSING for the compound het and are not affected by disease
    def get_number_missing_unaffected(self):
        # Rule
        # A missing call for either variant plus a carrier
        count = 0
        for m in self.members:
            patient = self.members[m]
            if not patient.has_disease() and self.member_statuses[m] == 'MISSING':
                count += 1
        self.unaffected_missing_count = count
        return count

    # Determine the number of family members who are UNCERTAIN for the compound het and are affected by disease
    def get_number_uncertain_affected(self):
        # Rule
        # An uncertain call plus a carrier
        count = 0
        for m in self.members:
            patient = self.members[m]
            if patient.has_disease() and self.member_statuses[m] == 'UNCERTAIN':
                count += 1
        self.affected_uncertain_count = count
        return count

    # Determine the number of family members who are UNCERTAIN for disease and are not affected by disease
    def get_number_uncertain_unaffected(self):
        # Rule
        # An uncertain call plus a carrier
        count = 0
        for m in self.members:
            patient = self.members[m]
            if not patient.has_disease() and self.member_statuses[m] == 'UNCERTAIN':
                count += 1
        self.unaffected_uncertain_count = count
        return count

    # Determine the number of family members who are CARRIERs of the compound het and are affected by disease
    def get_number_affected_carriers(self):
        # Rule
        # Affected and carrier of both variants
        count = 0
        for m in self.members:
            patient = self.members[m]
            if patient.has_disease() and self.member_statuses[m] == 'CARRIER':
                count += 1
        self.affected_carriers_count = count
        return count

    # Determine the nuber of family members who are CARRIERs of the compound het and are not affected by disease
    def get_number_unaffected_carriers(self):
        # Rule
        # Affected and carrier of both variants
        count = 0
        for m in self.members:
            patient = self.members[m]
            if not patient.has_disease() and self.member_statuses[m] == 'CARRIER':
                count += 1
        self.unaffected_carriers_count = count
        return count

    # Determine the number of family members who are NONCARRIERs of the compound het and are affected by disease
    def get_number_noncarrier_affected(self):
        # Rule
        # Affected and carrier of both variants
        count = 0
        for m in self.members:
            patient = self.members[m]
            if patient.has_disease() and self.member_statuses[m] == 'NONCARRIER':
                count += 1
        return count

    # Determine the number of family members who are NONCARRIERs of the compound het and are unaffected by disease
    def get_number_noncarrier_unaffected(self):
        # Rule
        # Affected and carrier of both variants
        count = 0
        for m in self.members:
            patient = self.members[m]
            if not patient.has_disease() and self.member_statuses[m] == 'NONCARRIER':
                count += 1
        return count

    # Get the intersection of two lists.  If there are duplicates in the first list that will be represented the result
    def intersection(self, list_a, list_b):
        intersection = []
        for i in list_a:
            for n in range(0, list_b.count(i)):
                intersection.append(i)
        return intersection

    # Determine carrier status for each family member according the the rules below
        # Classification Rules
        # State 1   State 2   Classification
        # 1. Carrier	Carrier	Carrier
        # 2. Noncarrier	Noncarrier	Noncarrier
        # 3. Carrier	Missing	Missing
        # 4. Carrier	Uncertain	Uncertain
        # 5. Missing	Carrier	Missing
        # 6. Uncertain	Carrier	Uncertain
        # 7. Noncarrier	Missing	Noncarrier
        # 8. Noncarrier	Uncertain	Noncarrier
        # 9. Missing	Noncarrier	Noncarrier
        # 10. Uncertain	Noncarrier	Noncarrier
        # 11. Missing   Missing Missing
    def classify_family_members(self):
        carrier = 'CARRIER'
        noncarrier = 'NONCARRIER'
        missing = 'MISSING'
        uncertain = 'UNCERTAIN'

        member_statuses = {}
        for m in self.members:
            statuses = []
            for v in self.statuses:
                for s in self.statuses[v]:

                    if m in self.statuses[v][s]:
                        statuses.append(s)

            # 1. Carrier	Carrier	Carrier

            mother_statuses = ['MOTHER_AFFECTED', 'MOTHER_UNAFFECTED']
            father_statuses = ['FATHER_AFFECTED', 'FATHER_UNAFFECTED']
            parent_inheritance = ['FATHER_AFFECTED', 'MOTHER_AFFECTED', 'FATHER_UNAFFECTED', 'MOTHER_UNAFFECTED']
            missing_statuses = ['MISSING_AFFECTED', 'MISSING_UNAFFECTED']
            uncertain_statuses = ['UNCERTAIN_AFFECTED', 'UNCERTAIN_UNAFFECTED', 'UNKNOWN_AFFECTED', 'UNKNOWN_UNAFFECTED']

            # todo - there are many cases of mismatches between affected and unaffected.  Look into this
            if len(self.intersection(mother_statuses, statuses)) == 1 and len(self.intersection(father_statuses, statuses)) == 1:
                member_statuses[m] = carrier

            # 2. Noncarrier	Noncarrier	Noncarrier
            elif statuses.count(None) == 2:
                member_statuses[m] = noncarrier

            # 3 & 5
            # 3. Carrier	Missing	Missing
            elif len(self.intersection(parent_inheritance, statuses)) == 1 and len(self.intersection(missing_statuses, statuses)) == 1:
                member_statuses[m] = missing

            # 4 & 6
            # 4. Carrier	Uncertain	Uncertain
            elif len(self.intersection(parent_inheritance, statuses)) == 1 and len(self.intersection(uncertain_statuses, statuses)) == 1:
                member_statuses[m] = uncertain

            # 7 & 9
            # 7. Noncarrier	Missing	Noncarrier
            elif None in statuses and \
                    (len(self.intersection(missing_statuses, statuses)) == 1 or
                            len(self.intersection(uncertain_statuses, statuses)) == 1):
                member_statuses[m] = noncarrier

            # 11. Missing   Missing Missing
            elif len(self.intersection(missing_statuses, statuses)) == 2:
                member_statuses[m] = missing

            elif len(self.intersection(uncertain_statuses, statuses)) == 2:
                member_statuses[m] = uncertain

            elif statuses.count(None) == 1 and len(self.intersection(parent_inheritance, statuses)) == 1:
                member_statuses[m] = noncarrier

            else:
                print('Unable to determine status!')
                print(statuses)
                print(m)
                self.members[m].print_patient()
                #exit(1)  #todo error catch
                member_statuses[m] = 'UNDETERMINED'

        self.member_statuses = member_statuses

    # Determine the fraction of disease affected family members that are carriers of the compound het
    def get_fraction_of_affected(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.affected_family_count:
            return None
        return round(self.affected_carriers_count / self.affected_family_count, 2)

    # Determine the fraction of family members that are unaffected by disease that are carriers of the compound het
    def get_fraction_of_unaffected(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.unaffected_family_count:
            return None
        return round(self.unaffected_carriers_count / self.unaffected_family_count, 2)

    # Determine the fraction of family members that are affected by disease and carriers of the compound het.
    #   Adjust total count by the number of disease affected family members that have a MISSING classification for the
    #   het.
    def get_fraction_of_affected_missing_adj(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.affected_family_count:
            return None
        return round(self.affected_carriers_count / (self.affected_family_count - self.affected_missing_count), 2)

    # Determine the fraction of family members that are unaffected by disease that are carriers of the compound het
    #   Adjust total count by the number of disease unaffected family members that have a MISSING classification for the
    #   het.
    def get_fraction_of_unaffected_missing_adj(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.unaffected_family_count:
            return None
        return round(self.unaffected_carriers_count / (self.unaffected_family_count - self.unaffected_missing_count), 2)

    # Determine the fraction of family members that are affected by disease and carriers of the compound het.
    #   Adjust total count by the number of disease affected family members that have a UNCERTAIN classification for the
    #   het.
    def get_fraction_of_affected_uncertain_adj(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.affected_family_count:
            return None
        return round(self.affected_carriers_count / (self.affected_family_count - self.affected_uncertain_count), 2)

    # Determine the fraction of family members that are unaffected by disease that are carriers of the compound het
    #   Adjust total count by the number of disease unaffected family members that have a UNCERTAIN classification for the
    #   het.
    def get_fraction_of_unaffected_uncertain_adj(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.unaffected_family_count:
            return None
        return round(self.unaffected_carriers_count / (self.unaffected_family_count - self.unaffected_uncertain_count), 2)

    # Determine the fraction of family members that are affected by disease and carriers of the compound het.
    #   Adjust total count by the number of disease affected family members that have a MISSING or UNCERTAIN
    #   classification for the het.
    def get_fraction_of_affected_all_adj(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.affected_family_count:
            return None
        return round(self.affected_carriers_count / (self.affected_family_count - self.affected_missing_count- self.affected_uncertain_count), 2)

    # Determine the fraction of family members that are unaffected by disease that are carriers of the compound het
    #   Adjust total count by the number of disease unaffected family members that have a MISSING or UNCERTAIN
    #   classification for the het.
    def get_fraction_of_unaffected_all_adj(self):
        if not self.initialized:
            raise Exception("Variants uninitialized!")
        if not self.unaffected_family_count:
            return None
        return round(self.unaffected_carriers_count / (self.unaffected_family_count - self.unaffected_missing_count - self.unaffected_uncertain_count), 2)

    # Reset all the counts and variants for family
    def reset_variant_counts(self):
        self.statuses = {}
        self.initialized = False
        self.affected_carriers_count = None
        self.affected_family_count = None
        self.unaffected_family_count = None
        self.unaffected_carriers_count = None
        self.affected_missing_count = None
        self.unaffected_missing_count = None
        self.affected_uncertain_count = None
        self.unaffected_uncertain_count = None
        self.member_statuses = None

# CompoundHet class
# Store information about a single compound het.  Stored are the patient the compound het was identified in, the two
#   variants, and the family dictionary
class CompoundHet(object):
    def __init__(self, patient_id, maternal, paternal, families):
        self.patient_id = patient_id
        self.maternal = maternal
        self.paternal = paternal
        self.families = families

    # Get the family ID from the patient ID
    def patient_id_to_family(self, id):
        return id[0:6]

    # Print the compound het info to the aggregate file
    def print_aggregate(self, output):
            patient_id = self.patient_id
            paternal_variant = self.paternal
            maternal_variant = self.maternal
            gene = paternal_variant.gene
            maternal_id = maternal_variant.id
            maternal_consequence = "|".join(maternal_variant.consequence)
            paternal_id = paternal_variant.id
            paternal_consequence = "|".join(paternal_variant.consequence)
            joint_consequence = "%s-%s" % (maternal_consequence, paternal_consequence)
            family_id = self.patient_id_to_family(patient_id)
            family = self.families[family_id]

            children_count = family.child_count()
            affected_count = family.affected_count()
            unaffected_count = family.unaffected_count()

            family.initialize_variants(paternal_variant, maternal_variant)

            affected_carrier_count = family.get_number_affected_carriers()
            unaffected_carrier_count = family.get_number_unaffected_carriers()
            affected_missing_count = family.get_number_missing_affected()
            unaffected_missing_count = family.get_number_missing_unaffected()
            affected_uncertain_count = family.get_number_uncertain_affected()
            unaffected_uncertain_count = family.get_number_uncertain_unaffected()
            affected_noncarrier_count = family.get_number_noncarrier_affected()
            unaffected_noncarrier_count = family.get_number_noncarrier_unaffected()
            fraction_affected = family.get_fraction_of_affected()
            fraction_unaffected = family.get_fraction_of_unaffected()
            fraction_affected_missing_adj = family.get_fraction_of_affected_missing_adj()
            fraction_unaffected_missing_adj = family.get_fraction_of_unaffected_missing_adj()
            fraction_affected_uncertain_adj = family.get_fraction_of_affected_uncertain_adj()
            fraction_unaffected_uncertain_adj = family.get_fraction_of_unaffected_uncertain_adj()
            fraction_affected_all_adj = family.get_fraction_of_affected_all_adj()
            fraction_unaffected_all_adj = family.get_fraction_of_affected_all_adj()

            file = open(output, 'a')
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
                  (gene,	#Gene
                     maternal_id,	#Maternal_VarID
                     maternal_consequence,	#Maternal_Var_CSQ
                     paternal_id,	#Paternal_VarID
                     paternal_consequence,	#Paternal_Var_CSQ
                     joint_consequence,	#Maternal_Var_CSQ-Paternal_Var_CSQ
                     family_id,	#family
                     children_count,	#fam_n_children
                     affected_count,	#fam_n_aff
                     unaffected_count,	#fam_n_unaff
                     affected_missing_count,	#fam_n_missing_aff
                     unaffected_missing_count,	#fam_n_missing_unaff
                     affected_uncertain_count,	#n_uncertain_aff
                     unaffected_uncertain_count,	#n_uncertain_unaff
                     affected_carrier_count,	#n_aff_carriers
                     unaffected_carrier_count,	#n_unaff_carriers
                     affected_noncarrier_count,	#n_noncarrier_aff
                     unaffected_noncarrier_count,	#n_noncarrier_unaff
                     fraction_affected,	#frac_of_aff [n_aff_carriers/fam_n_aff]
                     fraction_unaffected,	#frac_of_unaff [n_unaff_carriers/fam_n_unaff]
                     fraction_affected_missing_adj,	#frac_of_aff_missing_adjusted [n_aff_carriers/(fam_n_aff - n_missing_aff)]
                     fraction_unaffected_missing_adj,	#frac_of_unaff_missing_adjusted [n_unaff_carriers/(fam_n_unaff - n_missing_unaff)]
                     fraction_affected_uncertain_adj,	#frac_of_aff_uncertain_adjusted [n_aff_carriers/(fam_n_aff - n_uncertain_aff)]
                     fraction_unaffected_uncertain_adj,	#frac_of_unaff_uncertain_adjusted [n_unaff_carriers/(fam_n_unaff - n_uncertain_unaff)]
                     fraction_affected_all_adj,	#frac_of_aff_carriers_missing_uncertain_adj [n_aff_carriers/(fam_n_aff - n_missing_aff - n_uncetain_aff)]
                     fraction_unaffected_all_adj,	#frac_of_unaff_carriers_missing_uncertain_adj [n_unaff_carriers/(fam_n_unaff - n_missing_unaff - n_uncetain_unaff)]
                     paternal_variant.info,	#Info_Variant1
                     maternal_variant.info	#Info_Variant2
                  ), file=file)
            file.close()
            family.reset_variant_counts()

    # Print the compound het info for individual families
    def print_family(self, output):
        paternal_variant = self.paternal
        maternal_variant = self.maternal

        gene = paternal_variant.gene

        maternal_consequence = "|".join(maternal_variant.consequence)
        paternal_consequence = "|".join(paternal_variant.consequence)
        joint_consequence = "%s-%s" % (maternal_consequence, paternal_consequence)

        patient_id = self.patient_id
        family_id = self.patient_id_to_family(patient_id)

        family = self.families[family_id]
        family.initialize_variants(paternal_variant, maternal_variant)

        file = open(output + '/' + family_id + '.family_file.tsv', 'a')
        for m in family.member_statuses:
            if family.member_statuses[m] == 'CARRIER':
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
                          (
                              gene,  #Gene
                              joint_consequence,  #Maternal_Var_CSQ-Paternal_Var_CSQ
                              family_id,  #family
                              m,  #child_id
                              family.members[m].has_disease(),  #is_aff
                              'FROM_MOTHER', #Inheritance_Variant1
                              maternal_variant.chr,  #Chr
                              maternal_variant.pos,  #Position
                              maternal_variant.ref,  #Ref
                              maternal_variant.alt,  #Alt
                              maternal_variant.id,  #VariantID
                              maternal_variant.esp_freq,  #esp6500siv2_all
                              maternal_variant.exac_freq,  #ExAC_ALL
                              maternal_variant.k1g_freq,  #ThousandGenomes_2014oct_all
                              maternal_variant.complete_freq,  #cg46
                              maternal_variant.cadd,  #CADD
                              maternal_variant.cadd_phred,  #CADD_Phred
                              maternal_variant.polyphen_hdiv_score,  #Polyphen2_HDIV_score
                              maternal_variant.polyphen_hdiv_pred,  #Polyphen2_HDIV_pred
                              maternal_variant.polyphen_hvar_score,  #Polyphen2_HVAR_score
                              maternal_variant.polyphen_hvar_pred,  #Polyphen2_HVAR_pred
                              maternal_consequence,  #consequence
                              'FROM_FATHER',  #Inheritance_Variant2
                              paternal_variant.chr,  #Chr
                              paternal_variant.pos,  #Position
                              paternal_variant.ref,  #Ref
                              paternal_variant.alt,  #Alt
                              paternal_variant.id,  #VariantID
                              paternal_variant.esp_freq,  #esp6500siv2_all
                              paternal_variant.exac_freq,  #ExAC_ALL
                              paternal_variant.k1g_freq,  #ThousandGenomes_2014oct_all
                              paternal_variant.complete_freq,  #cg46
                              paternal_variant.cadd,  #CADD
                              paternal_variant.cadd_phred,  #CADD_Phred
                              paternal_variant.polyphen_hdiv_score,  #Polyphen2_HDIV_score
                              paternal_variant.polyphen_hdiv_pred,  #Polyphen2_HDIV_pred
                              paternal_variant.polyphen_hvar_score,  #Polyphen2_HVAR_score
                              paternal_variant.polyphen_hvar_pred,  #Polyphen2_HVAR_pred
                              paternal_consequence,  #consequence
                              maternal_variant.info,  #Info_Variant1
                              paternal_variant.info  #Info_Variant2
                          ),
                      file=file
                )
        file.close()
        family.reset_variant_counts()
