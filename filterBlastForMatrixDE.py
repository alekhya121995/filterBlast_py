########################################################################
# Program: filterBlastForMatrixDE.py
# Name: Alekhya Akkunuri
# Date: 12-11-2017
# Description: This program parses the BLAST and DE file and filters the
# results. It performs transcript to protein lookup and displays the 
# entries that meet the criteria
#########################################################################

# define BLAST class
class Blast():
    def __init__(self, line):
        line_split = line.rstrip("\n").split("\t")
        # store each required element into a variable
        self.transcriptId_isoform = line_split[0]
        self.swissprotId = line_split[1]
        self.identity = line_split[2]
        # split the variables further
        self.subsplit1 = self.transcriptId_isoform.split("|")
        self.subsplit2 = self.swissprotId.split("|")
        self.transcript = self.subsplit1[0]
        self.swissprot = self.subsplit2[3]
      

# reading in the BLAST file and mapping it
blast_file = open("blastp.outfmt6").readlines()
blast_list = map(Blast, blast_file)

# function to accept and filter BLAST object, only those with pident > 95 will be considered
def filter_blast(blast_object):
    if(float(blast_object.identity) > 95):
       return True

# filtering the BLAST list
filtered_blast = filter(filter_blast, blast_list)

# dictionary comprehension 
transcript_swissprot = {blast_object.transcript:blast_object for blast_object in filtered_blast}

# define Matrix class
class Matrix():
    def __init__(self, record):
        matrix_split = record.rstrip("\n").split("\t")
        (self.matrix_transcript, self.sp_ds, self.sp_hs, self.sp_log, self.sp_plat) = matrix_split       
        # check for BLAST matches
        if(self.matrix_transcript in transcript_swissprot):
          self.protein = transcript_swissprot.get(self.matrix_transcript).swissprot
        else:
           self.protein = self.matrix_transcript

# reading and mapping de file
matrix_file = open("diffExpr.P1e-3_C2.matrix").readlines()
matrix_list = map(Matrix, matrix_file)

# function to accept tuple and return 
def formatted_tuple(obj):
    return("\t".join(obj))

# printing the results
for matrix_object in matrix_list:
    print(formatted_tuple((matrix_object.protein, matrix_object.sp_ds, matrix_object.sp_hs, matrix_object.sp_log, matrix_object.sp_plat)))
    
