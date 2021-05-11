#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys, os
sys.path.append(".")
from download_relevant_GCF import downloadGCF
from selectHumanGeneID import selectHumanGeneID
from addGenomes import addGenome




def helpMessage():
    print ("This program build ortholog fasta files \
usage: python buildOrtholog.py assembly_summary.tsv gene_orthologs.tsv core_species.list output_fasta_folder \n\
-gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs: \n\
   tax_id GeneID  relationship    Other_tax_id    Other_GeneID (following NCBI convention: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene-ortholog.gz )\n \
-core_species.list should contain one taxon id per line, for a core made of Homo Mus and Canis it will be (human should always be first):\n\
9606\n\
10090\n\
9615")

# ajouter origine summary file
#-------------------------------------------------------------------------------------

# G
def buildFasta (GCF_list,resFolder):
    with open (GCF_list, 'r') as GCF_list_File :
        for line in GCF_list_File :
            taxInfo = line.strip().split("\t")
            taxName = taxInfo[0]
            GCF_CDS_File = taxInfo[1]
            GCF_Prot_File = taxInfo[2]
            taxID = taxInfo[3]
            addGenome("humanGeneID.list", orthologFile, GCF_CDS_File, GCF_Prot_File, taxID, taxName, resFolder)


#-------------------------------------------------------------------------------------

summaryFile = sys.argv[1]
coreTaxonList = sys.argv[2]
orthologFile = sys.argv[3]
outputFastaFolder = sys.argv[4]


#downloadGCF(summaryFile, "GCF.list")
#selectHumanGeneID(orthologFile, coreTaxonList,"humanGeneID.list")
buildFasta("GCF.list", outputFastaFolder)