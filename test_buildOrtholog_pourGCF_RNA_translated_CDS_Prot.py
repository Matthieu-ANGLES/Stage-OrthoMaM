#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys, os
#sys.path.append(".")
# Pourquoi sys.path ? Pour signifier la localisation des imports ?

from test_downloadRelevantGCF_pourGCF_RNA_translated_CDS_Prot import downloadGCF
from test_selectHumanGeneID_pourGCF_RNA_translated_CDS_Prot import selectHumanGeneID
from test_addGenomes_pourGCF_RNA_translated_CDS_Prot import addGenome
#from checkSequences import 


#-------------------------------------------------------------------------------------
'''
  Help Message
'''

def usage():
    '''
    print ("This program build ortholog fasta files \
    Usage: python buildOrtholog.py assembly_summary.tsv gene_orthologs.tsv core_species.list output_fasta_folder \n\
    -gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs: \n\
    tax_id GeneID  relationship    Other_tax_id    Other_GeneID (following NCBI convention: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene-ortholog.gz )\n \
    -core_species.list should contain one taxon id per line, for a core made of Homo Mus and Canis it will be (human should always be first):\n\
    9606\n\
    10090\n\
    9615")
    '''

    sys.stderr.write('''

    buildOrtholog.py
    
    This program build ortholog fasta files.
    
    Usage : python3 buildOrtholog.py assembly_summary.tsv core_species.list gene_orthologs.tsv output_fasta_folder
    
    - assembly_summary.tsv : is supposed to have 22 fields including :
        assembly_accession  refseq_category taxid   organisme_name  ftp_path (to download the GCF files)
        (following NCBI convention: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt)

    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id GeneID  relationship    Other_tax_id    Other_GeneID 
        (following NCBI convention: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene-ortholog.gz )

    -core_species.list should contain one taxon id per line, for a core made of Homo Sapiens, Mus musculus and Canis Lupus
    familiaris it will be (human should always be first):
        9606
        10090
        9615


    ''')


#-------------------------------------------------------------------------------------
'''
  
'''

def buildFasta (GCF_list,outputFastaFolder, outputInfoFolder,humanGeneIDList):
    GCFFolder = os.path.dirname(GCF_list)
    with open (GCF_list, 'r') as GCF_list_File :
        for line in GCF_list_File :
            print("---------------")
            taxInfo = line.strip().split("\t")
            taxName = taxInfo[0]
            print (taxName)
            GCF_CDS_File = GCFFolder+os.path.sep+taxInfo[1]
            print (GCF_CDS_File)
            GCF_Prot_File = GCFFolder+os.path.sep+taxInfo[2]
            print (GCF_Prot_File)
            #---ajout/modif---
            #taxID = taxInfo[3]
            GCF_RNA_File = GCFFolder+os.path.sep+taxInfo[3]
            print (GCF_RNA_File)
            GCF_Translate_File = GCFFolder+os.path.sep+taxInfo[4]
            print (GCF_Translate_File)
            GCF_RNA_GBFF_File = GCFFolder+os.path.sep+taxInfo[5]
            print (GCF_RNA_GBFF_File)
            taxID = taxInfo[6]
            print (taxID)
            #addGenome(humanGeneIDList, orthologFile, GCF_CDS_File, GCF_Prot_File, taxID, taxName, outputFastaFolder, outputInfoFolder)
            addGenome(humanGeneIDList, orthologFile, GCF_CDS_File, GCF_Prot_File, GCF_RNA_File, GCF_Translate_File, GCF_RNA_GBFF_File, taxID, taxName, outputFastaFolder, outputInfoFolder)

#-------------------------------------------------------------------------------------
'''
  Main
'''

if __name__ == "__main__" :
    if len(sys.argv) != 5 :
        usage()
        sys.exit()
    else :
        summaryFile = sys.argv[1]       # small_resume_assembly
        coreTaxonList = sys.argv[2]     # HomoMusCanis.id
        orthologFile = sys.argv[3]      # gene_orthologs
        outputFolder = sys.argv[4]      # verificationFichiers_210526/ data_runComplet

        outputFastaFolder = outputFolder+os.path.sep+"FASTA"
        outputInfoFolder = outputFolder+os.path.sep+"INFO"
        outputGCFFolder = outputFolder+os.path.sep+"GCF"
        outputGCFList = outputGCFFolder+os.path.sep+"GCF.list"
        humanGeneIDList = outputFolder+os.path.sep+"humanGeneID.list"
        
        os.system("mkdir "+outputFastaFolder)
        os.system("mkdir "+outputInfoFolder)
        os.system("mkdir "+outputGCFFolder)

        print ("Download GCF :")
        downloadGCF(summaryFile,outputGCFFolder, outputGCFList)
        print ("select Human GeneID :")
        selectHumanGeneID(orthologFile, coreTaxonList,humanGeneIDList)
        print ("build Fasta :")
        buildFasta(outputGCFList, outputFastaFolder,outputInfoFolder, humanGeneIDList)