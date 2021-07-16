#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys, os

from v1_downloadRelevantGCF import downloadGCF
from v1_selectHumanGeneID import selectHumanGeneID
from v1_addGenomes import addGenome
from v1_checkSequences import checkDownloadedSequences

#-------------------------------------------------------------------------------------
'''
  Help Message
'''

def usage():

    sys.stderr.write('''

    Usage : Launch error.
    'buildOrtholog.py'
  
    This program build ortholog fasta files of orthologous genes using the human gene identifier as a cross reference and three core taxa.
    

    - assembly_summary.tsv : is supposed to have 22 fields including :
        assembly_accession  refseq_category taxid   organisme_name  ftp_path (to download the GCF files)
        (following NCBI convention: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/ assembly_summary.txt)
        WARNING : REMEMBER TO DELETE THE HYBRID TAXON (30522) FROM THE FILE (Bos indicus x Bos taurus)

    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id GeneID  relationship    Other_tax_id    Other_GeneID 
        (following NCBI convention: http://ftp.ncbi.nlm.nih.gov/gene/DATA/ gene-ortholog.gz )

    - core_species.list : should contain one taxon id per line, for a core made of Homo Sapiens, Mus musculus and Canis Lupus
    familiaris it will be (human should always be first):
        9606        (Homo sapiens)
        10090       (Mus musculus)
        9615        (Canis Lupus familiaris)

    Usefull links : 
      - https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice
      - https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.ncbi_s_annotation_displayed_on


    Usage :  (Requires five parameters)
    python3 buildOrtholog.py assembly_summary.txt core_species.list gene_orthologs.txt outputFolder(= work directory)

    ''')

#-------------------------------------------------------------------------------------
'''
  Used the 'GCF.list' file to launch the 'addGenome' function
'''

def buildFasta (GCF_list,outputFastaFolder, outputInfoFolder,humanGeneIDList):
    GCFFolder = os.path.dirname(GCF_list)
    with open (GCF_list, 'r') as GCF_list_File : # GCF.list parsing
        for line in GCF_list_File :
            taxInfo = line.strip().split("\t") # recovery of the differents fields of interest
            taxName = taxInfo[0]
            GCF_CDS_File = GCFFolder+os.path.sep+taxInfo[1]
            GCF_Prot_File = GCFFolder+os.path.sep+taxInfo[2]
            GCF_RNA_File = GCFFolder+os.path.sep+taxInfo[3]
            GCF_Translate_File = GCFFolder+os.path.sep+taxInfo[4]
            GCF_RNA_GBFF_File = GCFFolder+os.path.sep+taxInfo[5]
            taxID = taxInfo[6]
            taxClassif = taxInfo[7]
            
            print("---------------")
            print (taxName)
            print (GCF_CDS_File)
            print (GCF_Prot_File)
            print (GCF_RNA_File)
            print (GCF_Translate_File)
            print (GCF_RNA_GBFF_File)
            print (taxID)
            print (taxClassif)

            addGenome(humanGeneIDList, orthologFile, GCF_CDS_File, GCF_Prot_File, GCF_RNA_File, GCF_Translate_File, GCF_RNA_GBFF_File, taxID, taxName, taxClassif, outputFastaFolder, outputInfoFolder)

#-------------------------------------------------------------------------------------
'''
  MAIN FUNTION OF THE PROGRAM
    Call the following functions :
        - downloadGCF (in downloadRelevantGCF.py)
        - selectHumanGeneID (in selectHumanGeneID.py)
        - addGenome (in addGenomes.py) 
        - checkDownloadedSequences (checkSequences.py)
'''

if __name__ == "__main__" :
    if len(sys.argv) != 5 : # Requires five parameters
        usage()
        sys.exit()
    else :
        summaryFile = sys.argv[1]       # assembly_summary.txt (NCBI)
        coreTaxonIDFile = sys.argv[2]   # HomoMusCanis.id (core_species.list)
        orthologFile = sys.argv[3]      # gene_orthologs (NCBI)
        outputFolder = sys.argv[4]      # work directory

        outputFastaFolder = outputFolder+os.path.sep+"FASTA"
        outputInfoFolder = outputFolder+os.path.sep+"INFO"
        outputGCFFolder = outputFolder+os.path.sep+"GCF"
        outputGCFList = outputGCFFolder+os.path.sep+"GCF.list"
        humanGeneIDList = outputFolder+os.path.sep+"humanGeneID.list"
        
        os.system("mkdir "+outputFastaFolder)
        os.system("mkdir "+outputInfoFolder)
        os.system("mkdir "+outputGCFFolder)

        print ("Download GCF :")
        downloadGCF(summaryFile, outputGCFFolder, outputGCFList)
        print ("select Human GeneID :")
        selectHumanGeneID(orthologFile, coreTaxonIDFile,humanGeneIDList)
        print ("build Fasta :")
        buildFasta(outputGCFList, outputFastaFolder,outputInfoFolder, humanGeneIDList)
        print ("checkDownloadedSequences :")
        checkDownloadedSequences(outputFolder)


        # WARNING : REMEMBER TO DELETE THE HYBRID TAXON (30522) FROM THE 'ASSEMBLY_SUMMARY' FILE
        # (Bos indicus x Bos taurus)