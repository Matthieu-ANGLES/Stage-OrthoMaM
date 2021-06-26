#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import os, sys
from os import path as os_path

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''

    Usage : Launch error.
    'downloadRalevantGCF.py'

    Use the 'assembly_summary.txt' retrieved at the bottom of the list at this address
      - https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/

    
    This module generate the "GCF.list" file used by 'addGenomes.py' in the GCF File.
    It has 8 fields per line in tabulated format :
      - tax_name
      - GCF_CDS_File 
      - GCF_Prot_File 
      - GCF_RNA_File
      - GCF_Feature_Table_File
      - GCF_RNA_GBFF_File
      - genometaxID
      - taxonClassif
    
    This module also generates the 'GCF_inexistent.list' filein case of malfunction during the download.
    It has 9 fields per line in tabulated format : (add refseq FTP link)
      - tax_name
      - genomeFTPFolder     (refseq FTP link)
      - GCF_CDS_File 
      - GCF_Prot_File 
      - GCF_RNA_File
      - GCF_Feature_Table_File
      - GCF_RNA_GBFF_File
      - genometaxID
      - taxonClassif
    
    Usefull links : 
      - https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice
      - https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.ncbi_s_annotation_displayed_on


    Usage :  (Requires three parameters)
    python3  downloadRalevantGCF.py  assembly_summary.txt  outputGCFFolder(= GCF folder in work directory)

    ''')

#------------------------------------------------------------------------------
'''
  Parsing of 'assembly_summary.txt' file.
  GCF.list and "GCF_inexistent.list" creations.
'''
def downloadGCF(summaryfile, outputGCFFolder, outFile) :
    
    # Classification recovery in 'taxonForOrthoMaM_v11.list' file :

    PATH = os_path.abspath(os_path.split(__file__)[0])
    classifFilePath = PATH+os.path.sep+"taxonForOrthoMaM_v11.list"

    with open (classifFilePath,'r') as classif :
        taxonIDtoClassif = {}
        for line in classif :
            parse = line.strip().split("\t")
            taxonID = parse[1]
            taxonClassif = parse[2]
            taxonIDtoClassif[taxonID]=taxonClassif  # making a dictionnary

    extract = outFile.split(".")
    inexistentGCF_file = extract[0]+"_inexistent.list"  # making a path
    GCFlist_File_GZ_extention = []
    GCFlist_File = []

    with open (summaryfile, 'r') as inputFile, open (outFile, 'w') as GCF_list, open (inexistentGCF_file, 'w') as inexistent_GCF_list :
        for line in inputFile :
            if line.startswith("#"): # jumps to the header lines of 'assembly_summary.txt' file
                continue
            genomeInfo = line.strip().split("\t") # recovery of the differents fields of interest
            genomeCategory = genomeInfo[4].strip()
            genometaxID = genomeInfo[5].strip()
            taxonName = genomeInfo[7].strip().replace(" ","_")
            if genometaxID in taxonIDtoClassif.keys(): # condition in case of no classification for a taxon
                taxonClassif = taxonIDtoClassif[genometaxID]
            else :
                taxonClassif = "NA"
            genomeFTPFolder = genomeInfo[19].strip() # FTP link of index files of genomes
            tmp = genomeFTPFolder.split("/")
            GCF_CDS_URL = genomeFTPFolder + "/" + tmp[-1] + "_cds_from_genomic.fna.gz"
            GCF_Prot_URL = genomeFTPFolder + "/" + tmp[-1] + "_protein.faa.gz"
            GCF_RNA_URL = genomeFTPFolder + "/" + tmp[-1] + "_rna.fna.gz"
            GCF_Feature_Table_URL = genomeFTPFolder + "/" + tmp[-1] + "_feature_table.txt.gz"
            GCF_RNA_GBFF_URL = genomeFTPFolder + "/" + tmp[-1] + "_rna.gbff.gz  "

            if (genomeCategory == "representative genome") or (genomeCategory == "reference genome") :
                GCF_CDS_File = tmp[-1] + "_cds_from_genomic.fna"
                GCF_Prot_File = tmp[-1] + "_protein.faa" 
                GCF_RNA_File = tmp[-1] + "_rna.fna"
                GCF_Feature_Table_File = tmp[-1] + "_feature_table.txt"
                GCF_RNA_GBFF_File = tmp[-1] + "_rna.gbff"
                cmd1 = "cd "+outputGCFFolder+"; wget "+GCF_CDS_URL+"; gzip -d "+GCF_CDS_File
                cmd2 = "cd "+outputGCFFolder+"; wget "+GCF_Prot_URL+"; gzip -d "+GCF_Prot_File
                cmd3 = "cd "+outputGCFFolder+"; wget "+GCF_RNA_URL+"; gzip -d "+GCF_RNA_File
                # other style possible : cmd3 = f"cd {outputGCFFolder}; wget {GCF_RNA_URL}; gzip -d {GCF_RNA_File}"
                cmd4 = "cd "+outputGCFFolder+"; wget "+GCF_Feature_Table_URL+"; gzip -d "+GCF_Feature_Table_File
                cmd5 = "cd "+outputGCFFolder+"; wget "+GCF_RNA_GBFF_URL+"; gzip -d "+GCF_RNA_GBFF_File
                os.system(cmd1)
                #rescmd1 = os.popen(cmd1).readlines()
                os.system(cmd2)
                #rescmd2 = os.popen(cmd2).readlines()               
                os.system(cmd3)
                #rescmd3 = os.popen(cmd3).readlines()
                os.system(cmd4)
                #rescmd4 = os.popen(cmd4).readlines()
                os.system(cmd5)
                #rescmd5 = os.popen(cmd5).readlines()

                GCFlist_File_GZ_extention.append(GCF_CDS_File+".gz")
                GCFlist_File_GZ_extention.append(GCF_Prot_File+".gz")
                GCFlist_File_GZ_extention.append(GCF_RNA_File+".gz")
                GCFlist_File_GZ_extention.append(GCF_Feature_Table_File+".gz")
                GCFlist_File_GZ_extention.append(GCF_RNA_GBFF_File+".gz")

                GCFlist_File.append(GCF_CDS_File)
                GCFlist_File.append(GCF_Prot_File)
                GCFlist_File.append(GCF_RNA_File)
                GCFlist_File.append(GCF_Feature_Table_File)
                GCFlist_File.append(GCF_RNA_GBFF_File)


                # Vérification GCF presence
                GCF_Ok = True
                for nameFile in GCFlist_File :
                    cmdCheck = "cd "+outputGCFFolder+";find -name "+nameFile
                    os.system(cmdCheck)
                    rescmdCheck = os.popen(cmdCheck).readlines()
                    if rescmdCheck == '':
                        # add in inexistant_GCF.list
                        inexistent_GCF_list.write(taxonName+"\t"+ genomeFTPFolder+"\t"+ GCF_CDS_File+"\t"+GCF_Prot_File+"\t"+GCF_RNA_File+"\t"+GCF_Feature_Table_File+"\t"+GCF_RNA_GBFF_File+"\t"+genometaxID+"\t"+taxonClassif+"\n")
                        GCF_Ok = False
                        break
                # add in GCF.list
                if GCF_Ok == True :
                    GCF_list.write(taxonName+"\t"+ GCF_CDS_File+"\t"+GCF_Prot_File+"\t"+GCF_RNA_File+"\t"+GCF_Feature_Table_File+"\t"+GCF_RNA_GBFF_File+"\t"+genometaxID+"\t"+taxonClassif+"\n")


#------------------------------------------------------------------------------


if __name__ == "__main__" :
    if len(sys.argv) != 3 : # Requires three parameters
        usage()
        sys.exit()
    summaryfile = sys.argv[1]       # assembly_summary.txt (NCBI)
    outputGCFFolder = sys.argv[2]   # GCF file Path in work drectory
    downloadGCF(summaryfile,outputGCFFolder, "GCF.list")
    