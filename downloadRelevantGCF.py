#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import os, sys
from os import path as os_path
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''

    Usage : Launch error.
    'downloadRalevantGCF.py'

    Use the 'assembly_summary_refseq.txt' retrieved at the bottom of the list at this address
      - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt

    
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

    os.system( "cd "+outputGCFFolder)

    with open (summaryfile, 'r') as inputFile, open (outFile, 'w') as GCF_list, open (inexistentGCF_file, 'w') as inexistent_GCF_list :
        for line in inputFile :
            if line.startswith("#"): # jumps to the header lines of 'assembly_summary.txt' file
                continue
            genomeInfo = line.strip().split("\t") # recovery of the differents fields of interest
            genomeCategory = genomeInfo[4].strip()
            genometaxID = genomeInfo[5].strip()

            if (genomeCategory == "representative genome") or (genomeCategory == "reference genome") :
            
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

                urls_to_get = []
                urls_to_get.append(GCF_CDS_URL)
                urls_to_get.append(GCF_Prot_URL)
                urls_to_get.append(GCF_RNA_URL)
                urls_to_get.append(GCF_Feature_Table_URL)
                urls_to_get.append(GCF_RNA_GBFF_URL)


                GCF_Ok = True

                for url in urls_to_get:
                    print("Web site " + url)
                    try:
                        req = urlopen(url)
                    except HTTPError as e:
                        print('The server couldn\'t fulfill the request.', flush=True)
                        print('Error code: ', e.code, flush=True)
                        GCF_Ok = False
                        break
                    except URLError as e:
                        print('We failed to reach a server.', flush=True)
                        print('Reason: ', e.reason, flush=True)
                        GCF_Ok = False
                        break
                    except ValueError as e:
                        print('Invalid ftp url.', flush=True)
                        GCF_Ok = False
                        break
                    else:
                        print("All good", flush=True)
                          
                if GCF_Ok == True :
                    GCF_CDS_File = tmp[-1] + "_cds_from_genomic.fna"
                    GCF_Prot_File = tmp[-1] + "_protein.faa" 
                    GCF_RNA_File = tmp[-1] + "_rna.fna"
                    GCF_Feature_Table_File = tmp[-1] + "_feature_table.txt"
                    GCF_RNA_GBFF_File = tmp[-1] + "_rna.gbff"
                    cmd1 = "cd "+outputGCFFolder+"; wget "+GCF_CDS_URL+"; gzip -d "+GCF_CDS_File
                    cmd2 = "cd "+outputGCFFolder+"; wget "+GCF_Prot_URL+"; gzip -d "+GCF_Prot_File
                    cmd3 = "cd "+outputGCFFolder+"; wget "+GCF_RNA_URL+"; gzip -d "+GCF_RNA_File
                    cmd4 = "cd "+outputGCFFolder+"; wget "+GCF_Feature_Table_URL+"; gzip -d "+GCF_Feature_Table_File
                    cmd5 = "cd "+outputGCFFolder+"; wget "+GCF_RNA_GBFF_URL+"; gzip -d "+GCF_RNA_GBFF_File
                    os.system(cmd1)
                    os.system(cmd2)
                    os.system(cmd3)
                    os.system(cmd4)
                    os.system(cmd5)

                    GCF_list.write(taxonName+"\t"+ GCF_CDS_File+"\t"+GCF_Prot_File+"\t"+GCF_RNA_File+"\t"+GCF_Feature_Table_File+"\t"+GCF_RNA_GBFF_File+"\t"+genometaxID+"\t"+taxonClassif+"\n")
                    print('Taxon with required files: ', taxonName + "\t"+genometaxID+"\t"+taxonClassif+"\n", flush=True)
                else :
                    inexistent_GCF_list.write(taxonName+"\t"+ genomeFTPFolder+"\t"+genometaxID+"\t"+taxonClassif+"\n")
                    print('Taxon without required files: ', taxonName +"\t"+genometaxID+"\t"+taxonClassif+"\n", flush=True)



#------------------------------------------------------------------------------


if __name__ == "__main__" :
    if len(sys.argv) != 3 : # Requires three parameters
        usage()
        sys.exit()
    summaryfile = sys.argv[1]       # assembly_summary.txt (NCBI)
    outputGCFFolder = sys.argv[2]   # GCF file Path in work drectory
    downloadGCF(summaryfile,outputGCFFolder, "GCF.list")
    
