#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import os, sys

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''

    downloadRalevantGCF.py
    
    This module generates the "GCF.list" file used by addGenomes.py
    It has 3 fields :
    tax_name    GCF_file_name   tax_id
    
    Usage :  python3 downloadGCF 
    
    ''')

#------------------------------------------------------------------------------

def downloadGCF(summaryfile, outputGCFFolder, outFile) :
    with open (summaryfile, 'r') as inputFile, open (outFile, 'w') as GCF_liste :
        for line in inputFile :
            if line.startswith("#"):
                continue
            genomeInfo = line.strip().split("\t")
            genomeCategory = genomeInfo[4].strip()
            genometaxID = genomeInfo[5].strip()
            taxonName = genomeInfo[7].strip().replace(" ","_")
            genomeFTPFolder = genomeInfo[19].strip()
            tmp = genomeFTPFolder.split("/")
            GCF_CDS_URL = genomeFTPFolder + "/" + tmp[-1] + "_cds_from_genomic.fna.gz"
            GCF_Prot_URL = genomeFTPFolder + "/" + tmp[-1] + "_protein.faa.gz"
            #---ajout---
            GCF_RNA_URL = genomeFTPFolder + "/" + tmp[-1] + "_rna.fna.gz"
            GCF_Translated_URL = genomeFTPFolder + "/" + tmp[-1] + "_translated_cds.faa.gz"
            GCF_RNA_GBFF_URL = genomeFTPFolder + "/" + tmp[-1] + "_rna.gbff.gz  "
            #-----------  

            if (genomeCategory == "representative genome") or (genomeCategory == "reference genome") :
                GCF_CDS_File = tmp[-1] + "_cds_from_genomic.fna"
                GCF_Prot_File = tmp[-1] + "_protein.faa" 
                #---ajout---
                GCF_RNA_File = tmp[-1] + "_rna.fna"
                GCF_Translated_File = tmp[-1] + "_translated_cds.faa"
                GCF_RNA_GBFF_File = tmp[-1] + "_rna.gbff"
                #-----------                             
                cmd1 = "cd "+outputGCFFolder+"; wget "+GCF_CDS_URL+"; gzip -d "+GCF_CDS_File
                cmd2 = "cd "+outputGCFFolder+"; wget "+GCF_Prot_URL+"; gzip -d "+GCF_Prot_File
                os.system(cmd1)
                os.system(cmd2)
                #---ajout---
                cmd3 = "cd "+outputGCFFolder+"; wget "+GCF_RNA_URL+"; gzip -d "+GCF_RNA_File
                cmd4 = "cd "+outputGCFFolder+"; wget "+GCF_Translated_URL+"; gzip -d "+GCF_Translated_File
                cmd5 = "cd "+outputGCFFolder+"; wget "+GCF_RNA_GBFF_URL+"; gzip -d "+GCF_RNA_GBFF_File
                os.system(cmd3)
                os.system(cmd4)
                os.system(cmd5)
                #----------- 
                
                #GCF_liste.write(taxonName+"\t"+ GCF_CDS_File+"\t"+GCF_Prot_File+"\t"+genometaxID+"\n")
                GCF_liste.write(taxonName+"\t"+ GCF_CDS_File+"\t"+GCF_Prot_File+"\t"+GCF_RNA_File+"\t"+GCF_Translated_File+"\t"+GCF_RNA_GBFF_File+"\t"+genometaxID+"\n")


if __name__ == "__main__" :
    if len(sys.argv) != 2 :
        usage()
        sys.exit()
    summaryfile = sys.argv[1]       # resume_assembly_summary.txt (NCBI)
    downloadGCF(summaryfile, "GCF.list")
    