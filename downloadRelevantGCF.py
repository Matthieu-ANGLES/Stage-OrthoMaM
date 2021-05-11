#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import os, sys

#------------------------------------------------------------------------------

'''
Fonction Usage
'''
def usage() :
    print ("\nUsage :")
    print ("For just one taxon")
    print ("./addGenome [listHumanTaxon*] [gene_orthologs(NCBI)] [GCF File(NCBI)] [taxon ID] [taxon name] [Output directory]")
    print("* : Get ListHumanTaxon with addGenome.py")

#------------------------------------------------------------------------------

def downloadGCF(summaryfile, outFile) :
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

            if (genomeCategory == "representative genome") or (genomeCategory == "reference genome") :
                GCF_CDS_File = tmp[-1] + "_cds_from_genomic.fna"
                GCF_Prot_File = tmp[-1] + "_protein.faa"                
                cmd1 = "wget "+GCF_CDS_URL+"; gzip -d "+GCF_CDS_File
                cmd2 = "wget "+GCF_Prot_URL+"; gzip -d "+GCF_Prot_File
                os.system(cmd1)
                os.system(cmd2)

                GCF_liste.write(taxonName+"\t"+ GCF_CDS_File+"\t"+GCF_Prot_File+"\t"+genometaxID+"\n")

if __name__ == "__main__" :
    if len(sys.argv) != 2 :
        usage()
        sys.exit()
    summaryfile = sys.argv[1]       # resume_assembly_summary.txt (NCBI)
    downloadGCF(summaryfile, "GCF.list")
    