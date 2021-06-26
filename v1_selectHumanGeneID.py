#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import sys

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''
    
    Usage : Launch error.
    'selectHumanGeneID.py'
    
    This module parses the 'gene_orthologs.tsv' file and generates a list of orthologous genes (human identifiers : 'humanGeneID.list')
    between the core species list (HomoMusCanis.id = core_species.list)
        
    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id  GeneID  relationship  Other_tax_id  Other_GeneID 
        (following NCBI convention: http://ftp.ncbi.nlm.nih.gov/gene/DATA/ gene-ortholog.gz )
    
    - core_species.list : should contain one taxon id per line, for a core made of Homo Sapiens, Mus musculus and Canis Lupus
    familiaris it will be (human should always be first):
        9606        (Homo sapiens)
        10090       (Mus musculus)
        9615        (Canis Lupus familiaris)

    Usefull links : 
      - https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice
      - https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.ncbi_s_annotation_displayed_on


    Usage :  (Requires four parameters)
    python3 selectHumanGeneID.py gene_orthologs.tsv gene_orthologs.txt outputFolder(= work directory)

    ''')

#------------------------------------------------------------------------------
'''
  Used in selectHumanGeneID() function
  Parsing 'HomoMusCanis.id' (core_species.list) file and generates a list.
'''

def readSpeciesCoreIds (coreTaxonIDFile):
    IDlist = []
    with open (coreTaxonIDFile,'r') as input :
        for line in input :
            parse = line.strip().split("\t")
            geneID = parse[0]
            IDlist.append(geneID)
    return (IDlist)

#------------------------------------------------------------------------------
'''
  Parsing function of the 'gene_orthologs.tsv' file to select the intersection of ortologous human genesID within the core species 
  Generate 'humanGeneID.list' file.
'''

def selectHumanGeneID(orthologFile, coreTaxonIDFile, humanGeneIDList):
    coresTaxonIDList=readSpeciesCoreIds(coreTaxonIDFile)
    orthologsPerTax= {}
    speRef = coresTaxonIDList[0] # human ID
    for i in range(1,len(coresTaxonIDList)): # Homo sapiens, Mus musculus and Canis Lupus familiaris keys.
        orthologsPerTax[coresTaxonIDList[i]]=set() # non-redundant values with a set()

    with open (orthologFile, 'r') as inputFile :
        for line in inputFile :
            ortho_info = line.strip().split("\t")
            spe1= ortho_info[0]
            spe2= ortho_info[3]
            if(spe1==speRef and spe2 in orthologsPerTax.keys()):
                orthologsPerTax[spe2].add(ortho_info[1]) 
    
    coreSet=orthologsPerTax[coresTaxonIDList[1]]
    for spe, setIds in orthologsPerTax.items(): # make intersection of sets values with each core taxon
        coreSet = coreSet.intersection(setIds) 
        print(spe, len(setIds))

    with open (humanGeneIDList, 'a+') as output :
        for ids in coreSet:
            output.write(ids+"\n")
       
#------------------------------------------------------------------------------
'''
  Main
'''

if __name__ == "__main__" :
    if len (sys.argv) != 4: # Requires four parameters
        usage()
        sys.exit()
    else :
        orthologFile=sys.argv[1]            # gene_ortholog (NCBI)
        coreTaxonIDFile=sys.argv[2]         # HomoMusCanis.id (core_species.list)
        humanGeneIDList=sys.argv[3]         # 'humanGeneID.list' (path and name file)

        selectHumanGeneID(orthologFile, coreTaxonIDFile, humanGeneIDList)           