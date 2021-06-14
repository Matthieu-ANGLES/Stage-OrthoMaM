#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import sys

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''

    selectHumanGeneID.py
    
    This module checks the similarity between a CDS reference and 
    its associated protein reference.
    
    Usage :  python3 selectHumanGeneID.py gene_orthologs.tsv output_fasta_folder
    
    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id GeneID  relationship    Other_tax_id    Other_GeneID 
        (following NCBI convention: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene-ortholog.gz )
    
    -core_species.list should contain one taxon id per line, for a core made of Homo Sapiens, Mus musculus and Canis Lupus
    familiaris it will be (human should always be first):
        9606
        10090
        9615
    ''')

#------------------------------------------------------------------------------
'''
  
'''

def readSpeciesCoreIds (core_species_file):
    IDlist = []
    with open (core_species_file,'r') as input :
        for line in input :
            parse = line.strip().split("\t")
            geneID = parse[0]
            IDlist.append(geneID)
    return (IDlist)

#------------------------------------------------------------------------------
'''

'''

def selectHumanGeneID(ortholog_file, core_species_file, res_file):
    #print (core_species_file)
    cores_species_ids=readSpeciesCoreIds(core_species_file)
    orthologs_per_spe= {}
    spe_ref = cores_species_ids[0]
    for i in range(1,len(cores_species_ids)):
        orthologs_per_spe[cores_species_ids[i]]=set()

    with open (ortholog_file, 'r') as inputFile :
        for line in inputFile :
            ortho_info = line.strip().split("\t")
            spe1= ortho_info[0]
            spe2= ortho_info[3]
            if(spe1==spe_ref and spe2 in orthologs_per_spe.keys()):
                #print (spe2)
                orthologs_per_spe[spe2].add(ortho_info[1]) 
    

    # Voir Vincent pour vérif compréhension du code 
    core_set=orthologs_per_spe[cores_species_ids[1]]
    # on prend les val pour un taxon puis on fait l'intersection avec toutes les val du dico ?

    for spe, setIds in orthologs_per_spe.items():
        core_set = core_set.intersection(setIds) 
        print(spe, len(setIds))

    with open (res_file, 'a+') as output :
        for ids in core_set:
            output.write(ids+"\n")
       
#------------------------------------------------------------------------------
'''
  Main
'''

if __name__ == "__main__" :
    if len (sys.argv) != 4:
        usage()
        sys.exit()
    else :
        ortholog_file=sys.argv[1]           # gene_ortholog (NCBI)
        core_species_file=sys.argv[2]       # HomoMusCanis.id
        res_file=sys.argv[3]                # répertoire de travail (destination)

        selectHumanGeneID(ortholog_file, core_species_file, res_file)           