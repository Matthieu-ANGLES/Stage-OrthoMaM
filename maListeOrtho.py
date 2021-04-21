#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import os, sys, re, getopt
# os pour lancement cmd bash
# sys pour l'argument en entrée
# re pour regex
# getopt pour la gestions des options et erreurs

import zipfile
# pour accéder à un fichier zip

from colorama import Fore, Style, Back
# pour couleur dans terminal


#----------------------------------------------------------------------
"""
Récupération des geneID du génome humain
"""

def humanGeneID ():
    cmd_Human = "./datasets download genome taxon human --exclude-gene --exclude-protein --exclude-rna"
    print ("Lancement commande bash :\n","./datasets download genome taxon human --exclude-gene --exclude-protein --exclude-rna")
    os.system(cmd_Human)

# tail -n+2 GCF_000001405.39_GRCh38.p13_feature_table.txt | cut -f15,16 |sort -u > gene.list
# tail -n+2 GCF_000001405.39_GRCh38.p13_feature_table.txt | cut -f15,16,20 |sort -u > gene2.list
# tail -n+2 GCF_000001405.39_GRCh38.p13_feature_table.txt | cut -f1,15,16,20 |sort -u > gene3.list

#----------------------------------------------------------------------
"""
Lancement datasets pour recherche des orthologues (en fonction de la liste des gènes récupérée)
"""

def datasetsOrthologues (symbole, taxon):

    print("Download Ortholog Data for : ", Fore.GREEN,symbole,Style.RESET_ALL)
    cmd_Ortho = "./datasets download ortholog symbol " + symbole + " --taxon " + taxon + " --exclude-gene --exclude-protein --exclude-rna --filename " + symbole + "_dataset.zip"
    #print ("Lancement commande bash : ","./datasets download ortholog symbol " +  symbole+ " --taxon " + taxon + " --exclude-gene --exclude-protein --exclude-rna --filename " + symbole + "_dataset.zip")
    os.system(cmd_Ortho)

# datasetsOrthologues ("gapdh", "mouse")

#----------------------------------------------------------------------
"""
LIRE un TSV contenu dans un fichier Zip (ici un zip de téléchargement d'orthologues sur NCBI)
"""

def parseZipFindTSV (fzip): # en parametre = nom du fichier zip ortho (gene en cours)

    #path = "mouse_dataset.zip"

    zfile = zipfile.ZipFile(fzip,'r')

    for filename in zfile.namelist():
        # afficher le contenu
        # print (filename)
        # extraire et réécrire un tsv
        if (filename == "ncbi_dataset/data/data_table.tsv"):
            extract = zfile.read(filename)
            with open ("data.tsv", "w+b") as data:
                data.write(extract)
    os.remove(fzip)


#----------------------------------------------------------------------
'''
Fonction parseOrthoTSV : 
'''

def parseOrthoTSV ():

    # parseZipFindTSV()

    with open ("data.tsv","r") as input :

        with open ("maListeOtho.tsv","a+") as output :
            gene_taxon_precedent = "none"
            for line in input :
                #print(line)
                parse = line.split("\t")
                gene_id = parse[0]
                gene_symbol = parse[1]
                gene_taxon = parse[3]
                if (gene_taxon!=gene_taxon_precedent):
                    if (gene_taxon== "Homo sapiens") or (gene_taxon== "Mus musculus") or (gene_taxon== "Canis lupus familiaris"):
                        output.write(gene_taxon+"\t"+gene_symbol+"\t"+gene_id+"\n")

                gene_id_precedent = gene_id
                gene_symbol_precedent = gene_symbol
                gene_taxon_precedent = gene_taxon
    os.remove("data.tsv")



#----------------------------------------------------------------------
"""
Diagramme de Venn
"""

def diag_venn2 ():
    '''
    set1 = set(['A', 'B', 'C'])
    set2 = set(['A', 'B', 'D'])
    set3 = set(['A', 'E', 'F'])

    venn3([set1, set2, set3], ('Group1', 'Group2', 'Group3'))
    '''
    homo = set(['Gene1', 'Gene2', 'Gene3'])
    mus = set(['Gene1', 'Gene2', 'Gene4'])
    canis = set(['Gene1', 'Gene5', 'Gene3'])

    venn3([homo, mus, canis], ('homo', 'mus', 'canis'))

    

    plt.show()  


#----------------------------------------------------------------------
"""
main
"""
def main(argv):

    # Module getopt pour gestion des args
    # short_options = "hi:o:f:"
    # long_options = ["help", "input=", "output="]
    # rmv = False
    
    try:
        # humanGeneID()

        with open ("gene3.list",'r') as geneList:
            for line in geneList :
                parse = line.split("\t")
                # POUR gene1.list ET gene2.list :
                #gene_symbol = parse[0]
                #gene_attribute = parse[2]
                #gene_product_acc = parse[3]
                
                # POUR gene3.list :
                gene_feature = parse[0]
                gene_symbol = parse[1]
                gene_id = parse[2]
                gene_attribute = parse[3]
                
                if ("-" not in gene_symbol) and (gene_attribute < "32") and (gene_feature == "CDS"):
                    datasetsOrthologues (gene_symbol, "human")
                    fzip = gene_symbol+"_dataset.zip"
                    if os.path.isfile(fzip):
                        # if pour eviter l'arrêt du script lorsque "Error: no valid NCBI gene identifiers, exiting"
                        # exemple   Download Data for :  ACOT1
                        #           Error: no valid NCBI gene identifiers, exiting
                        parseZipFindTSV(fzip)
                        parseOrthoTSV()

    

    # Couple try/except avec Module getopt pour gestion des args
    except getopt.error as err:
        print (str(err))
        # usage()
        #sys.exit(2)

#----------------------------------------------------------------------

if __name__ == "__main__" :
    main(sys.argv[1:])


# gestion des erreurs à revoir 
#   car "Error: no valid NCBI gene identifiers, exiting"