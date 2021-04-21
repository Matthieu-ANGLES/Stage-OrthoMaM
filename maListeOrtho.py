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

    with open ("data.tsv","r") as input, open ("maListeOtho.tsv","a+") as output :
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

def diagVenn ():
    
    with open ("maListeOtho.tsv",'r') as listeOrtho :
        
        # Lecture du fichier tsv généré et créations de Dictionnaires Primaires
        dp_gene_symbole= {}
        dp_gene_id = {}
        
        for line in listeOrtho :
            parse = line.split("\t")
            gene_taxon = parse[0]
            gene_symbol = parse[1].upper()
            gene_id = parse[2]

            if gene_taxon not in dp_gene_symbole.keys() :
                dp_gene_symbole[gene_taxon] = [gene_symbol]
                dp_gene_id[gene_taxon] = [gene_id]
            else :
                if gene_symbol not in dp_gene_symbole.values() :
                    dp_gene_symbole[gene_taxon].append(gene_symbol)
                    dp_gene_id[gene_taxon].append(gene_id)
        #print (dp_gene_symbole)

        for k, v in dp_gene_symbole.items():
            #print (k, v)
            if k == "Homo sapiens" :
                humain = set(dp_gene_symbole[k])
            if k == "Mus musculus" :
                souris = set(dp_gene_symbole[k])
            if k == "Canis lupus familiaris" :
                chien = set(dp_gene_symbole[k])

        venn3([humain, souris, chien], ('Homo sapiens', 'Mus musculus', 'Canis lupus'))
        
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
            nb_gene_total = 0
            nb_gene_symbol_with_attribute_or_not_CDS = 0
            nb_gene_pb_zipfile = 0
            nb_gene_used = 0
            
            for line in geneList :
                nb_gene_total +=1
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
                    nb_gene_symbol_with_attribute_or_not_CDS += 1
                    
                    datasetsOrthologues (gene_symbol, "human")
                    fzip = gene_symbol+"_dataset.zip"
                    if os.path.isfile(fzip) or zipfile.is_zipfile(fzip):
                        # if os.path pour eviter l'arrêt du script lorsque "Error: no valid NCBI gene identifiers, exiting"
                        # exemple   Download Data for :  ACOT1
                        #           Error: no valid NCBI gene identifiers, exiting
                        # if is_zipfile pour le cas où le zip est un fichier txt vide
                        # exemple   zipfile.BadZipFile : File is not a zip File
                        nb_gene_used +=1
                        parseZipFindTSV(fzip)
                        parseOrthoTSV()
                    else :
                        nb_gene_pb_zipfile += 1

            print ("nb gènes recherchés : ",nb_gene_used,"n", "nb gènes de la liste = ",nb_gene_total,"n", "nb gènes avec symboles exclus = ",nb_gene_symbol_with_attribute_or_not_CDS,"n", "nb gènes avec problème de zip = ",nb_gene_pb_zipfile)
            
            diagVenn()

    # Couple try/except avec Module getopt pour gestion des args
    except getopt.error as err:
        print (str(err))
        # usage()
        #sys.exit(2)

#----------------------------------------------------------------------

if __name__ == "__main__" :
    main(sys.argv[1:])