#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import sys, getopt
# sys pour l'argument en entrée
# getopt pour la gestions des options et erreurs (non utilisé)

from colorama import Fore, Style, Back
# pour couleur dans terminal

import matplotlib.pyplot as plt
# pour visu
from matplotlib_venn import venn3
# pour Diagramme de Venn


#----------------------------------------------------------------------
"""
Parcours Fichier gene_orthologs
Recherches orthologies 1:1 pour 5 taxons.
"""

with open ("gene_orthologs", 'r') as inputFile, open ("listes.Resultats",'a+') as output :

    # Listes pour Homo Sapiens
    listeGeneID_HomoSapiens = [] 
    listeGeneID_Homo_Mus = []
    listeOther_gene_id_Homo_Mus = []
    listeGeneID_Homo_Canis = []
    listeOther_gene_id_Homo_Canis = []
    listeGeneID_Homo_Loxodonta = []
    listeOther_gene_id_Homo_Loxodonta = []
    listeGeneID_Homo_Choloepus = []
    listeOther_gene_id_Homo_Choloepus = []

    # Liste pour Mus Musculus
    listeGeneID_MusMusculus = []
    listeGeneID_Mus_Canis = []
    listeOther_gene_id_Mus_Canis = []
    listeGeneID_Mus_Loxodonta = []
    listeOther_gene_id_Mus_Loxodonta = []
    listeGeneID_Mus_Choloepus = []
    listeOther_gene_id_Mus_Choloepus = []  

    # Liste pour Canis Lupus
    listeGeneID_CanisLupus = []
    listeGeneID_Canis_Loxodonta = []
    listeOther_gene_id_Canis_Loxodonta = []
    listeGeneID_Canis_Choloepus = []
    listeOther_gene_id_Canis_Choloepus = []  

    # Liste pour Loxodonta Africana
    listeGeneID_LoxodontaAfricana = []
    listeGeneID_Loxodonta_Choloepus = []
    listeOther_gene_id_Loxodonta_Choloepus = []  

    # Liste pour Choloepus Didactylus
    listeGeneID_CholoepusDidactylus = []

#----------------------------------------------
    # Compteur de lignes
    cptLigne = 0

    # Compteur de taxons 
    dicoCptTaxons = {}

#----------------------------------------------

    '''
    Création d'une listes de gènes ortho Homo et 
    de deux autres listes par couple Homo-autretaxon
    Listes de tailles égales pour récupération other-gene-id 
    en fonction de la position du gene_id dans la liste)
        première liste = gene-id Homo
        deuxième liste = other_gene_id de l'autre taxon
    '''

    for line in inputFile :
        parse = line.strip().split("\t")
        tax_id = parse[0]
        gene_id = parse[1]
        relationship = parse[2]
        other_tax_id = parse[3]
        other_gene_id = parse[4]

        if relationship == "Ortholog" :

        # Compteur lignes
            cptLigne += 1

        # Compteurs taxons et nombre de lignes où il apparait
            if tax_id not in dicoCptTaxons.keys():
                dicoCptTaxons[tax_id] = 1
            else :
                dicoCptTaxons[tax_id] += 1

        # Première partie : Recherche des Othologues 1:1 de Homo vs 4 autres taxons

            if tax_id == "9606" :    # Homo Sapiens
                if gene_id not in listeGeneID_HomoSapiens :
                    listeGeneID_HomoSapiens.append(gene_id)
                if other_tax_id == "10090" :   # Mus Musculus
                    if gene_id not in listeGeneID_Homo_Mus :
                        listeGeneID_Homo_Mus.append(gene_id)
                        listeOther_gene_id_Homo_Mus.append(other_gene_id)
                if other_tax_id == "9615" :   # Canis Lupus
                    if gene_id not in listeGeneID_Homo_Canis :
                        listeGeneID_Homo_Canis.append(gene_id)
                        listeOther_gene_id_Homo_Canis.append(other_gene_id)
                if other_tax_id == "9785" :   # Loxodonta Africana
                    if gene_id not in listeGeneID_Homo_Loxodonta :
                        listeGeneID_Homo_Loxodonta.append(gene_id)
                        listeOther_gene_id_Homo_Loxodonta.append(other_gene_id)
                if other_tax_id == "27675" :   # Choloepus Didactylus
                    if gene_id not in listeGeneID_Homo_Choloepus :
                        listeGeneID_Homo_Choloepus.append(gene_id)
                        listeOther_gene_id_Homo_Choloepus.append(other_gene_id)


    # Deuxième partie : Recherche des Othologues 1:1 de Mus Musculus vs Canis Lupus, Loxodonta Africana et Choloepus Didactylus

            if tax_id == "10090" :    # Mus Musculus
                if gene_id not in listeGeneID_MusMusculus :
                    listeGeneID_MusMusculus.append(gene_id)
                if other_tax_id == "9615" :   # Canis Lupus
                    if gene_id not in listeGeneID_Mus_Canis :
                        listeGeneID_Mus_Canis.append(gene_id)
                        listeOther_gene_id_Mus_Canis.append(other_gene_id)
                if other_tax_id == "9785" :   # Loxodonta Africana
                    if gene_id not in listeGeneID_Mus_Loxodonta :
                        listeGeneID_Mus_Loxodonta.append(gene_id)
                        listeOther_gene_id_Mus_Loxodonta.append(other_gene_id)
                if other_tax_id == "27675" :   # Choloepus Didactylus
                    if gene_id not in listeGeneID_Mus_Choloepus :
                        listeGeneID_Mus_Choloepus.append(gene_id)
                        listeOther_gene_id_Mus_Choloepus.append(other_gene_id)


    # Troisième partie : Recherche des Othologues 1:1 de Canis Lupus vs Loxodonta Africana et Choloepus Didactylus

            if tax_id == "9615" :    # Canis Lupus
                if gene_id not in listeGeneID_CanisLupus:
                    listeGeneID_CanisLupus.append(gene_id)
                if other_tax_id == "9785" :   # Loxodonta Africana
                    if gene_id not in listeGeneID_Canis_Loxodonta :
                        listeGeneID_Canis_Loxodonta.append(gene_id)
                        listeOther_gene_id_Canis_Loxodonta.append(other_gene_id)
                if other_tax_id == "27675" :   # Choloepus Didactylus
                    if gene_id not in listeGeneID_Canis_Choloepus :
                        listeGeneID_Canis_Choloepus.append(gene_id)
                        listeOther_gene_id_Canis_Choloepus.append(other_gene_id)


    # Quatrième partie : Recherche des Othologues 1:1 de Loxodonta Africana vs Choloepus Didactylus

            if tax_id == "10090" :    # Loxodonta Africana
                if gene_id not in listeGeneID_LoxodontaAfricana :
                    listeGeneID_LoxodontaAfricana.append(gene_id)
                if other_tax_id == "27675" :   # Choloepus Didactylus
                    if gene_id not in listeGeneID_Loxodonta_Choloepus :
                        listeGeneID_Loxodonta_Choloepus.append(gene_id)
                        listeOther_gene_id_Loxodonta_Choloepus.append(other_gene_id)


    # Cinquième partie : Recherche des représentants de Choloepus Didactylus

            if tax_id == "10090" :    # Choloepus Didactylus
                if gene_id not in listeGeneID_CholoepusDidactylus :
                    listeGeneID_CholoepusDidactylus.append(gene_id)


#-----------------------------------------------------------------------------------

    # Nombre de lignes dans le fichier
    print(Fore.CYAN,"Nombre de lignes dans le fichier",Style.RESET_ALL)
    print(cptLigne)

    # Nombre de lignes par taxon
    print(Fore.CYAN,"Nombre de lignes par taxons",Style.RESET_ALL)
    total = 0
    for k , v in dicoCptTaxons.items() :
        total += v
        print("Taxon : ",k," =",v," lignes")
    print ("Total : ",total)

    # Print des listes pour Homo Sapiens
    print (Fore.CYAN,"listeGeneID_HomoSapiens",Style.RESET_ALL)
    print (len(listeGeneID_HomoSapiens))
    print ("listeGeneID_Homo_Mus")
    print (len(listeGeneID_Homo_Mus))
    print ("listeOther_gene_id_Homo_Mus")
    print (len(listeOther_gene_id_Homo_Mus))
    print ("listeGeneID_Homo_Canis")
    print (len(listeGeneID_Homo_Canis))
    print ("listeOther_gene_id_Homo_Canis")
    print (len(listeOther_gene_id_Homo_Canis))
    print ("listeGeneID_Homo_Loxodonta")
    print (len(listeGeneID_Homo_Loxodonta))
    print ("listeOther_gene_id_Homo_Loxodonta")
    print (len(listeOther_gene_id_Homo_Loxodonta))
    print ("listeGeneID_Homo_Choloepus")
    print (len(listeGeneID_Homo_Choloepus))
    print ("listeOther_gene_id_Homo_Choloepus")
    print (len(listeOther_gene_id_Homo_Choloepus))
    print ("====================================")

    # Print des listes pour Mus Musculus
    print (Fore.CYAN,"listeGeneID_MusMusculus",Style.RESET_ALL)
    print (len(listeGeneID_MusMusculus))
    print ("listeGeneID_Mus_Canis")
    print (len(listeGeneID_Mus_Canis))
    print ("listeOther_gene_id_Mus_Canis")
    print (len(listeOther_gene_id_Mus_Canis))
    print ("listeGeneID_Mus_Loxodonta")
    print (len(listeGeneID_Mus_Loxodonta))
    print ("listeOther_gene_id_Mus_Loxodonta")
    print (len(listeOther_gene_id_Mus_Loxodonta))
    print ("listeGeneID_Mus_Choloepus")
    print (len(listeGeneID_Mus_Choloepus))
    print ("listeOther_gene_id_Mus_Choloepus")
    print (len(listeOther_gene_id_Mus_Choloepus))
    print ("====================================")

    # Print des listes pour Canis Lupus
    print (Fore.CYAN,"listeGeneID_CanisLupus",Style.RESET_ALL)
    print (len(listeGeneID_CanisLupus))
    print ("listeGeneID_Canis_Loxodonta")
    print (len(listeGeneID_Canis_Loxodonta))
    print ("listeOther_gene_id_Canis_Loxodonta")
    print (len(listeOther_gene_id_Canis_Loxodonta))
    print ("listeGeneID_Canis_Choloepus")
    print (len(listeGeneID_Canis_Choloepus))
    print ("listeOther_gene_id_Canis_Choloepus")
    print (len(listeOther_gene_id_Canis_Choloepus))
    print ("====================================")

    # Print des listes pour Loxodonta Africana
    print (Fore.CYAN,"listeGeneID_LoxodontaAfricana",Style.RESET_ALL)
    print (len(listeGeneID_LoxodontaAfricana))
    print ("listeGeneID_Loxodonta_Choloepus")
    print (len(listeGeneID_Loxodonta_Choloepus))
    print ("listeOther_gene_id_Loxodonta_Choloepus")
    print (len(listeOther_gene_id_Loxodonta_Choloepus))
    print ("====================================")

    # Print de la liste pour Choloepus Didactylus
    print (Fore.CYAN,"listeGeneID_CholoepusDidactylus",Style.RESET_ALL)
    print (len(listeGeneID_CholoepusDidactylus))

#-----------------------------------------------------------------------------------

    # Ecriture des listes pour Homo Sapiens
    output.write("listeGeneID_HomoSapiens"+"\n")
    output.write(str(len(listeGeneID_HomoSapiens))+"\n")
    output.write("listeGeneID_Homo_Mus"+"\n")
    output.write(str(len(listeGeneID_Homo_Mus))+"\n")
    output.write("listeOther_gene_id_Homo_Mus"+"\n")
    output.write(str(len(listeOther_gene_id_Homo_Mus))+"\n")
    output.write("listeGeneID_Homo_Canis"+"\n")
    output.write(str(len(listeGeneID_Homo_Canis))+"\n")
    output.write("listeOther_gene_id_Homo_Canis"+"\n")
    output.write(str(len(listeOther_gene_id_Homo_Canis))+"\n")
    output.write("listeGeneID_Homo_Loxodonta"+"\n")
    output.write(str(len(listeGeneID_Homo_Loxodonta))+"\n")
    output.write("listeOther_gene_id_Homo_Loxodonta"+"\n")
    output.write(str(len(listeOther_gene_id_Homo_Loxodonta))+"\n")
    output.write("listeGeneID_Homo_Choloepus"+"\n")
    output.write(str(len(listeGeneID_Homo_Choloepus))+"\n")
    output.write("listeOther_gene_id_Homo_Choloepus"+"\n")
    output.write(str(len(listeOther_gene_id_Homo_Choloepus))+"\n")
    output.write("===================================="+"\n")

    # Ecriture des listes pour Mus Musculus
    output.write("listeGeneID_MusMusculus"+"\n")
    output.write(str(len(listeGeneID_MusMusculus))+"\n")
    output.write("listeGeneID_Mus_Canis"+"\n")
    output.write(str(len(listeGeneID_Mus_Canis))+"\n")
    output.write("listeOther_gene_id_Mus_Canis"+"\n")
    output.write(str(len(listeOther_gene_id_Mus_Canis))+"\n")
    output.write("listeGeneID_Mus_Loxodonta"+"\n")
    output.write(str(len(listeGeneID_Mus_Loxodonta))+"\n")
    output.write("listeOther_gene_id_Mus_Loxodonta"+"\n")
    output.write(str(len(listeOther_gene_id_Mus_Loxodonta))+"\n")
    output.write("listeGeneID_Mus_Choloepus"+"\n")
    output.write(str(len(listeGeneID_Mus_Choloepus))+"\n")
    output.write("listeOther_gene_id_Mus_Choloepus"+"\n")
    output.write(str(len(listeOther_gene_id_Mus_Choloepus))+"\n")
    output.write("===================================="+"\n")

    # Ecriture des listes pour Canis Lupus
    output.write("listeGeneID_CanisLupus"+"\n")
    output.write(str(len(listeGeneID_CanisLupus))+"\n")
    output.write("listeGeneID_Canis_Loxodonta"+"\n")
    output.write(str(len(listeGeneID_Canis_Loxodonta))+"\n")
    output.write("listeOther_gene_id_Canis_Loxodonta"+"\n")
    output.write(str(len(listeOther_gene_id_Canis_Loxodonta))+"\n")
    output.write("listeGeneID_Canis_Choloepus"+"\n")
    output.write(str(len(listeGeneID_Canis_Choloepus))+"\n")
    output.write("listeOther_gene_id_Canis_Choloepus"+"\n")
    output.write(str(len(listeOther_gene_id_Canis_Choloepus))+"\n")
    output.write("===================================="+"\n")

    # Ecriture des listes pour Loxodonta Africana
    output.write("listeGeneID_LoxodontaAfricana"+"\n")

    output.write("listeGeneID_Loxodonta_Choloepus"+"\n")
    output.write(str(len(listeGeneID_Loxodonta_Choloepus))+"\n")
    output.write("listeOther_gene_id_Loxodonta_Choloepus"+"\n")
    output.write(str(len(listeOther_gene_id_Loxodonta_Choloepus))+"\n")
    output.write("===================================="+"\n")

    # Ecriture de la liste pour Choloepus Didactylus
    output.write("listeGeneID_CholoepusDidactylus"+"\n")
    output.write(str(len(listeGeneID_CholoepusDidactylus))+"\n")


#----------------------------------------------------------------------
    """
    Diagramme de Venn
    """

    Homo = set(listeGeneID_HomoSapiens)
    Mus = set(listeGeneID_Homo_Mus)
    Canis = set(listeGeneID_Homo_Canis)
    Loxodonta = set(listeGeneID_Homo_Loxodonta)
    #Choloepus = set(listeGeneID_Homo_Choloepus)

    venn3([Homo, Mus, Canis], ('Homo sapiens', 'Mus musculus', 'Canis lupus'))
    #venn4([Homo, Mus, Canis, Loxodonta], ('Homo sapiens', 'Mus musculus', 'Canis lupus','Loxodonta Africana'))
    
    
    
    plt.show()  




