#!/usr/bin/env python3
#_*_ coding:utf-8 _*_


import sys
import matplotlib.pyplot as plt
# pour visu
from matplotlib_venn import venn3
# pour méthode 1 et 2
#pip install matplotlib-venn

#----------------------------------------
'''
Fonction de création du Diagramme de Venn (sur un fichier en paramètre)
'''

def diagrammeDeVenn (inputFile):
    '''
    ########################### EXEMPLE ###########################
    homo = set(['Gene1', 'Gene2', 'Gene3'])
    mus = set(['Gene1', 'Gene2', 'Gene4'])
    canis = set(['Gene1', 'Gene5', 'Gene3'])

    venn3([homo, mus, canis], ('homo', 'mus', 'canis'))
    ###############################################################
    '''

    with open (inputFile,'r') as listeOrtho :
        
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


#Lancement fonction :
diagrammeDeVenn(sys.argv[1])