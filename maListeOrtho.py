#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

from Bio import Entrez
from Bio import SeqIO

#----------------------------------------------------------------------
"""
Recherche orthologues
"""

Entrez.email = "matthieu.angles@hotmail.fr"

# Récupération d'une liste d'ID de séquences résultant d'une requête : fichier XML
ma_req = Entrez.esearch(db="homologene", term="Homo sapiens[Orgn]")
    # requete à améliorer pour obtenir les orthologues 
    # choix db ??

# mon_res devient un dictionnaire
mon_res = Entrez.read(ma_req)

#print("Count : ", mon_res["Count"])     # nbre de séquences
#print("Liste Id : ", mon_res["IdList"]) # liste des id de séquences

#----------------------------------------------------------------------
"""
Récupération identifiants
"""
list_id = mon_res["IdList"]
print("Nombre d'Id : ", len(list_id))
print ("Liste Id : ", list_id)

ma_req.close()
#print("\n")

#----------------------------------------------------------------------
"""
Ecriture des seq au format genBank
Cette parti est aussi à améliorer
"""

# On récupère dans la banque de données protein les séquences de la liste au format Genbank
fic_seq = Entrez.efetch(db="protein", id=list_id, rettype="gb")

# Génère un ensemble de objet seqRecord
mes_seq = SeqIO.parse(fic_seq,"gb")

# Ecriture de l'ensemble dans un fichier output au format fasta
SeqIO.write(mes_seq,"maListeOrtho.fasta", "fasta")
fic_seq.close()