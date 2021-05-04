#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys


def usage() :
    print ("Ajouter fichier pour analyse (ex: res.canis")
    # améliorer

if len(sys.argv) != 7 :
    usage()
    sys.exit()

# ajouter exemple de commande
#------------------------------------------------------------------------------


listIdentOrthoMaM = sys.argv[1]
ortholog_file = sys.argv[2] #
genomeFile = sys.argv[3]               # fasta courant
taxonID = sys.argv[4]
taxonName = sys.argv[5]
orthoFolder = sys.argv[6]           # pour fichiers fastas

# recupérer list ID dans un ensemble
def readGeneIds (geneIdsFile):
    IDlist = set()
    with open (core_species_file,'r') as input :
        for line in input :
            parse = line.strip().split("\t")
            geneID = parse[0]
            IDlist.add(geneID)
    return (IDlist)

def addSeq (seqName, seq, taxonToHuman, orthoFolder) :
    geneID = getGeneID(seqName)
    if geneID != None and geneID in taxonToHuman.keys() :
        humanGene = taxonToHuman[geneID]
        if geneID in taxonToHuman.keys() :
            humanGene = taxonToHuman[geneID]
            fileName = orthoFolder+"/"+humanGene+".fasta"
            with open (fileName,'a+') as fastaFile :
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+seqName[1:]
                fastaFile.write(newSeqName+"\n"+seq+"\n")
            
def getGeneID (seqName):
    geneID = None
    info = seqName.split("[db_xref=GeneID:")
    if len(info)==2 :
        info2 = info[1].split("]")
        #geneID = re.search("dbxref="]  FAIRE REGEX !!
        geneID = info2[0]
    return (geneID)


# Création dictionnaire d'orthologies Humains vs espèce considéré
taxonToHuman = {}
spe_ref = "9606"

with open (ortholog_file, 'r') as inputFile :
    for line in inputFile :
        ortho_info = line.strip().split("\t")
        spe1= ortho_info[0]
        spe2= ortho_info[3]
        idRef = ortho_info[1]
        idTaxon = ortho_info[4]
        if(spe1==spe_ref and spe2==taxonID):
            taxonToHuman[idTaxon] = idRef

# Parse Fasta
with open (genomeFile,'r') as fastaFile :
    seq = ""
    seqName = ""
    first = True
    for line in fastaFile :
        if line.startswith(">") :
            if first == False :
                addSeq(seqName,seq,taxonToHuman,orthoFolder)
            seqName = line.strip()
            seq = ""
            first = False
        else :
            seq += line.strip()
    if first == False :
        addSeq(seqName,seq,taxonToHuman,orthoFolder)
