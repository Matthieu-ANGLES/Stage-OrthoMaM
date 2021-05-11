#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys

#------------------------------------------------------------------------------

'''
Fonction Usage
'''
def usage() :
    print ("\nUsage :")
    print ("For just one taxon")
    print ("./addGenome [listHumanTaxon*] [gene_orthologs(NCBI)] [GCF File(NCBI)] [taxon ID] [taxon name] [Output directory]")
    print("* : Get ListHumanTaxon with addGenome.py")

#---------------------------------------------------------------------------------        
'''
Fonction qui transforme le fichier obtenue par addGeneID.py en un ensemble de valeurs
'''
def getRefID (refIDFile) :
    refID = set()
    with open (refIDFile,'r') as input :
        for line in input :
            refID.add(line.strip())
    return (refID)

#---------------------------------------------------------------------------------        
'''
Retourne un dictionnaire de header des plus long CDS
'''
def getLongestCDS (GCFFile):
    geneToLongestCDS = {}
    with open (GCFFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    updateLongestCDS(geneToLongestCDS, seq, seqName)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            updateLongestCDS(geneToLongestCDS, seq, seqName)
    
    return (geneToLongestCDS)

#---------------------------------------------------------------------------------        
'''
fonction pour getLongestCDS
'''
def updateLongestCDS(geneToLongestCDS, seq, seqName):
    geneID = getGeneID(seqName)
    if geneID in geneToLongestCDS.keys():
        name_lg = geneToLongestCDS[geneID]
        if len(seq) > name_lg[1]:
            geneToLongestCDS[geneID] = [seqName,len(seq)]
    else :
        geneToLongestCDS[geneID] = [seqName,len(seq)]

#---------------------------------------------------------------------------------        
'''
Fonction d'accession au GeneID d'une ligne header d'une séquence fasta
'''
def getGeneID (seqName):
    pattern = re.search("(.+)(GeneID:)(\d+)([\,\]])",seqName)
    if pattern == None :
        print(seqName,"\n")
        return (None)
    geneID = pattern.group(3)

    return (geneID)

#------------------------------------------------------------------------------
'''
Fonction de création et d'écriture des séquences fasta (par gène) dans un répertoire
'''
def addSeq (seqName, seq, taxonToHuman, orthoFolder, geneToLongestCDS, taxonName, taxonID) :
    geneID = getGeneID(seqName)
    if geneID != None and geneID in taxonToHuman.keys() :
        longestCDSName = geneToLongestCDS[geneID][0]
        if seqName == longestCDSName :
            humanGene = taxonToHuman[geneID]
            fastaFileName = orthoFolder+"/"+humanGene+".fasta"
            infoFileName = orthoFolder+"/"+humanGene+".info"

            # Voir Vincent : Inversion format fasta et info ??
            
            with open (fastaFileName,'a+') as fastaFile :
                newSeqName = ">"+taxonName
                fastaFile.write(newSeqName+"\n"+seq+"\n")
            with open (infoFileName,'a+') as infoFile :
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+seqName[1:]
                infoFile.write(newSeqName+"\n"+seq+"\n")

# voir pour ajouter prot

#---------------------------------------------------------------------------------        

########################   MAIN   ########################


def addGenome(listIdentOrthoMaM, orthologFile, GCFcdsFile, GCFproteinFile, taxonID, taxonName, orthoFolder):
    # Création dictionnaire d'orthologies Humains vs espèce considéré 
    taxonToHuman = {}
    spe_ref = "9606"
    refID = getRefID(listIdentOrthoMaM)
    #msg = "@"+taxonID+"@"
    #print (msg)
    with open (orthologFile, 'r') as inputFile :   # gene_orthologs (NCBI)
        for line in inputFile :
            ortho_info = line.strip().split("\t")
            spe1= ortho_info[0]
            spe2= ortho_info[3]
            idRef = ortho_info[1]
            idTaxon = ortho_info[4]
            if (taxonID ==spe_ref) and (idRef in refID):
                taxonToHuman[idRef] = idRef
            else :
                if(spe1==spe_ref) and (spe2==taxonID) and (idRef in refID):
                    taxonToHuman[idTaxon] = idRef
    
    # Parse Fasta (GCFcdsFile = GCF)
    geneToLongestCDS = getLongestCDS(GCFcdsFile)

    with open (GCFcdsFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addSeq(seqName,seq,taxonToHuman,orthoFolder,geneToLongestCDS, taxonName, taxonID)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addSeq(seqName,seq,taxonToHuman,orthoFolder, geneToLongestCDS, taxonName, taxonID)


if __name__ == "__main__" :
    if len(sys.argv) != 8 :
        usage()
        sys.exit()
    else :
        # Récupération Paramètres :
        listIdentOrthoMaM = sys.argv[1]     # obtenu avec selectHumanGeneID.py
        orthologFile = sys.argv[2]         # fichier gene_orthologs (NCBI)
        GCFcdsFile = sys.argv[3]          # fichier GCF (NCBI)
        GCFproteinFile = sys.argv[4]         # fichier GCF (NCBI)
        taxonID = sys.argv[5]               # ex 9606 Homo Sapiens
        taxonName = sys.argv[6]             # ex Homo_Sapiens
        orthoFolder = sys.argv[7]           # répertoire de destination pour les fichiers fasta générés (ex : orthologFasta/  -> pour orthologFasta/999.fasta)
        addGenome(listIdentOrthoMaM, orthologFile, GCFcdsFile, GCFproteinFile, taxonID, taxonName, orthoFolder)
