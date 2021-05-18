#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys, os

#------------------------------------------------------------------------------

'''
  Help Message
'''
def usage() :
    
    sys.stderr.write('''
    
    addGenomes.py

    This program build ortholog fasta files.
    
    Usage : python3 addGenomes humanGeneID.list gene_orthologs.tsv GCF_File taxon_ID taxon_name output_fasta_folder
    
    - humanGeneID.list : obtained with the program selectHumanGeneID.py

    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id GeneID  relationship    Other_tax_id    Other_GeneID 
        (following NCBI convention: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene-ortholog.gz )

    - GCF_File : download address in the file assembly_summary.tsv : is supposed to have 22 fields including :
        assembly_accession  refseq_category taxid   organisme_name  ftp_path (to download the GCF_files)
        (following NCBI convention: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt)

    ''')


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
    protID = getProtID(seqName)

    if geneID in geneToLongestCDS.keys():
        name_lg = geneToLongestCDS[geneID]
        if len(seq) > name_lg[1]:
            geneToLongestCDS[geneID] = [seqName,len(seq), protID]
    else :
        geneToLongestCDS[geneID] = [seqName,len(seq), protID]
        
#---------------------------------------------------------------------------------        
'''
Fonction d'accession au GeneID d'une ligne header d'une séquence fasta
'''
def getGeneID (seqName):
    pattern = re.search("(.+)(GeneID:)(\d+)([\,\]])",seqName)
    if pattern == None :
        #print(seqName,"\n")
        return (None)
    geneID = pattern.group(3)

    return (geneID)

#---------------------------------------------------------------------------------        
'''
Fonction d'accession à ACCPROT ! d'une ligne header d'une séquence fasta
'''
def getProtID (seqName):
    pattern = re.search("(.+)(protein_id=)([^\,\]]+)([\,\]])",seqName)
    if pattern == None :
        #print(seqName,"\n")
        return (None)
    protID = pattern.group(3)

    return (protID)

#------------------------------------------------------------------------------
'''
Fonction de création et d'écriture des séquences CDS au format fasta (par gène) dans un répertoire
'''
def addSeq (seqName, seq, taxonToHuman, outputFastaFolder, outputInfoFolder, geneToLongestCDS, taxonName, taxonID) :
    geneID = getGeneID(seqName)
    protID = getProtID(seqName)
    if geneID != None and geneID in taxonToHuman.keys() :
        longestCDSName = geneToLongestCDS[geneID][0]
        if seqName == longestCDSName :
            humanGene = taxonToHuman[geneID]

            fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_NT.fasta"
            infoFileName = outputInfoFolder+os.path.sep+humanGene+".info"
        
            with open (fastaFileName,'a+') as fastaFile :
                newSeqName = ">"+taxonName
                #newSeqName = ">"+taxonName+" [protein_id="+protID+"]"
                fastaFile.write(newSeqName+"\n"+seq+"\n")
            with open (infoFileName,'a+') as infoFile :
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
                infoFile.write(newSeqName+"\n"+seq+"\n")
  
#------------------------------------------------------------------------------
'''
Fonction de création et d'écriture des séquences protéiques au format fasta (par gène) dans un répertoire
'''
def addProtSeq (seqName, seq, outputFastaFolder, protIDtogeneID, taxonName) :

    patternProtID = re.search("(>)([^ ]+)",seqName)
    if patternProtID == None :
        currentProtID = None
    else :
        currentProtID = patternProtID.group(2)

    if (currentProtID != None) and currentProtID in protIDtogeneID.keys() :
        humanGene = protIDtogeneID[currentProtID]
        fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_AA.fasta"
        
        with open (fastaFileName,'a+') as fastaFile :
            newSeqName = ">"+taxonName
            #newSeqName = ">"+taxonName+" [protein_id="+currentProtID+"]"
            fastaFile.write(newSeqName+"\n"+seq+"\n")

#---------------------------------------------------------------------------------        

########################   MAIN   ########################


def addGenome(listIdentOrthoMaM, orthologFile, GCFcdsFile, GCFproteinFile, taxonID, taxonName,outputFastaFolder, outputInfoFolder):
    # Création dictionnaire d'orthologies Humains vs espèce considéré 
    taxonToHuman = {}
    spe_ref = "9606"
    refID = getRefID(listIdentOrthoMaM)

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
                    addSeq(seqName,seq,taxonToHuman,outputFastaFolder, outputInfoFolder,geneToLongestCDS, taxonName, taxonID)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addSeq(seqName,seq,taxonToHuman,outputFastaFolder, outputInfoFolder, geneToLongestCDS, taxonName, taxonID)
    
   
    protIDtogeneID = {}
    for geneID, geneInfo in geneToLongestCDS.items() :
        if geneID in taxonToHuman.keys():
            humanGeneID = taxonToHuman[geneID]
            protID = geneInfo[2]
            protIDtogeneID[protID] = humanGeneID

    with open (GCFproteinFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addProtSeq(seqName,seq,outputFastaFolder,protIDtogeneID, taxonName)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addProtSeq(seqName,seq,outputFastaFolder, protIDtogeneID, taxonName)


if __name__ == "__main__" :
    if len(sys.argv) != 8 :
        usage()
        sys.exit()
    else :
        # Récupération Paramètres :
        listIdentOrthoMaM = sys.argv[1]     # obtenu avec selectHumanGeneID.py
        orthologFile = sys.argv[2]          # fichier gene_orthologs (NCBI)
        GCFcdsFile = sys.argv[3]            # fichier GCF (NCBI)
        GCFproteinFile = sys.argv[4]        # fichier GCF (NCBI)
        taxonID = sys.argv[5]               # ex 9606 Homo Sapiens
        taxonName = sys.argv[6]             # ex Homo_Sapiens
        orthoFolder = sys.argv[7]           # répertoire de destination pour les fichiers fasta générés (ex : orthologFasta/  -> pour orthologFasta/999.fasta)
        addGenome(listIdentOrthoMaM, orthologFile, GCFcdsFile, GCFproteinFile, taxonID, taxonName, orthoFolder)
