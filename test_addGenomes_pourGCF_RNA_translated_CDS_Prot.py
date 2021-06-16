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
    geneName = getGeneName(seqName) #utile pour lien geneID vs geneName (rna)

    if geneID in geneToLongestCDS.keys():
        name_lg = geneToLongestCDS[geneID]
        if len(seq) > name_lg[1]:
            geneToLongestCDS[geneID] = [seqName,len(seq), protID, geneName]
    else :
        geneToLongestCDS[geneID] = [seqName,len(seq), protID, geneName]
        

#--------------------------------------------------------------------------------- 
'''       
'''
#Retourne un dictionnaire de header des plus long ARN  => NE RECUPERE PAS LE BON RNA
'''
def getLongestRNA (GCFFile):
    geneToLongestRNA = {}
    with open (GCFFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            #print (line)
            if line.startswith(">") :
                if first == False :
                    #print (seqName)
                    updateLongestRNA(geneToLongestRNA, seq, seqName)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            updateLongestRNA(geneToLongestRNA, seq, seqName)
    
    return (geneToLongestRNA)


#---------------------------------------------------------------------------------        
'''
#fonction pour getLongestRNA  => => NE RECUPERE PAS LE BON RNA
'''
def updateLongestRNA(geneToLongestRNA, seq, seqName):
    #print (seqName)
    patternSeqID = re.search("(>)(\w+.\d)([^ ]+)",seqName)
    if patternSeqID == None :
        currentSeqID = None
    else :
        currentSeqID = patternSeqID.group(2)+patternSeqID.group(3)
    #print (seqName)
    geneName = getGeneNamebis(seqName)
    #print (geneName)

    if geneName in geneToLongestRNA.keys():
        name_lg = geneToLongestRNA[geneName]
        if len(seq) > name_lg[1]:
            geneToLongestRNA[geneName] = [seqName,len(seq),currentSeqID]
    else :
        geneToLongestRNA[geneName] = [seqName,len(seq),currentSeqID]
        
'''
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
Fonction d'accession au GeneName d'une ligne header d'une séquence fasta
'''
def getGeneName (seqName):
    pattern = re.search("(.+)(gene=)(\w+)([\,\]])",seqName)
    if pattern == None :
        #print(seqName,"\n")
        return (None)
    geneName = pattern.group(3)

    return (geneName)

#---------------------------------------------------------------------------------        
'''
Fonction d'accession au GeneName d'une ligne header d'une séquence fasta
'''
def getGeneNamebis (seqName):
    pattern = re.search("(.+)(\()(\w+)([\)])",seqName)
    if pattern == None :
        #print(seqName,"\n")
        return (None)
    geneName = pattern.group(3)

    return (geneName)

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

#---------------------------------------------------------------------------------        
'''
Fonction d'accession à ACCRNA ! d'une ligne header d'une séquence fasta
'''
def getRNAacc(seqName):
    pattern = re.search("(>)(\w+)(\.\d+)([^.+])",seqName)
    if pattern == None :
        return (None)
    RNAacc = pattern.group(2)+pattern.group(3)

    return (RNAacc)

#------------------------------------------------------------------------------
'''
Fonction de création et d'écriture des séquences CDS au format fasta (par gène) dans un répertoire
'''
def addSeq (seqName, seq, taxonToHuman, outputFastaFolder, outputInfoFolder, geneToLongestCDS, taxonName, taxonID) :
    geneID = getGeneID(seqName)
    if geneID != None and geneID in taxonToHuman.keys() :
        longestCDSName = geneToLongestCDS[geneID][0]
        if seqName == longestCDSName :
            humanGene = taxonToHuman[geneID]

            fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_NT.fasta"
            infoFileName = outputInfoFolder+os.path.sep+humanGene+".info"
        
            with open (fastaFileName,'a+') as fastaFile :
                #newSeqName = ">"+taxonName
                #newSeqName = ">"+taxonName+" [protein_id="+protID+"]"
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]

                fastaFile.write(newSeqName+"\n"+seq+"\n")
            with open (infoFileName,'a+') as infoFile :
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
                infoFile.write(newSeqName+"\n"+seq+"\n")
  
#------------------------------------------------------------------------------
'''
Fonction de création et d'écriture des séquences protéiques au format fasta (par gène) dans un répertoire
'''
def addProtSeq (seqName, seq, outputFastaFolder, protIDtogeneID, taxonName, taxonID) :

    patternProtID = re.search("(>)([^ ]+)",seqName)
    if patternProtID == None :
        currentProtID = None
    else :
        currentProtID = patternProtID.group(2)

    if (currentProtID != None) and currentProtID in protIDtogeneID.keys() :
        humanGene = protIDtogeneID[currentProtID]
        fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_AA.fasta"
        
        with open (fastaFileName,'a+') as fastaFile :
            #newSeqName = ">"+taxonName
            #newSeqName = ">"+taxonName+" [protein_id="+currentProtID+"]"
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
            fastaFile.write(newSeqName+"\n"+seq+"\n")

#------------------------------------------------------------------------------
''' ECRIT AVEC LA METHODE QUI NE RECUPERE PAS LE BON RNA
'''
#Fonction de création et d'écriture des séquences rna au format fasta (par gène) dans un répertoire
'''
def addRNASeq (seqName, seq, outputFastaFolder, geneIDtoNameGene, geneToLongestRNA, taxonName, taxonID) : 

    patternRNAGeneName = re.search("(.+)(\()(\w+)([\)])",seqName)
    if patternRNAGeneName == None :
        geneName = None
    else :
        geneName = patternRNAGeneName.group(3)
    #print (seqName)
    #print(geneIDtoNameGene.values())
    if (geneName != None) and geneName in geneIDtoNameGene.values() :
        longestRNAName = geneToLongestRNA[geneName][0]
        #print(seqName)
        #print(longestRNAName)
        if seqName == longestRNAName :
            key = [k for(k, v) in geneIDtoNameGene.items() if v == geneName]
            humanGene = key[0]
            fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_RNA.fasta"
                
            with open (fastaFileName,'a+') as fastaFile :
                #newSeqName = ">"+taxonName
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
                fastaFile.write(newSeqName+"\n"+seq+"\n")
'''


#------------------------------------------------------------------------------
'''
Fonction de récupération et d'écriture des séquences protéiques au format fasta (par gène) dans un répertoire
'''
def addRNASeq (seqName, seq, protIDtogeneID, protIDtoRNAacc, outputFastaFolder, outputInfoFolder, taxonName, taxonID) : 

    currentRNAacc = getRNAacc(seqName)
    key_list = [k for (k, v) in protIDtoRNAacc.items() if v[0] == currentRNAacc]
    
    if (len(key_list) == 0):
        pass
    else:
        currentProtID = key_list[0]
        curentValues = protIDtoRNAacc[currentProtID]
        #print (curentValues)
        markerCDS = curentValues[1].split("..")
        if (len(markerCDS) != 2) or (markerCDS[0].isdigit()== False)  or (markerCDS[1].isdigit()== False) :
            print("Bornes CDS complexes pour : ",currentRNAacc,markerCDS)
            pass
        else :
            startMarker = int(markerCDS[0])-1
            endMarker = int(markerCDS[1])
    
            if currentProtID != None and currentProtID in protIDtogeneID.keys():
                humanGene = protIDtogeneID[currentProtID]
                if currentProtID != "" and currentProtID in protIDtoRNAacc.keys() :
                    fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_RNA.fasta"
                    fastaFileNameforCkeck = outputFastaFolder+os.path.sep+humanGene+"_RNAtoCDS.fasta"
                    infoFileName = outputInfoFolder+os.path.sep+humanGene+"_RNAtoCDS.info"

                    #print (currentRNAacc)
                    #print (currentProtID)
                    #print (humanGene)
                    #print("--------")

                    # Construction RNA.fasta
                    with open (fastaFileName,'a+') as fastaFile :
                        #newSeqName = ">"+taxonName
                        newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+" [protein_id="+currentProtID+"] "+seqName[1:]
                        fastaFile.write(newSeqName+"\n"+seq+"\n")

                    # Construction CDS à partir du fichier GBFF et de la sequence RNA
                    with open (fastaFileNameforCkeck,'a+') as otherFastaFile :
                        newSeq = seq[startMarker:endMarker]
                        newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" (CDS maked with RNA gbff informations of"+currentRNAacc+")"
                        otherFastaFile.write(newSeqName+"\n"+newSeq+"\n")

                    with open (infoFileName,'a+') as infoFile :
                        #newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
                        newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" (CDS maked with RNA gbff informations of"+currentRNAacc+")"
                        infoFile.write(newSeqName+"\n"+seq+"\n")
  



#------------------------------------------------------------------------------
'''
Fonction de création et d'écriture des séquences protéiques au format fasta (par gène) dans un répertoire
'''
def addTranslatedSeq (seqName, seq, outputFastaFolder, protIDtogeneID, taxonName, taxonID) :

    protID = getProtID(seqName)

    if protID == None :
        currentProtID = None
    else :
        currentProtID = protID

    if (currentProtID != None) and currentProtID in protIDtogeneID.keys() :
        humanGene = protIDtogeneID[currentProtID]
        fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_CDSTranslated.fasta"
        
        with open (fastaFileName,'a+') as fastaFile :
            #newSeqName = ">"+taxonName
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
            fastaFile.write(newSeqName+"\n"+seq+"\n")

#---------------------------------------------------------------------------------        
       

########################   MAIN   ########################


def addGenome(listIdentOrthoMaM, orthologFile, GCFcdsFile, GCFproteinFile, GCFRNAFile, GCFTranslatedFile, GCFRNA_GBFFFile, taxonID, taxonName, outputFastaFolder, outputInfoFolder):
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
    #print (geneToLongestCDS)
        
    # For CDS
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
    geneIDtoNameGene = {}
    for geneID, geneInfo in geneToLongestCDS.items() :
        if geneID in taxonToHuman.keys():
            humanGeneID = taxonToHuman[geneID]
            protID = geneInfo[2]
            protIDtogeneID[protID] = humanGeneID
            
            '''
            # ajout pour récup geneName du GCF RNA => fonctionne avec la version qui ne récupère pas le bon RNA
            geneName = geneInfo[3]
            geneIDtoNameGene[humanGeneID] = geneName
            '''

    # For Prot
    with open (GCFproteinFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addProtSeq(seqName,seq,outputFastaFolder,protIDtogeneID, taxonName, taxonID)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addProtSeq(seqName,seq,outputFastaFolder, protIDtogeneID, taxonName, taxonID)

    '''
    # Cette méthode ne récupère pas le bon RNA !!

    # For RNA
    geneToLongestRNA = getLongestRNA(GCFRNAFile)
    #print (geneToLongestRNA)
    with open (GCFRNAFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            #print (line)
            if line.startswith(">") :
                if first == False :
                    #print (seqName)
                    addRNASeq(seqName, seq, outputFastaFolder, geneIDtoNameGene,geneToLongestRNA, taxonName, taxonID)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addRNASeq(seqName, seq, outputFastaFolder, geneIDtoNameGene,geneToLongestRNA, taxonName, taxonID)
    '''

    # Pour RNA récupération RNA grace au prot ID du CDS initial confronté au fichier feature table
    with open (GCFRNA_GBFFFile,'r') as gbffFile :
        protIDtoRNAacc = {}
        for line in gbffFile:
            #if line.startswith("LOCUS"):
                #parse = line.strip().split("       ")     # 7 espaces !
            if line.startswith("VERSION"):
                parse = line.strip().split("     ")     # 5 espaces !
                RNAacc = parse[1]
                #print (RNAacc)
                #print (parse[1], parse[2])

            if line.startswith("     "):                  # 5 espaces !
                parse = line.strip().split("     ")
                if parse[0] == "CDS" :
                    CDS = parse[2].strip()
                    #print ("CDS : ",CDS)
            if line.startswith("                     "):  # 21 espaces !
                parse = line.strip().split("                     ")
                parse = parse[0]
                if parse.startswith("/protein_id="):
                    subline = parse.split("\"")
                    protID = subline[1]
                    #print (protID)
                    #print ("------------------")
                    if protID not in protIDtoRNAacc.keys():
                        protIDtoRNAacc[protID] = [RNAacc,CDS]

    with open (GCFRNAFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            #print (line)
            if line.startswith(">") :
                if first == False :
                    #print (seqName)
                    addRNASeq(seqName, seq, protIDtogeneID, protIDtoRNAacc, outputFastaFolder, outputInfoFolder taxonName, taxonID)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addRNASeq(seqName, seq, protIDtogeneID, protIDtoRNAacc, outputFastaFolder, outputInfoFolder taxonName, taxonID)




    # For Translated
    with open (GCFTranslatedFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addTranslatedSeq(seqName,seq,outputFastaFolder,protIDtogeneID, taxonName, taxonID)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addTranslatedSeq(seqName,seq,outputFastaFolder, protIDtogeneID, taxonName, taxonID)

if __name__ == "__main__" :
    if len(sys.argv) != 10 :
        usage()
        sys.exit()
    else :
        # Récupération Paramètres :
        listIdentOrthoMaM = sys.argv[1]     # obtenu avec selectHumanGeneID.py
        orthologFile = sys.argv[2]          # fichier gene_orthologs (NCBI)
        GCFcdsFile = sys.argv[3]            # fichier GCF (NCBI)
        GCFproteinFile = sys.argv[4]        # fichier GCF (NCBI)
        GCFRNAFile = sys.argv[5] 
        GCFTranslatedFile = sys.argv[6] 
        GCFRNA_GBFFFile = sys.argv[7]
        taxonID = sys.argv[8]               # ex 9606 Homo Sapiens
        taxonName = sys.argv[9]             # ex Homo_Sapiens
        orthoFolder = sys.argv[10]           # répertoire de destination pour les fichiers fasta générés (ex : orthologFasta/  -> pour orthologFasta/999.fasta)
        
        outputFastaFolder = orthoFolder+os.path.sep+"FASTA"
        outputInfoFolder = orthoFolder+os.path.sep+"INFO"
        os.system("mkdir "+outputFastaFolder)
        os.system("mkdir "+outputInfoFolder)
        
        addGenome(listIdentOrthoMaM, orthologFile, GCFcdsFile, GCFproteinFile, GCFRNAFile, GCFTranslatedFile, GCFRNA_GBFFFile, taxonID, taxonName, outputFastaFolder, outputInfoFolder)
