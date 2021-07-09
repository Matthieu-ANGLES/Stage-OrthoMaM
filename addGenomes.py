#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import re, sys, os

#------------------------------------------------------------------------------

'''
  Help Message
'''
def usage() :
    
    sys.stderr.write('''

    Usage : Launch error.
    'addGenomes.py'

    This program build fasta files of the orthologous genes selected in the 'humanGeneID.list'.
    From the gene identifier, the module browses all the referenced CDS and retrived the largest one
    as the reference one.
    Thanks to the identifiers of these CDS, tu module also recovers the corresponding proteins and
    generates new CDS sequences from the RNA sequences whose bounds an others informations are recovered 
    in the .bgff and features tables files.
    (see the usefull links for more informations)

    - humanGeneID.list : obtained with the program selectHumanGeneID.py

    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id  GeneID  relationship  Other_tax_id  Other_GeneID 
        (following NCBI convention: http://ftp.ncbi.nlm.nih.gov/gene/DATA/ gene-ortholog.gz )
    
    - core_species.list : should contain one taxon id per line, for a core made of Homo Sapiens, Mus musculus and Canis Lupus
    familiaris it will be (human should always be first):
        9606        (Homo sapiens)
        10090       (Mus musculus)
        9615        (Canis Lupus familiaris)

    Usefull links : 
      - https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice
      - https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.ncbi_s_annotation_displayed_on


    Usage :  (Requires twelve parameters)
    python3 addGenomes.py humanGeneIDList orthologFile GCF_CDS_File GCF_Prot_File GCF_RNA_File GCF_Translate_File GCF_RNA_GBFF_File taxID taxName taxClassif outputFolder

    ''')


#---------------------------------------------------------------------------------        
'''
Function that transforms the 'humanGeneID.list' file into a set of values
'''
def getRefID (humanGeneIDList) :
    refID = set()
    with open (humanGeneIDList,'r') as input :
        for line in input :
            refID.add(line.strip())
    return (refID)

#---------------------------------------------------------------------------------        
'''
Function that return a header dictionary of the longest CDS
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
Subfunction for getLongestCDS() function
'''
def updateLongestCDS(geneToLongestCDS, seq, seqName):
    geneID = getGeneID(seqName)
    protID = getProtID(seqName)
    geneName = getGeneName(seqName) # useful for the geneID vs geneName link

    if geneID in geneToLongestCDS.keys():
        CDSinformations = geneToLongestCDS[geneID]
        if len(seq) > CDSinformations[1]: # = len(seq)
            geneToLongestCDS[geneID] = [seqName, len(seq), protID, geneName]
    else :
        geneToLongestCDS[geneID] = [seqName, len(seq), protID, geneName]
        
#---------------------------------------------------------------------------------        
'''
REGEX function to access the geneID of a header line of the fasta file.
'''
def getGeneID (seqName):
    pattern = re.search("(.+)(GeneID:)(\d+)([\,\]])",seqName)
    if pattern == None :
        return (None)
    geneID = pattern.group(3)

    return (geneID)

#---------------------------------------------------------------------------------        
'''
REGEX function to access the geneName of a header line of the fasta file.
Used in updateLongestCDS() function.
'''
def getGeneName (seqName):
    pattern = re.search("(.+)(gene=)(\w+)([\,\]])",seqName)
    if pattern == None :
        return (None)
    geneName = pattern.group(3)

    return (geneName)

#---------------------------------------------------------------------------------        
'''
REGEX function to access the protein accession number of a header line of the fasta file.
Used in updateLongestCDS() function.
'''
def getProtID (seqName):
    pattern = re.search("(.+)(protein_id=)([^\,\]]+)([\,\]])",seqName)
    if pattern == None :
        return (None)
    protID = pattern.group(3)

    return (protID)

#---------------------------------------------------------------------------------        
'''
REGEX function to access the RNA accession number of a header line of the fasta file.
Used in addRNASeq() function.
'''
def getRNAacc(seqName):
    pattern = re.search("(>)(\w+)(\.\d+)([^.+])",seqName)
    if pattern == None :
        return (None)
    RNAacc = pattern.group(2)+pattern.group(3)

    return (RNAacc)

#------------------------------------------------------------------------------
'''
Fasta files creation (by human geneID) and writing of CDS sequences in a dedicated directory
'''
def addSeq (seqName, seq, taxonToHuman, outputFastaFolder, outputInfoFolder, geneToLongestCDS, taxonName, taxonID, taxClassif) :
    geneID = getGeneID(seqName)
    if geneID != None and geneID in taxonToHuman.keys() :
        longestCDSName = geneToLongestCDS[geneID][0]
        if seqName == longestCDSName :
            humanGene = taxonToHuman[geneID]

            fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_NT.fasta"
            infoFileName = outputInfoFolder+os.path.sep+humanGene+"_NT.info"
        
            with open (fastaFileName,'a+') as fastaFile :
                #newSeqName = ">"+taxonName
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]+" [taxonClassif="+taxClassif+"]"
                fastaFile.write(newSeqName+"\n"+seq+"\n")
                
            with open (infoFileName,'a+') as infoFile :
                newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]+" [taxonClassif="+taxClassif+"]"
                infoFile.write(newSeqName+"\n"+seq+"\n")
  
#------------------------------------------------------------------------------
'''
Fasta files creation (by human geneID) and writing of Proteins sequences in a dedicated directory
'''
def addProtSeq (seqName, seq, outputFastaFolder, protIDtogeneID, taxonName, taxonID, taxClassif) :

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
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]+" [taxonClassif="+taxClassif+"]"
            fastaFile.write(newSeqName+"\n"+seq+"\n")


#------------------------------------------------------------------------------
'''
Fasta files creation (by human geneID) and writing of RNA sequences in a dedicated directory
'''
def addRNASeq (seqName, seq, protIDtogeneID, protIDtoRNAacc, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif) : 

    currentRNAacc = getRNAacc(seqName)
    key_list = [k for (k, v) in protIDtoRNAacc.items() if v[0] == currentRNAacc]
    
    if (len(key_list) == 0):
        pass
    else:
        currentProtID = key_list[0]
        curentValues = protIDtoRNAacc[currentProtID] # = start and stop end points CDS
        if currentProtID != None and currentProtID in protIDtogeneID.keys():
            humanGene = protIDtogeneID[currentProtID]
            if currentProtID != "" and currentProtID in protIDtoRNAacc.keys() :
                markerCDS = curentValues[1].split("..")
                if (len(markerCDS) != 2) or (markerCDS[0].isdigit()== False)  or (markerCDS[1].isdigit()== False) :
                    print("Bornes CDS complexes pour : ",currentRNAacc,markerCDS)
                    #pass

                    # Cas du 'join' (exemple : markerCDS = ['join(4','26,35','1986)'])
                    if (len(markerCDS) == 3) :
                        joinStart = markerCDS[0].split("(")         # = 'join(4' sur l'exemple
                        joinStart = int(joinStart[1])-1             # = 4 sur l'exemple
                        extractMiddle = markerCDS[1].split(",")     # = 26,35 sur l'exemple
                        joinStartToMiddle = int(extractMiddle[0])   # = 26
                        joinMiddleToEnd = int(extractMiddle[1])-1   # = 35
                        joinStop = markerCDS[2].split(")")          # = '1986)' sur l'exemple
                        joinStop = int(joinStop[0])                 # = 1986

                        newSeq = seq[joinStart:joinStartToMiddle]+seq[joinMiddleToEnd:joinStop]
                        createRNAFastaFiles (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)
                        print (currentRNAacc, " ajoutée aux séquences fasta")

                    # Cas du '<' (example : markerCDS = ['<1', '1986'])
                    # ATTENTION
                    elif (markerCDS[0].isdigit()== False) and (markerCDS[1].isdigit()== True):
                        extractMarker = markerCDS[0].split("<")     # ['', '1']
                        startMarker = int(extractMarker[1])-1       # = '1' (-1)
                        endMarker = int(markerCDS[1])

                        newSeq = seq[startMarker:endMarker]
                        excludedPartialSequenceRNAtoCDS (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)
                        print (currentRNAacc, " (RNAtoCDS) ajoutée dans dossier INFO/")

                    # Cas du '>' (example : markerCDS = ['1', '>1986'])
                    elif (markerCDS[0].isdigit()== True) and (markerCDS[1].isdigit()== False):
                        extractMarker = markerCDS[1].split(">")     # ['', '1986']
                        startMarker = int(markerCDS[0])-1           
                        endMarker = int(extractMarker[1])           # = '1986'

                        newSeq = seq[startMarker:endMarker]
                        excludedPartialSequenceRNAtoCDS (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)
                        print (currentRNAacc, " (RNAtoCDS) ajoutée dans dossier INFO/")

                    # Cas du '<' et '>' (example : markerCDS = ['<1', '>1986'])
                    elif (markerCDS[0].isdigit()== False) and (markerCDS[1].isdigit()== False):
                        extractMarker = markerCDS[0].split("<")     # ['', '1']
                        startMarker = int(extractMarker[1])-1       # = '1' (-1)
                        extractMarker = markerCDS[1].split(">")     # ['', '1986']
                        endMarker = int(extractMarker[1])           # = '1986'

                        newSeq = seq[startMarker:endMarker]
                        excludedPartialSequenceRNAtoCDS (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)
                        print (currentRNAacc, " (RNAtoCDS) ajoutée dans dossier INFO/")
                        
                else : # Cas 'normal' (exemple : markerCDS = ['86','1986']
                    startMarker = int(markerCDS[0])-1
                    endMarker = int(markerCDS[1])

                    newSeq = seq[startMarker:endMarker]
                    createRNAFastaFiles (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)

#------------------------------------------------------------------------------
'''
  RNA and RNAtoCDS writing function in FASTA and INFO file
  Used in addRNASeq() function.
'''

def createRNAFastaFiles (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif):

    fastaFileName = outputFastaFolder+os.path.sep+humanGene+"_RNA.fasta"
    fastaFileNameforCkeck = outputFastaFolder+os.path.sep+humanGene+"_RNAtoCDS.fasta"
    infoFileName = outputInfoFolder+os.path.sep+humanGene+"_RNAtoCDS.info"

    # RNAaccToGeneFeatureInformation[RNAacc]=[geneID, chromosome, strand, chromStart, chromStop, CCDS]
    geneID = RNAaccToGeneFeatureInformation[currentRNAacc][0]
    chromosome = RNAaccToGeneFeatureInformation[currentRNAacc][1]
    strand = RNAaccToGeneFeatureInformation[currentRNAacc][2]
    chromStart = RNAaccToGeneFeatureInformation[currentRNAacc][3]
    chromStop = RNAaccToGeneFeatureInformation[currentRNAacc][4]
    if len(RNAaccToGeneFeatureInformation[currentRNAacc]) == 6 : # not all sequences have a CCDS accession
        CCDS = RNAaccToGeneFeatureInformation[currentRNAacc][5]

    # RNA.fasta construction
    with open (fastaFileName,'a+') as fastaFile :
        #newSeqName = ">"+taxonName
        newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+" [protein_id="+currentProtID+"] "+seqName[1:]+" [taxonClassif="+taxClassif+"]"
        fastaFile.write(newSeqName+"\n"+seq+"\n")

    # Construction CDS à partir du fichier GBFF et de la sequence RNA
    with open (fastaFileNameforCkeck,'a+') as otherFastaFile :
        #newSeqName = ">"+taxonName
        newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" (CDS maked with RNA gbff informations from "+currentRNAacc+")"+" [taxonClassif="+taxClassif+"]"
        otherFastaFile.write(newSeqName+"\n"+newSeq+"\n")

    with open (infoFileName,'a+') as infoFile :
        if len(RNAaccToGeneFeatureInformation[currentRNAacc]) == 6 :
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" [RNAacc="+currentRNAacc+"]"+" [geneID="+geneID+"]"+" [chromosome="+chromosome+"]"+" [strand="+strand+"]"+" [start="+chromStart+"]"+" [stop="+chromStop+"]"+" [CCDS="+CCDS+"]"+" [taxonClassif="+taxClassif+"]"
        else :
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" [RNAacc="+currentRNAacc+"]"+" [geneID="+geneID+"]"+" [chromosome="+chromosome+"]"+" [strand="+strand+"]"+" [start="+chromStart+"]"+" [stop="+chromStop+"]"+" [CCDS="+"NA"+"]"+" [taxonClassif="+taxClassif+"]"
        
        infoFile.write(newSeqName+"\n")

#------------------------------------------------------------------------------
'''
  Writting RNA and RNAtoCDS in fastaFile or in infoFile 
  Used in addRNASeq function
'''

def excludedPartialSequenceRNAtoCDS (seqName, seq, newSeq, humanGene, currentRNAacc, currentProtID, curentValues, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif):

    infoFileName = outputInfoFolder+os.path.sep+humanGene+"_excludedRNAtoCDS.info"

    # RNAaccToGeneFeatureInformation[RNAacc]=[geneID, chromosome, strand, chromStart, chromStop, CCDS]
    geneID = RNAaccToGeneFeatureInformation[currentRNAacc][0]
    chromosome = RNAaccToGeneFeatureInformation[currentRNAacc][1]
    strand = RNAaccToGeneFeatureInformation[currentRNAacc][2]
    chromStart = RNAaccToGeneFeatureInformation[currentRNAacc][3]
    chromStop = RNAaccToGeneFeatureInformation[currentRNAacc][4]
    if len(RNAaccToGeneFeatureInformation[currentRNAacc]) == 6 :
        CCDS = RNAaccToGeneFeatureInformation[currentRNAacc][5]

    with open (infoFileName,'w') as infoFile :
        #newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+seqName[1:]
        if len(RNAaccToGeneFeatureInformation[currentRNAacc]) == 6 :
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" [RNAacc="+currentRNAacc+"]"+" [geneID="+geneID+"]"+" [chromosome="+chromosome+"]"+" [strand="+strand+"]"+" [start="+chromStart+"]"+" [stop="+chromStop+"]"+" [CCDS="+CCDS+"]"+" [taxonClassif="+taxClassif+"]"
        else :
            newSeqName = ">[taxonID="+taxonID+"] "+"[taxonName="+taxonName+"] "+"[CDSLength="+str(len(seq))+"] "+"[CDSgbff="+curentValues[1]+"]"+" [RNAacc="+currentRNAacc+"]"+" [geneID="+geneID+"]"+" [chromosome="+chromosome+"]"+" [strand="+strand+"]"+" [start="+chromStart+"]"+" [stop="+chromStop+"]"+" [CCDS="+"NA"+"]"+" [taxonClassif="+taxClassif+"]"
        
        #infoFile.write(newSeqName+"\n"+newSeq+"\n")
        infoFile.write(newSeqName+"\n")


#---------------------------------------------------------------------------------          

########################   MAIN function   ########################
'''
  CORE OF THE PROGRAM :
  Function that generates the set of sequences in fasta format (and info format)
    - matches the list of human geneIDs with the taxon geneIDs
    - retrieves the longest CDS (in CDS GCF file)
    - recovers proteins according to the information of these CDS
    - recovers the RNA according to the information of these Proteins
    - makes new CDS (RNAtoCDS) from the informations obtained in the feature table and the RNA.gbff files.
'''

def addGenome(humanGeneIDList, orthologFile, GCFcdsFile, GCFproteinFile, GCFRNAFile, GCFFeatureTable, GCFRNA_GBFFFile, taxonID, taxonName, taxClassif, outputFastaFolder, outputInfoFolder):
#---------------------------------------------------------------------------------          
    # Creation of a dictionary of human orthologs gene name vs species considered
    # =(Recovery of the corresponding human geneID for the orthologous geneID of the taxon)
    taxonToHuman = {}
    speRef = "9606"
    refID = getRefID(humanGeneIDList)

    with open (orthologFile, 'r') as inputFile :   # gene_orthologs (NCBI)
        for line in inputFile :
            ortho_info = line.strip().split("\t")
            spe1= ortho_info[0]
            spe2= ortho_info[3]
            idRef = ortho_info[1]
            idTaxon = ortho_info[4]
            if (taxonID ==speRef) and (idRef in refID):  
                taxonToHuman[idRef] = idRef
            else :
                if(spe1==speRef) and (spe2==taxonID) and (idRef in refID):
                    taxonToHuman[idTaxon] = idRef
    
#---------------------------------------------------------------------------------          
    # Parse CDS GCF file to recover the longest CDS for each gene
    geneToLongestCDS = getLongestCDS(GCFcdsFile)
        
    # CDS sequences files contruction : CDS GCF file parsing
    with open (GCFcdsFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addSeq(seqName,seq,taxonToHuman,outputFastaFolder, outputInfoFolder,geneToLongestCDS, taxonName, taxonID, taxClassif)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addSeq(seqName,seq,taxonToHuman,outputFastaFolder, outputInfoFolder, geneToLongestCDS, taxonName, taxonID, taxClassif)
    
    # ProtID to human geneID correspondence dictionary :
    protIDtogeneID = {}
    for geneID, geneInfo in geneToLongestCDS.items() :
        if geneID in taxonToHuman.keys():
            humanGeneID = taxonToHuman[geneID]
            protID = geneInfo[2]
            protIDtogeneID[protID] = humanGeneID
            
#---------------------------------------------------------------------------------          
    # Proteins sequences files construction : Protein GCF file parsing
    with open (GCFproteinFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addProtSeq(seqName,seq,outputFastaFolder,protIDtogeneID, taxonName, taxonID, taxClassif)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addProtSeq(seqName,seq,outputFastaFolder, protIDtogeneID, taxonName, taxonID, taxClassif)

#---------------------------------------------------------------------------------          
    # RNA sequences files construction :
    # NB : Recovery of the right RNA thanks to the protID of the initial CDS compared to the RNA.gbff ans the feature_table files
    
    # RNA Feature_table parsing and build informations dictionary : (per RNA accession number)
    with open (GCFFeatureTable,'r') as featureTableFile : 
        RNAaccToGeneFeatureInformation = {}
        for line in featureTableFile :
            parse = line.strip().split("\t")
            if parse[0] == "mRNA": # recovery of the differents fields of interest (mRNA lines only)
                RNAacc = parse[10]
                chromosome = parse[5]
                genomicAcc = parse[6]
                chromStart = parse[7]
                chromStop = parse[8]
                chromStrand = parse[9]
                geneID = parse[15]
                RNAaccToGeneFeatureInformation[RNAacc]=[geneID, chromosome, chromStrand, chromStart, chromStop]
    
    # RNA GCF.gbff file parsing and build an other information dictionary
    with open (GCFRNA_GBFFFile,'r') as gbffFile :
        protIDtoRNAacc = {}
        # Dictionary : protIDtoRNAacc[protID]=['RNAacc','startCDS..endCDS']
        
        for line in gbffFile:
            if line.startswith("VERSION"):
                parse = line.strip().split("     ")         # 5 spaces !
                RNAacc = parse[1]
            if line.startswith("     "):                    # 5 spaces !
                parse = line.strip().split("     ")
                if parse[0] == "CDS" :
                    CDS = parse[2].strip()
            if line.startswith("                     "):    # 21 spaces !
                parse = line.strip().split("                     ")
                parse = parse[0]
                if parse.startswith("/protein_id="):
                    subline = parse.split("\"")
                    protID = subline[1]
                    if protID not in protIDtoRNAacc.keys():
                        protIDtoRNAacc[protID] = [RNAacc,CDS]
                if parse.startswith("/db_xref=\"CCDS:"):
                    subline = parse.split(":")
                    CCDS = subline[1][:-1]
                    RNAaccToGeneFeatureInformation[RNAacc].append(CCDS)

                # add: if startswith (/transl_except=)   when mutation out of start codon
                # exemple /transl_except=(pos:1277..1279,aa:Sec)

                # add: if startswith (/note=..)   when start codon mutation
                # exemple /note="non-AUG (CUG) translation initiation codon; TFIIIA;factor A"

    # So we have : RNAaccToGeneFeatureInformation[RNAacc]=[geneID, chromosome, chromStrand, chromStart, chromStop, CCDS] 


        # RNA and RNAtoCDS construction (+ infofile)
    with open (GCFRNAFile,'r') as fastaFile :
        seq = ""
        seqName = ""
        first = True
        for line in fastaFile :
            if line.startswith(">") :
                if first == False :
                    addRNASeq(seqName, seq, protIDtogeneID, protIDtoRNAacc, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)
                seqName = line.strip()
                seq = ""
                first = False
            else :
                seq += line.strip()
        if first == False :
            addRNASeq(seqName, seq, protIDtogeneID, protIDtoRNAacc, RNAaccToGeneFeatureInformation, outputFastaFolder, outputInfoFolder, taxonName, taxonID, taxClassif)



#---------------------------------------------------------------------------------        

if __name__ == "__main__" :
    if len(sys.argv) != 12 : # Requires twelve parameters
        usage()
        sys.exit()
    else :
        humanGeneIDList = sys.argv[1]       # 'humanGeneID.list' (path and name file)
        orthologFile = sys.argv[2]          # gene_orthologs file (NCBI)
        GCFcdsFile = sys.argv[3]            # GCF file for CDS sequences (NCBI)
        GCFproteinFile = sys.argv[4]        # GCF file for proteins sequences (NCBI)
        GCFRNAFile = sys.argv[5]            # GCF file for RNA sequences (NCBI)
        GCFFeatureTable = sys.argv[6]       # GCF file for gene informations (NCBI)
        GCFRNA_GBFFFile = sys.argv[7]       # GCF file for RNA transformation in CDS sequences (NCBI)
        taxonID = sys.argv[8]               # ex 9606 Homo Sapiens
        taxonName = sys.argv[9]             # ex Homo_Sapiens
        taxonClassif = sys.argv[10]         # ex Euarchontoglires
        outputFolder = sys.argv[10]         # work directory
        
        outputFastaFolder = outputFolder+os.path.sep+"FASTA"
        outputInfoFolder = outputFolder+os.path.sep+"INFO"
        os.system("mkdir "+outputFastaFolder)
        os.system("mkdir "+outputInfoFolder)
        
        addGenome(humanGeneIDList, orthologFile, GCFcdsFile, GCFproteinFile, GCFRNAFile, GCFFeatureTable, GCFRNA_GBFFFile, taxonID, taxonName, taxClassif, outputFastaFolder, outputInfoFolder)
