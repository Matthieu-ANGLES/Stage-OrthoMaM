#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import sys, os, re

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''
    
    This module checks the similarity between a RNAtoCDS reference and 
    its associated protein reference in OrthoMaM pipeline and make checks commands bash
    and others verifications (stats)
    
    Usage :  ./checkSequences.py outputFolder (contain FASTA folder)
    
    ''')

#------------------------------------------------------------------------------
'''
  Make list per gene for a 'type' of fasta files in Fasta Files Directory (NT, AA, RNA..)
'''

def makeListOfFastaFile (FastaFilesDirectory,typeFile) :
  
  listOfFastaFile = []
  if typeFile == "_excludedRNAtoCDS":
    arg = typeFile+".info"
  else :
    arg = typeFile+".fasta"

  listFiles = os.listdir(FastaFilesDirectory)

  for fileName in listFiles :
    ext = fileName.split("_")
    if ext[-1] == arg :
      listOfFastaFile.append(fileName)

  return (listOfFastaFile)

#------------------------------------------------------------------------------
'''
  MACSE translate CDS fasta files
'''

# CHEMINS A REVOIR SELON UTILISATION !!!! (chemin MACSE et répertoires scripts et d'enregistrementdes fastas)

def translateCDSMACSE(outputFolder):
  path = os.getcwd()
  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,"NT")
  cpt = 0
  for fastaFile in listFastaFiles :
    cpt+=1
    ext = fastaFile.split("_")
    geneID = ext[0]
    pathFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+fastaFile
    pathOutputFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+geneID+"_NTtranslate.fasta"
    cmdMACSE = "java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA -out_AA " + pathOutputFile + " -seq " + pathFile
    print ("Translate ",str(cpt),"/ ",str(len(listFastaFiles))," :")
    os.system(cmdMACSE)


#------------------------------------------------------------------------------
'''
  MACSE translate RNAtoCDS fasta files
'''

# CHEMINS A REVOIR SELON UTILISATION !!!! (chemin MACSE et répertoires scripts et d'enregistrementdes fastas)

def translateRNAtoCDSMACSE(outputFolder):
  path = os.getcwd()
  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,"RNAtoCDS")
  cpt = 0
  for fastaFile in listFastaFiles :
    cpt+=1
    ext = fastaFile.split("_")
    geneID = ext[0]
    pathFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+fastaFile
    pathOutputFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+geneID+"_RNAtoCDStranslated.fasta"
    cmdMACSE = "java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA -out_AA " + pathOutputFile + " -seq " + pathFile
    print (cmdMACSE)
    print ("Translate ",str(cpt),"/ ",str(len(listFastaFiles))," :")
    os.system(cmdMACSE)


#------------------------------------------------------------------------------
  #######################################################################################
  ### DEBUT COMMANDES BASH ###
'''
  Bash Commands Function
'''

def bashCommands (outputFolder):
  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"+os.path.sep

  print("\n////////////////////////////////////////////////////////////////////////\n")
  
  ### Contrôle sur les fichiers RNAtoCDS
  
  
  print ("Commandes Bash (analyses primaires) :\n")

  ### cmd0 et cmd1 ###
  cmd0 = "ls "+fastaFilesDirectory+"*_RNAtoCDStranslated.fasta | wc -l"
  print (cmd0)
  os.system(cmd0)
  for tax in "Canis_lupus_familiaris", "Mus_musculus", "Homo_sapiens" :
    cmd1 = "cat "+fastaFilesDirectory+"*_RNAtoCDStranslated.fasta | grep -c "+tax
    print (cmd1)
    os.system(cmd1)

  print("------------------------------------------------------------------------")
  ### cmd2 et cmd3 ###
  cptSeqWithoutCoreTaxon = 0
  for tax in "Canis_lupus_familiaris", "Mus_musculus", "Homo_sapiens" :
    cmd2 = "grep -c "+tax+" "+fastaFilesDirectory+"*_RNAtoCDStranslated.fasta | grep -c :0 "
    print (cmd2)
    os.system(cmd2)
    cmd3 = "grep -c "+tax+" "+fastaFilesDirectory+"*_RNAtoCDStranslated.fasta| grep :0 > "+"fastaFileWithout_"+tax+".list"
    os.system(cmd3)
    rescmd2 = os.popen(cmd2).readlines()
    cptSeqWithoutCoreTaxon += int(rescmd2[0])
  print ("Total : ", str(cptSeqWithoutCoreTaxon))
  print ("(Listes accessibles)")

  print("------------------------------------------------------------------------")
  ### Récupération informations de cmd3 ###
  fastaFileWithoutCanis = getGeneIDListInBashCommandFile("fastaFileWithout_Canis_lupus_familiaris.list")
  fastaFileWithoutMus = getGeneIDListInBashCommandFile("fastaFileWithout_Mus_musculus.list")
  fastaFileWithoutHomo = getGeneIDListInBashCommandFile("fastaFileWithout_Homo_sapiens.list")
  fastaFileWithoutCoreTaxons = fastaFileWithoutCanis + fastaFileWithoutMus + fastaFileWithoutHomo
  
  '''
  print ("CDS Files (GeneID) without Canis ("+str(len(fastaFileWithout_Canis_lupus_familiaris))+" ) : \n",fastaFileWithoutCanis,"\n")
  print ("CDS Files (GeneID) without Mus ("+str(len(fastaFileWithout_Mus_musculus))+" ) : \n",fastaFileWithoutMus,"\n")
  print ("CDS Files (GeneID) without Homo ("+str(len(fastaFileWithout_Homo_sapiens))+" ) : \n",fastaFileWithoutHomo,"\n")
  print ("Total : ",str(len(fastaFileWithoutCoreTaxons)))
  '''

  os.remove("fastaFileWithout_Canis_lupus_familiaris.list")
  os.remove("fastaFileWithout_Mus_musculus.list")
  os.remove("fastaFileWithout_Homo_sapiens.list")

  #print("------------------------------------------------------------------------")
 
  ### cmd4 et cmd5 ###
  cptSeqWithoutCoreTaxon = 0
  for tax in "Canis_lupus_familiaris", "Mus_musculus", "Homo_sapiens" :
    cmd4 = "grep -c "+tax+" "+fastaFilesDirectory+"*_AA.fasta | grep -c :0 "
    print (cmd4)
    os.system(cmd4)
    cmd5 = "grep -c "+tax+" "+fastaFilesDirectory+"*_AA.fasta| grep :0 > "+"fastaFileWithout_"+tax+".list"
    os.system(cmd5)
    rescmd4 = os.popen(cmd4).readlines()
    cptSeqWithoutCoreTaxon += int(rescmd4[0])
  print ("Total : ", str(cptSeqWithoutCoreTaxon))

  print("------------------------------------------------------------------------")
  ### Récupération informations de cmd5 ###
  fastaFileWithoutCanis = getGeneIDListInBashCommandFile("fastaFileWithout_Canis_lupus_familiaris.list")
  fastaFileWithoutMus = getGeneIDListInBashCommandFile("fastaFileWithout_Mus_musculus.list")
  fastaFileWithoutHomo = getGeneIDListInBashCommandFile("fastaFileWithout_Homo_sapiens.list")
  fastaFileWithoutCoreTaxons = fastaFileWithoutCanis + fastaFileWithoutMus + fastaFileWithoutHomo
  
  '''
  print ("CDS Files (GeneID) without Canis ("+str(len(fastaFileWithout_Canis_lupus_familiaris))+" ) : \n",fastaFileWithoutCanis,"\n")
  print ("CDS Files (GeneID) without Mus ("+str(len(fastaFileWithout_Mus_musculus))+" ) : \n",fastaFileWithoutMus,"\n")
  print ("CDS Files (GeneID) without Homo ("+str(len(fastaFileWithout_Homo_sapiens))+" ) : \n",fastaFileWithoutHomo,"\n")
  print ("Total : ",str(len(fastaFileWithoutCoreTaxons)))
  '''

  os.remove("fastaFileWithout_Canis_lupus_familiaris.list")
  os.remove("fastaFileWithout_Mus_musculus.list")
  os.remove("fastaFileWithout_Homo_sapiens.list")

 
  ### cmd6 ###
  infoFilesDirectory = outputFolder+os.path.sep+"INFO"+os.path.sep

  print ("cmd6")
  for tax in "Canis_lupus_familiaris", "Mus_musculus", "Homo_sapiens" :
    #cmd6 = "grep -c "+tax+" "+infoFilesDirectory+"*_excludedRNAtoCDS.info| grep :1 > "+"infoFileWith_"+tax+".list"
    cmd6 = "grep -c "+tax+" "+infoFilesDirectory+"*_excludedRNAtoCDS.info| grep :1 | wc -l"
    print (cmd6)
    os.system(cmd6)

  '''
  ### Récupération informations de cmd6 ###
  infoExcludedSeqFileWithCanis = getGeneIDListInBashCommandFile("infoFileWith_Canis_lupus_familiaris.list")
  infoExcludedSeqFileWithMus = getGeneIDListInBashCommandFile("infoFileWith_Mus_musculus.list")
  infoExcludedSeqFileWithHomo = getGeneIDListInBashCommandFile("infoFileWith_Homo_sapiens.list")
  excludedFileWithCoreTaxons = infoExcludedSeqFileWithCanis + infoExcludedSeqFileWithMus + infoExcludedSeqFileWithHomo

  print ("CDS Files (GeneID) with excluded sequence for Canis ("+str(len(infoExcludedSeqFileWithCanis))+" ) : \n",infoExcludedSeqFileWithCanis,"\n")
  print ("CDS Files (GeneID) with excluded sequence for Mus ("+str(len(infoExcludedSeqFileWithMus))+" ) : \n",infoExcludedSeqFileWithMus,"\n")
  print ("CDS Files (GeneID) with excluded sequence for Homo ("+str(len(infoExcludedSeqFileWithHomo))+" ) : \n",infoExcludedSeqFileWithHomo,"\n")
  print ("Total : ",str(len(excludedFileWithCoreTaxons)))
  
  #os.remove("infoFileWith_Canis_lupus_familiaris.list")
  #os.remove("infoFileWith_Mus_musculus.list")
  #os.remove("infoFileWith_Homo_sapiens.list")
  '''

  print("------------------------------------------------------------------------")
  ### Copie des fichiers FASTA '*_RNAtoCDStranslated.fasta' pour correction des condons stop ###
  ### Fichiers supplémentaires corrigés : '*_RNAtoCDStranslatedCorr.fasta' ###
  
  print ("Copie et correction des fichiers '*_RNAtoCDStranslated.fasta' :")

  cmd7 = "cd "+fastaFilesDirectory+";"+"find . -name '*_RNAtoCDStranslated.fasta' -print | xargs -i basename -s '.fasta' {} | xargs -i cp {}.fasta {}Corr.fasta ; sed -i \"s/*/X/g\" *Corr.fasta"
  os.system(cmd7)

  print ("The files '*_RNAtoCDStranslatedCorr.fasta' are created.")


  ### FIN COMMANDES BASH ###
  #######################################################################################


#------------------------------------------------------------------------------
'''
  REGEX function to access the geneID in the result of a bash command.
  Used in checkDownloadedSequences() function.
'''
def getGeneIDListInBashCommandFile(bashCommandFile):

  fastaFileWithoutTaxon = []
  with open(bashCommandFile, 'r') as file:
    for line in file:
      pattern = "(.+)(FASTA"+os.path.sep+")(\d+)([\_])"
      research = re.search(pattern,line)
      if research == None :
        continue
      geneID = research.group(3)
      fastaFileWithoutTaxon.append(geneID)

  return (fastaFileWithoutTaxon)

#------------------------------------------------------------------------------
'''
  REGEX function to access the taxonID of a header line of the fasta file.
  Used in differenceMaker() and toStringOutput() functions
'''
def getTaxonID(seqHeader):
  pattern = re.search("(.+)(taxonID=)([^\,\]]+)([\,\]])",seqHeader)    
  if pattern == None :
    return (None)
  taxonID = pattern.group(3)

  return (taxonID)
#------------------------------------------------------------------------------
'''
  Analysis function and creation of a header list and a sequence list from a given taxon list in a given fasta file.
  Used in toStringOutputFastaSequences() function.
'''
def headerAndSeqListsMaker (taxonIDList, fastaFilePath) :

  listSeq = []
  listHeader = []
  with open(fastaFilePath,'r') as inputFile :
    seeHeader = False
    for line in inputFile :
      if line.startswith(">"):
        currentTaxonID = getTaxonID(line)
        if currentTaxonID in taxonIDList :
          listHeader.append(line)
          seeHeader = True
      else :
        if seeHeader == True :
          listSeq.append(line)
          seeHeader = False

  return (listSeq, listHeader)

#------------------------------------------------------------------------------
'''
  Copy writing function.
  From a geneID and a taxonID passed as parameters, the function parses the fasta files 
  of the CDS sequences (RNAtoCDS) and the protein sequences and rewrites them into a new 
  output file whose name is also passed as parameter.
  Used checkDownloadedSequences() function.
'''

def toStringOutputFastaSequences(geneID, taxonIDList, outputFolder, outputNameFile):
  
  fastaFileDirectory = outputFolder+os.path.sep+"FASTA"+os.path.sep
  RNAtoCDSFastaFile = fastaFileDirectory+geneID+"_RNAtoCDStranslatedCorr.fasta"
  proteinFastaFile = fastaFileDirectory+geneID+"_AA.fasta"
  outputFile = outputFolder+os.path.sep+outputNameFile
  
  listSeqRNAtoCDS, listHeaderRNAtoCDS = headerAndSeqListsMaker (taxonIDList, RNAtoCDSFastaFile)
  listSeqProt, listHeaderProt = headerAndSeqListsMaker (taxonIDList, proteinFastaFile)

  with open (outputFile, 'a+') as output :
    output.write("# GeneID : "+geneID+"\n")

    for i in range (0,len(taxonIDList)-1) :
      output.write("-------------------------------------------------------------------\n")
      output.write("/TaxonID : "+taxonIDList[i]+"\n")
      output.write("RNAtoCDS sequence :\n")
      output.write(listHeaderRNAtoCDS[i])
      output.write(listSeqRNAtoCDS[i])
      output.write("Protein sequence :\n")
      output.write(listHeaderProt[i])
      output.write(listSeqProt[i])
    output.write("###################################################################\n")


#------------------------------------------------------------------------------
'''
  Make dictionnary for each sequence fastafile typeFile (k = geneFile , v = sequences)
'''

def parseFastaFile (outputFolder, typeFile):

  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,typeFile)

  seqDictionary = {}
  headerDictionary = {}
  for fastaFile in listFastaFiles :
    ext = fastaFile.split("_")
    geneID = ext[0]
    pathFile = fastaFilesDirectory+os.path.sep+fastaFile
    with open (pathFile,'r') as file :
      head = ""
      seq = ""
      first = True
      for line in file :
        if line.startswith(">") :
          head = line.strip()
          if geneID not in headerDictionary.keys():
            headerDictionary[geneID]=head+"$"
          else :
            headerDictionary[geneID]+=head+"$"
          if first == False :
            if geneID not in seqDictionary.keys():
              seqDictionary[geneID]=seq+"$"
            else :
              seqDictionary[geneID]+=seq+"$"
          seq = ""
          first = False
        else :
          seq += line.strip()
      if first == False :
        if geneID not in seqDictionary.keys():
          seqDictionary[geneID]=seq+"$"
        else :
          seqDictionary[geneID]+=seq+"$"

  return (seqDictionary, headerDictionary)

#------------------------------------------------------------------------------
'''
  Make a copy of fasta file
'''
def copyFastaFile (outputFastaFolder, type, targetGeneID, copyFolder):

  listFastaFiles = makeListOfFastaFile(outputFastaFolder,type)
  for fastaFile in listFastaFiles :
    ext = fastaFile.split("_")
    geneID = ext[0]
    targetFastaFile = outputFastaFolder+os.path.sep+geneID+"_"+type+".fasta"
    copyTargetFastaFile = copyFolder+os.path.sep+geneID+"_"+type+".fasta"
    if geneID == targetGeneID :
      os.system("cp "+targetFastaFile+" "+copyTargetFastaFile)

#------------------------------------------------------------------------------
'''
  Find the difference between two sequences and give position
'''
def differenceMaker(RNAtoCDSSequence, RNAtoCDSHeader, proteinSequence, proteinHeader, orthomamTaxonIDList):

  splitSeq = RNAtoCDSSequence.split("$")
  listRNAtoCDSSequence = splitSeq[:-1]
  splitHead = RNAtoCDSHeader.split("$")
  listRNAtoCDSHeader = splitHead[:-1]

  splitSeq = proteinSequence.split("$")
  listProteinSequence = splitSeq[:-1]
  splitHead = proteinHeader.split("$")
  listProteinHeader = splitHead[:-1]
  
  nbrSeqRNAtoCDS = len(listRNAtoCDSSequence)
  nbrSeqProt = len(listProteinSequence)
  nbrTaxonTotal = len(orthomamTaxonIDList)
 
  # Dictionnaries Informations :
  taxonPresenceInProtFile = {}
  #   taxonPresenceInProtFile[taxonID] = ['bool-presence']
  #   lenDict = 1   if 'presence'  = [True]
  #   lenDict = 1   if 'absence'   = [False]
  
  taxonPresenceInRNAtoCDSFile = {}
  #   taxonPresenceInRNAtoCDSFile[taxonID] = ['bool-presence','bool-stopCodonCorrectionNecessity', 'bool-equality']
  #   lenDict = 1   if 'absence'  = [False]
  #   lenDict = 3   if 'presence','No-stopCodonCorrection','equality'  = [True,False,True]
  #   lenDict = 3   if 'presence','stopCodonCorrection','equality'     = [True,True,True]
  #   lenDict = 4   if 'presence','No-stopCodonCorrection','No-equality','position diff' = [True,False,False,list-position]
  #   lenDict = 4   if 'presence','stopCodonCorrection','No-equality','position diff'    = [True,True,False,list-position]


  # AJOUTER RECUPERATION DES ACIDES AMINES CONCERNES POUR LES MUTATIONS

  
  cptRNAtoCDS = 0
  cptProt = 0
  cptTaxon = 0
  
  while cptRNAtoCDS < nbrSeqRNAtoCDS and cptProt < nbrSeqProt and cptTaxon < nbrTaxonTotal :
    #print("Entrée dans la boucle ", cptRNAtoCDS, cptProt, cptTaxon)
    currentRNAtoCDSTaxonID = getTaxonID(listRNAtoCDSHeader[cptRNAtoCDS])
    currentProtTaxonID = getTaxonID(listProteinHeader[cptProt])
    taxonID = orthomamTaxonIDList[cptTaxon]

    # Taxon in files
    if taxonID == currentRNAtoCDSTaxonID and taxonID == currentProtTaxonID :
      taxonPresenceInProtFile[taxonID] = [True]
      # Equality
      if listRNAtoCDSSequence[cptRNAtoCDS] == listProteinSequence[cptProt] :
        taxonPresenceInRNAtoCDSFile[taxonID] = [True,False,True]
        cptRNAtoCDS+=1
        cptProt +=1
        cptTaxon +=1
      # stopCodonCorrection
      else :
        if listRNAtoCDSSequence[cptRNAtoCDS].count('*') != 0 :
          listRNAtoCDSSequence[cptRNAtoCDS] = listRNAtoCDSSequence[cptRNAtoCDS].replace('*','X')
          position = getPositionsDifferences(listRNAtoCDSSequence[cptRNAtoCDS],listProteinSequence[cptProt])
          # Equality with stopCodonCorrection
          if listRNAtoCDSSequence[cptRNAtoCDS] == listProteinSequence[cptProt] :
            taxonPresenceInRNAtoCDSFile[taxonID] = [True,True,True]
            cptRNAtoCDS +=1
            cptProt +=1
            cptTaxon +=1
          else :
          # Inequality with stopCodonCorrection
            taxonPresenceInRNAtoCDSFile[taxonID] = [True,True,False,position]
            cptRNAtoCDS +=1
            cptProt +=1
            cptTaxon +=1
        else :
          # Inegality without stopCodonCorrection necessity
          position = getPositionsDifferences(listRNAtoCDSSequence[cptRNAtoCDS],listProteinSequence[cptProt])        
          taxonPresenceInRNAtoCDSFile[taxonID] = [True,False,False,position]
          cptRNAtoCDS +=1
          cptProt +=1
          cptTaxon +=1
        
    else :
      # Taxon in Prot File only
      if taxonID != currentRNAtoCDSTaxonID and taxonID == currentProtTaxonID :
        taxonPresenceInRNAtoCDSFile[taxonID] = [False]
        taxonPresenceInProtFile[taxonID] = [True]
        cptProt +=1
        cptTaxon +=1
      # Taxon not in RNAtoCDS and Prot file
      elif taxonID != currentRNAtoCDSTaxonID and taxonID != currentProtTaxonID :
        taxonPresenceInRNAtoCDSFile[taxonID] = [False]
        taxonPresenceInProtFile[taxonID] = [False]
        cptTaxon +=1

  for taxonID in orthomamTaxonIDList :
    if taxonID not in taxonPresenceInRNAtoCDSFile :
      taxonPresenceInRNAtoCDSFile[taxonID] = [False]
    if taxonID not in taxonPresenceInProtFile :
      taxonPresenceInProtFile[taxonID] = [False]
  
  return (taxonPresenceInRNAtoCDSFile, taxonPresenceInProtFile)
  
#------------------------------------------------------------------------------
'''
  Use in differenceMaker function  
  Give a vector of positions of the difference between CDS sequence vs Prot sequence
'''
def getPositionsDifferences(CDS, Prot):
  i = 0
  positions = []
  while i<len(CDS) and i<len(Prot) :
    if CDS[i] == Prot[i] :
      i+=1
    else :
      positions.append(i+1)
      i+=1

  return (positions) # ajouter aa CDS et aa prot

#------------------------------------------------------------------------------
'''

'''

def getTaxonStats (taxonID, dicoTranslatedRNAtoCDS_WithDollar, geneIDtoRNAtoCDSFilesStats, geneIDtoProtFilesStats):
  presentSequencesInRNAtoCDSFile = []
  presentSequencesInProtFile = []
  goodSequencesWhenCodonStopCorrection = []
  goodSequencesConsideredWhenStartCodonMutationOnly = []
  badSequencesWithOtherMutations = []

  for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
    extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
    extractPresenceProt = geneIDtoProtFilesStats[geneID][taxonID][0]
    if extractPresenceProt == True :
      presentSequencesInProtFile.append(geneID)
    if extractPresenceRNAtoCDS == True :
      presentSequencesInRNAtoCDSFile.append(geneID)
      extractStopCodonCorrectionRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][1]
      extractEqualityRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][2]
      if extractStopCodonCorrectionRNAtoCDS == True and extractEqualityRNAtoCDS == True :
        goodSequencesWhenCodonStopCorrection.append(geneID)

      if extractEqualityRNAtoCDS == False :
        extractPositionMutation = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][3]
        if len(extractPositionMutation) == 1:
          if extractPositionMutation[0] == 1 :
            goodSequencesConsideredWhenStartCodonMutationOnly.append(geneID)
        if len(extractPositionMutation) != 0 and geneID not in goodSequencesConsideredWhenStartCodonMutationOnly :
          badSequencesWithOtherMutations.append(geneID)

  return (presentSequencesInRNAtoCDSFile, presentSequencesInProtFile, goodSequencesWhenCodonStopCorrection, goodSequencesConsideredWhenStartCodonMutationOnly, badSequencesWithOtherMutations)

#------------------------------------------------------------------------------
'''
  Compare RNAtoCDS sequence and AA sequences per gene File.
    Used methods : translateMACSE, makeListOfFastaFile, parseFastaFile 
'''

def checkDownloadedSequences(outputFolder):

  # Translate CDS fasta file (*_NT.fasta) with MACSE (create *_NTtranslate.fasta files) #####
  #translateCDSMACSE(outputFolder)

  # Translate RNAtoCDS fasta file (*_RNAtoCDS.fasta) with MACSE (create *_RNAtoCDStranslated.fasta files) #####
  #translateRNAtoCDSMACSE(outputFolder)

  # Bash primary stats (Homo, Mus, Canis)
  #bashCommands (outputFolder)

  # Open GCF.list for create three lists for taxons (taxonName taxonID and taxonClassif)
  GCFlistPath = outputFolder+os.path.sep+"GCF"+os.path.sep+"GCF.list"
  orthomamTaxonNameList = []
  orthomamTaxonIDList = []
  orthomamTaxonClassifList = []

  with open (GCFlistPath,'r') as GCFlist :
    for line in GCFlist :
      parse = line.strip().split("\t")
      orthomamTaxonNameList.append(parse[0])
      orthomamTaxonClassifList.append(parse[-1])
      orthomamTaxonIDList.append(parse[-2])
  
  print("\n////////////////////////////////////////////////////////////////////////\n")
  
  outputFile = outputFolder+os.path.sep+"resume_fastaSequences.stats"
  with open (outputFile, 'a+') as output :
    for i in range(0,len(orthomamTaxonNameList)):
      print (orthomamTaxonNameList[i]," / ",orthomamTaxonIDList[i]," / ",orthomamTaxonClassifList[i])
      output.write(str(orthomamTaxonNameList[i])+" / "+str(orthomamTaxonIDList[i])+" / "+str(orthomamTaxonClassifList[i])+"\n")
    output.write ("------------------------------------------------------------------------\n")

  # Dictionnary with '$' behind each sequences (AVEC DICO SEQ ET DICO HEADER)
  dicoCDS_WithDollar, dicoCDS_header_WithDollar = parseFastaFile(outputFolder, "NT")
  dicoProt_WithDollar, dicoProt_header_WithDollar = parseFastaFile(outputFolder, "AA")
  dicoRNA_WithDollar, dicoRNA_header_WithDollar = parseFastaFile(outputFolder, "RNA")
  dicoRNAtoCDS_WithDollar, dicoRNAtoCDS_header_WithDollar = parseFastaFile(outputFolder, "RNAtoCDS")
  dicoTranslatedRNAtoCDS_WithDollar, dicoTranslatedRNAtoCDS_header_WithDollar = parseFastaFile(outputFolder, "RNAtoCDStranslated")

  print("\n////////////////////////////////////////////////////////////////////////\n")


  # Stats : RNAtoCDStranslated.fasta sequence with differenceMaker() function (line 255)
  
  geneIDtoRNAtoCDSFilesStats = {}
  geneIDtoProtFilesStats = {}
  for gene in dicoTranslatedRNAtoCDS_WithDollar.keys():
    taxonPresenceInRNAtoCDSFile, taxonPresenceInProtFile = differenceMaker(dicoTranslatedRNAtoCDS_WithDollar[gene],dicoTranslatedRNAtoCDS_header_WithDollar[gene],dicoProt_WithDollar[gene],dicoProt_header_WithDollar[gene], orthomamTaxonIDList)
    geneIDtoRNAtoCDSFilesStats[gene] = taxonPresenceInRNAtoCDSFile
    geneIDtoProtFilesStats[gene] = taxonPresenceInProtFile
  
  # Dictionnaries Informations : in differenceMaker() function
  # ==> taxonPresenceInProtFile[taxonID] = ['bool-presence']
  #       lenDict = 1   if 'presence'  = [True]
  #       lenDict = 1   if 'absence'   = [False]
  
  # ==> taxonPresenceInRNAtoCDSFile[taxonID] = ['bool-presence','bool-stopCodonCorrectionNecessity', 'bool-equality']
  #       lenDict = 1   if 'absence'  = [False]
  #       lenDict = 3   if 'presence','No-stopCodonCorrection','equality'  = [True,False,True]
  #       lenDict = 3   if 'presence','stopCodonCorrection','equality'     = [True,True,True]
  #       lenDict = 4   if 'presence','No-stopCodonCorrection','No-equality','position diff' = [True,False,False,list-position]
  #       lenDict = 4   if 'presence','stopCodonCorrection','No-equality','position diff'    = [True,True,False,list-position]

  outputFile = outputFolder+os.path.sep+"resume_fastaSequences.stats"
  with open (outputFile, 'a+') as output :
    print("Parsing of the files and use of the dictionaries (verifications/stats) :")
    output.write("Number of geneID fasta files : "+ str(len(geneIDtoRNAtoCDSFilesStats))+"\n")
    print ("Number of geneID fasta files (RNAtoCDS) : ", len(geneIDtoRNAtoCDSFilesStats),"\n")
    output.write ("------------------------------------------------------------------------\n")

  # Fichier avec taxons pillier et avec tous les taxons :
  #print ("Fichier avec taxons pillier et avec tous les taxons :")
  geneIDWithAllTaxons = []
  geneIDWithCoreTaxons = []
  coreTaxonsList = ['9615','10090','9606']
  for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
    cptTaxons = []
    cptCoresTaxons = []
    for taxonID in orthomamTaxonIDList :
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        cptTaxons.append(taxonID)
        if taxonID in coreTaxonsList :
          cptCoresTaxons.append(geneID)

    if len(cptTaxons) == (len(orthomamTaxonIDList)):
      geneIDWithAllTaxons.append(geneID)
    if len(cptCoresTaxons) == 3 :
      geneIDWithCoreTaxons.append(gene)

  #print (len(geneIDWithCoreTaxons)," fichiers RNAtoCDS ont au moins les séquences pour Homo, Mus et Canis")
  #print (len(geneIDWithAllTaxons)," fichiers RNAtoCDS ont les séquences de tous les taxons")
  #print ("(sur ",len(geneIDtoRNAtoCDSFilesStats),")")
  #print("------------------------------------------------------------------------\n")


  # Nbr de séquences avec ou sans discordance avec la référence protéique associé :
  #print ("Nbr de séquences avec ou sans discordance avec la référence protéique associé :")
  geneIDWithAllSequenceOk = []
  geneIDWithSequenceNotOk = []
  geneIDWithMoreOfOnTaxonHaveSeqNotOk = []
  cptSeq = 0
  cptOk = 0
  cptNotOk = 0
  cptCodonStartMutation = 0
  dicoGeneIDtoTaxonIDWithSeqNotOk = {}
  dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned = {}
  for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
    for taxonID in orthomamTaxonIDList :
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        cptSeq +=1
        extractGeneIDWithEquality = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][2]
        if extractGeneIDWithEquality == True :
          cptOk +=1
        else :
          cptNotOk +=1
          if geneID not in dicoGeneIDtoTaxonIDWithSeqNotOk.keys():
            dicoGeneIDtoTaxonIDWithSeqNotOk[geneID] = [taxonID]
          else :
            dicoGeneIDtoTaxonIDWithSeqNotOk[geneID].append(taxonID)
          extractgeneIDWithMutationOnFirstAminoAcid = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][3]
          if len(extractgeneIDWithMutationOnFirstAminoAcid) == 1:
            if extractgeneIDWithMutationOnFirstAminoAcid[0] == 1 :
              cptCodonStartMutation +=1
              if geneID not in dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned.keys():
                dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned[geneID] = [taxonID]
              else :
                dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned[geneID].append(taxonID)
    #pour vérif
    if geneID in dicoGeneIDtoTaxonIDWithSeqNotOk.keys():
      if len(dicoGeneIDtoTaxonIDWithSeqNotOk[geneID]) > 1 :
        #print ("plusieurs taxons ont des erreurs pour ce gene")
        geneIDWithMoreOfOnTaxonHaveSeqNotOk.append(geneID)

  dicoGeneIDtoTaxonIDWithSeqNotOk_WithoutStartCodonMutation = {}
  for geneID, taxonList in dicoGeneIDtoTaxonIDWithSeqNotOk.items():
    if geneID not in dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned.keys():
      dicoGeneIDtoTaxonIDWithSeqNotOk_WithoutStartCodonMutation[geneID] = taxonList
  
  
  outputFile = outputFolder+os.path.sep+"resume_fastaSequences.stats"
  with open (outputFile, 'a+') as output :
    print("------------------------------------------------------------------------\n")
    output.write("Out of "+str(cptSeq)+" translated RNAtoCDS (total in all files) :\n")
    print ("Out of ",cptSeq," translated RNAtoCDS (total in all files) :")
    
    output.write("  - "+str(cptOk)+ " are similar to the corresponding protein.\n")
    print ("  - ",cptOk, " are similar to the corresponding protein.")
    
    output.write("  - "+str(cptNotOk)+ " are discordant. (of which, "+str(cptCodonStartMutation)+" concern the start codon only)\n")
    print ("  - ",cptNotOk, " are discordant. (of which, ",cptCodonStartMutation," concern the start codon only)\n")
    output.write("See 'sequencesWithMutations_ExceptOnlyFirstCodonMutation.seq'\n")
    output.write("See 'sequencesWithOnlyFirstCodonMutation.seq'\n")
    output.write("------------------------------------------------------------------------\n")



  # Output File
  cptSeq = 0
  for geneID in dicoGeneIDtoTaxonIDWithSeqNotOk_WithoutStartCodonMutation.keys():
    taxonIDList = dicoGeneIDtoTaxonIDWithSeqNotOk_WithoutStartCodonMutation[geneID]
    cptSeq += len(taxonIDList)
    outputNameFile = "sequencesWithMutations_ExceptOnlyFirstCodonMutation.seq"
    toStringOutputFastaSequences(geneID, taxonIDList, outputFolder, outputNameFile)
  print ("The file 'sequencesWithMutations_ExceptOnlyFirstCodonMutation.seq' was created.")
  print ("(",cptSeq," sequences in ",len(dicoGeneIDtoTaxonIDWithSeqNotOk_WithoutStartCodonMutation.keys())," fasta files)")

  cptSeq=0
  for geneID in dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned.keys():
    taxonIDList = dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned[geneID]
    cptSeq += len(taxonIDList)
    outputNameFile = "sequencesWithOnlyFirstCodonMutation.seq"
    toStringOutputFastaSequences(geneID, taxonIDList, outputFolder, outputNameFile)
  print ("The file 'sequencesWithOnlyFirstCodonMutation.seq' was created.")
  print ("(",cptSeq," sequences in ",len(dicoGeneIDtoTaxonIDWithSeqNotOk_startCodonConcerned.keys())," fasta files)")
  
  print("\n------------------------------------------------------------------------\n")


  ####### ANCIENNE VERSION DE STATS A PARTIR DE LA = CONFRONTATIONS/TESTS AVEC COMMANDES BASH #######
  '''
  # Taxons avec mutation unique en première position :
  #print("Taxons avec mutation unique en première position :")
  for taxonID in orthomamTaxonIDList :
    geneIDWithMutationOnFirstAminoAcid =[]
    for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        extractGeneIDWithEquality = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][2]
        if extractGeneIDWithEquality == False :
          extractgeneIDWithMutationOnFirstAminoAcid = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][3]
          if len(extractgeneIDWithMutationOnFirstAminoAcid) == 1:
            if extractgeneIDWithMutationOnFirstAminoAcid[0] == 1 :
              geneIDWithMutationOnFirstAminoAcid.append(geneID)
    
    #print ("Taxon ",taxonID," :  a une différence d'acide aminé en première position dans ",len(geneIDWithMutationOnFirstAminoAcid)," fichiers RNAtoCDS ")
  #print("------------------------------------------------------------------------\n")
  '''

  '''
  # Fichiers qui ont des séquences différentes après traduction Macse des RNAtoCDS (hors codonStart) :
  print ("Fichiers qui ont des séquences différentes après traduction Macse des RNAtoCDS (hors codonStart) :")
  geneIDWithAllSequenceOk = []
  geneIDWithSequenceNotOk = []
  dicoGeneIDStartCodonModif = {}
  for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
    cptSeq = 0
    cptTaxon = 0
    for taxonID in orthomamTaxonIDList :
      cptTaxon +=1
      #print ("geneID : ",geneID,"taxonID : ",taxonID, cptTaxon)
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        extractGeneIDWithEquality = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][2]
        if extractGeneIDWithEquality == True :
          cptSeq +=1
        if extractGeneIDWithEquality == False :
          extractgeneIDWithMutationOnFirstAminoAcid = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][3]
          if len(extractgeneIDWithMutationOnFirstAminoAcid) == 1:
            if extractgeneIDWithMutationOnFirstAminoAcid[0] == 1 :
              if geneID not in dicoGeneIDStartCodonModif.keys():
                dicoGeneIDStartCodonModif[geneID] = [taxonID]
              else:
                dicoGeneIDStartCodonModif[geneID].append(taxonID)
              cptSeq +=1
    if cptSeq == cptTaxon :
      geneIDWithAllSequenceOk.append(geneID)
    else :
      geneIDWithSequenceNotOk.append(geneID)
  
  cptExcptTrans = 0
  for k,v in dicoGeneIDStartCodonModif.items():
    extracTaxNb = len(v)
    cptExcptTrans += extracTaxNb

  print (len(geneIDWithAllSequenceOk)," fichiers RNAtoCDS ont toutes leurs séquences sans discordance avec les séquences protéiques")
  print ("(",cptExcptTrans,"'Exception translation' codon start comptées comme correctes)")
  print (len(geneIDWithSequenceNotOk)," fichiers RNAtoCDS ont des discordances (ailleurs)")
  print("------------------------------------------------------------------------\n")
  '''
  '''
  # Absence/Présence taxons dans les fichiers RNAtoCDS(_RNAtoCDStranslated.fasta) et Protéiques (_AA) :
  #print ("Absence taxons dans les fichiers RNAtoCDS(_RNAtoCDStranslated.fasta) et Protéiques (_AA) :")
  for taxonID in orthomamTaxonIDList :
    absenceRNAtoCDS = []
    absenceProt = []
    for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      extractPresenceProt = geneIDtoProtFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == False :
        absenceRNAtoCDS.append(geneID)
      if extractPresenceProt == False :
        absenceProt.append(geneID)
    #print ("Taxon ",taxonID," :  absent dans ",len(absenceRNAtoCDS)," fichiers RNAtoCDS et dans ",len(absenceProt)," fichiers Prot")
  #print("------------------------------------------------------------------------")
  '''

  '''
  # Présence taxons dans les fichiers RNAtoCDS(_RNAtoCDStranslated.fasta) et Protéiques (_AA) :
  for taxonID in orthomamTaxonIDList :
    presenceRNAtoCDS = []
    presenceProt = []
    for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      extractPresenceProt = geneIDtoProtFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        presenceRNAtoCDS.append(geneID)
      if extractPresenceProt == True :
        presenceProt.append(geneID)
    #print ("Taxon ",taxonID," :  présent dans ",len(presenceRNAtoCDS)," fichiers RNAtoCDS et dans ",len(presenceProt)," fichiers Prot")
  #print("------------------------------------------------------------------------\n")
  '''

  '''
  # Correction des codons stop :
  for taxonID in orthomamTaxonIDList :
    geneIDWithCorrection =[]
    geneIDWithoutCorrection =[]
    for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        extractStopCodonCorrection = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][1]
        if extractStopCodonCorrection == False :
          geneIDWithoutCorrection.append(geneID)
        if extractStopCodonCorrection == True :
          geneIDWithCorrection.append(geneID)

    #print ("Taxon ",taxonID," :  a necessité une correction 'codonStop' dans ",len(geneIDWithCorrection)," fichiers RNAtoCDS ")
    #print ("Taxon ",taxonID," :  n'a pas necessité de correction dans ",len(geneIDWithoutCorrection)," fichiers RNAtoCDS ")
    #print("------")
    #print ("Total : ",len(geneIDWithCorrection)+len(geneIDWithoutCorrection))
  #geneIDWithCorrection.sort()
  #print (geneIDWithCorrection)

  print("\n////////////////////////////////////////////////////////////////////////\n")
  '''


####### NOUVELLE VERSION A PARTIR DE LA #######

  '''
  On garde 15398  séquences.
  Pour Homo_sapiens (pour Mus etc... Il faut une section par espèces, et il faut tout faire tourner avec toutes les espèces)

    - X sont différentes quand on les traduit
    - Y sont égales quand on met X à la place des *
    - Z sont égales quand le codon start est considéré comme "translation exception" 
    - Q sont incomplètes (celles avec <, > dans le CDS)
    - P ont plus d'une mutation/pas que dans le codonStart
  '''
  # Utilisation de la fonction getTaxonStats () par taxonID :
  
  outputFile = outputFolder+os.path.sep+"resume_fastaSequences.stats"
  infoFilesDirectory = outputFolder+os.path.sep+"INFO"+os.path.sep

  with open (outputFile, 'a+') as output :
    i = 0
    for taxonID in orthomamTaxonIDList :
      presentSequencesInRNAtoCDSFile, presentSequencesInProtFile, goodSequencesWhenCodonStopCorrection, goodSequencesConsideredWhenStartCodonMutationOnly, badSequencesWithOtherMutations = getTaxonStats (taxonID, dicoTranslatedRNAtoCDS_WithDollar, geneIDtoRNAtoCDSFilesStats, geneIDtoProtFilesStats)
      
      cmdBash = "grep -c "+taxonID+" "+infoFilesDirectory+"*_excludedRNAtoCDS.info| grep :1 | wc -l"
      excludedSequences = os.popen(cmdBash).readlines()[0].strip()

      output.write("# Taxon : "+orthomamTaxonNameList[i]+"("+taxonID+")\n")
      print ("# Taxon : ",orthomamTaxonNameList[i],"(",taxonID,")")
      i+=1

      output.write("  We have : "+str(len(presentSequencesInRNAtoCDSFile))+"  RNAtoCDS sequences for "+str(len(presentSequencesInProtFile))+" protein sequences.\n")
      print ("  We have : ",len(presentSequencesInRNAtoCDSFile)," RNAtoCDS sequences for ",len(presentSequencesInProtFile)," protein sequences..")

      output.write("  "+str(len(goodSequencesWhenCodonStopCorrection))+" are concordant with the protein sequence when 'X' instead of '*'.\n")
      print("  ",str(len(goodSequencesWhenCodonStopCorrection))," are concordant with the protein sequence when 'X' instead of '*'.")

      output.write("  "+str(len(goodSequencesConsideredWhenStartCodonMutationOnly))+" have a single mutation in position 1 (start codon, 'non-AUG')."+"\n")
      print("  ",str(len(goodSequencesConsideredWhenStartCodonMutationOnly))," have a single mutation in position 1 (start codon, 'non-AUG').")

      output.write("  "+str(excludedSequences)+" are excluded because one of the CDS bounds is outside the RNA sequence (> or <)"+"\n")
      print("  ",str(excludedSequences)," are excluded because one of the CDS bounds is outside the RNA sequence (> or <)")

      output.write("  "+str(len(badSequencesWithOtherMutations))+" have more than one mutation (except single mutation in position 1).\n\n")
      print("  ",str(len(badSequencesWithOtherMutations))," have more than one mutation (except single mutation in position 1).\n")

      
  print ("The file 'resume_fastaSequences.stats' was created.\n")


  # TESTS

  '''
  print(geneIDtoRNAtoCDSFilesStats['147040']['9606'][3])
  
  print (geneIDtoRNAtoCDSFilesStats['147040']['9606'])
  print (geneIDtoProtFilesStats['147040']['9606'])
  
  print(dicoTranslatedRNAtoCDS_header_WithDollar['147040'])
  
  print(dicoTranslatedRNAtoCDS_WithDollar['147040'])

  print(dicoProt_WithDollar['147040'])
  
  # voir pourquoi quelques seq ne sont pas affichées dans le fichier sequencesWithOnlyFirstCodonMutation.seq
  '''


#------------------------------------------------------------------------------
'''
  Main
'''

if __name__ == "__main__" :
  if len (sys.argv) != 2:
      usage()
      sys.exit()
  else :
    outputFolder = sys.argv[1]  
    checkDownloadedSequences(outputFolder)

    # ATTENTION AUX LIGNES COMMENTEES SUIVANTES :
        # Debut fonction checkDownloadedSequences : Translate RNAtoCDS fasta file (*_RNAtoCDS.fasta) with MACSE
        # (à commenter si traductions déja faites)

        # Commande Bash (cmd7) : RNAtoCDStranslated.fasta >> RNAtoCDStranslatedCorr.fasta
        # (à commenter si les fichiers sont déja corrigés)