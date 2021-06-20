#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import sys, os, re
import Levenshtein as lev

#------------------------------------------------------------------------------
'''
Help message
'''

def usage() :
    
    sys.stderr.write('''
    
    This module checks the similarity between a CDS reference and 
    its associated protein reference in OrthoMaM pipeline and make checks commands bash
    
    Usage :  ./checkSequences.py outputFolder    (contain FASTA folder)
    
    ''')

#------------------------------------------------------------------------------
'''
  Levenshtein algo
  Distance: This is the measure of the minimum number of changes (insertions, deletions or substitutions) that must be made to change a sequence from one word to another.
  Ratio: Similarity between two words
'''

def levCalclulate(seq1, seq2):
  distance = lev.distance(seq1, seq2)
  ratio = lev.ratio(seq1, seq2)
  return (distance, ratio)


#------------------------------------------------------------------------------
'''
  Make list for NT files or AA files per gene of Fasta Files Directory
'''

def makeListOfFastaFile (FastaFilesDirectory,type) :
  
  listOfFastaFile = []
  arg = type+".fasta"
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
    pathOutputFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+geneID+"_RNAtoCDStranslate.fasta"
    cmdMACSE = "java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA -out_AA " + pathOutputFile + " -seq " + pathFile
    print ("Translate ",str(cpt),"/ ",str(len(listFastaFiles))," :")
    os.system(cmdMACSE)


#------------------------------------------------------------------------------
'''
  Bash Commands Function
'''

def bashCommands (outputFolder):
  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"+os.path.sep

  print("\n////////////////////////////////////////////////////////////////////////\n")
  
  '''
  ### Contrôle sur les fichiers RNAtoCDS
  '''
  
  print ("Commandes Bash (analyses primaires) :\n")

  ### cmd0 et cmd1 ###
  cmd0 = "ls "+fastaFilesDirectory+"*_RNAtoCDStranslate.fasta | wc -l"
  print (cmd0)
  os.system(cmd0)
  for tax in "Canis", "Mus", "Homo" :
    cmd1 = "cat "+fastaFilesDirectory+"*_RNAtoCDStranslate.fasta | grep -c "+tax
    print (cmd1)
    os.system(cmd1)

  print("------------------------------------------------------------------------")
  ### cmd2 et cmd3 ###
  cptSeqWithoutCoreTaxon = 0
  for tax in "Canis", "Mus", "Homo" :
    cmd2 = "grep -c "+tax+" "+fastaFilesDirectory+"*_RNAtoCDStranslate.fasta | grep -c :0 "
    print (cmd2)
    os.system(cmd2)
    cmd3 = "grep -c "+tax+" "+fastaFilesDirectory+"*_RNAtoCDStranslate.fasta| grep :0 > "+"fastaFileWithout_"+tax+".list"
    os.system(cmd3)
    rescmd2 = os.popen(cmd2).readlines()
    cptSeqWithoutCoreTaxon += int(rescmd2[0])
  print ("Total : ", str(cptSeqWithoutCoreTaxon))
  print ("(Listes accessibles)")

  print("------------------------------------------------------------------------")
  ### Récupération informations de cmd3 ###
  fastaFileWithoutCanis = getGeneIDListInBashCommandFile("fastaFileWithout_Canis.list")
  fastaFileWithoutMus = getGeneIDListInBashCommandFile("fastaFileWithout_Mus.list")
  fastaFileWithoutHomo = getGeneIDListInBashCommandFile("fastaFileWithout_Homo.list")
  fastaFileWithoutCoreTaxons = fastaFileWithoutCanis + fastaFileWithoutMus + fastaFileWithoutHomo
  
  '''
  print ("CDS Files (GeneID) without Canis ("+str(len(fastaFileWithoutCanis))+" ) : \n",fastaFileWithoutCanis,"\n")
  print ("CDS Files (GeneID) without Mus ("+str(len(fastaFileWithoutMus))+" ) : \n",fastaFileWithoutMus,"\n")
  print ("CDS Files (GeneID) without Homo ("+str(len(fastaFileWithoutHomo))+" ) : \n",fastaFileWithoutHomo,"\n")
  print ("Total : ",str(len(fastaFileWithoutCoreTaxons)))
  '''

  os.remove("fastaFileWithout_Canis.list")
  os.remove("fastaFileWithout_Mus.list")
  os.remove("fastaFileWithout_Homo.list")

  #print("------------------------------------------------------------------------")
 
  ### cmd4 et cmd5 ###
  cptSeqWithoutCoreTaxon = 0
  for tax in "Canis", "Mus", "Homo" :
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
  fastaFileWithoutCanis = getGeneIDListInBashCommandFile("fastaFileWithout_Canis.list")
  fastaFileWithoutMus = getGeneIDListInBashCommandFile("fastaFileWithout_Mus.list")
  fastaFileWithoutHomo = getGeneIDListInBashCommandFile("fastaFileWithout_Homo.list")
  fastaFileWithoutCoreTaxons = fastaFileWithoutCanis + fastaFileWithoutMus + fastaFileWithoutHomo
  '''
  print ("CDS Files (GeneID) without Canis ("+str(len(fastaFileWithoutCanis))+" ) : \n",fastaFileWithoutCanis,"\n")
  print ("CDS Files (GeneID) without Mus ("+str(len(fastaFileWithoutMus))+" ) : \n",fastaFileWithoutMus,"\n")
  print ("CDS Files (GeneID) without Homo ("+str(len(fastaFileWithoutHomo))+" ) : \n",fastaFileWithoutHomo,"\n")
  print ("Total : ",str(len(fastaFileWithoutCoreTaxons)))
  '''
  os.remove("fastaFileWithout_Canis.list")
  os.remove("fastaFileWithout_Mus.list")
  os.remove("fastaFileWithout_Homo.list")

#------------------------------------------------------------------------------
'''
  REGEX : extract geneID in bash command
  (cf cmd3 in checkDownloadedSequences()  > file.list)
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
  REGEX : extract taxonID in header sequence
  Used in differenceMaker function
'''
def getTaxonID(seqHeader):
  pattern = re.search("(.+)(taxonID=)([^\,\]]+)([\,\]])",seqHeader)    
  if pattern == None :
    return (None)
  taxonID = pattern.group(3)

  return (taxonID)
  

#------------------------------------------------------------------------------
'''
  Make dictionnary for each sequence fastafile type (k = geneFile , v = sequences)
'''

def parseFastaFile (outputFolder, type):

  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,type)

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
  A REVOIR

  Make codon stop correction in fastafile RNAtoCDS
'''
def correctedFastaFile (outputFastaFolder, type, targetGeneID):

  listFastaFiles = makeListOfFastaFile(outputFastaFolder,type)
  for fastaFile in listFastaFiles :
    ext = fastaFile.split("_")
    geneID = ext[0]
    targetFastaFile = outputFastaFolder+os.path.sep+geneID+"_"+type+".fasta"
    if geneID == targetGeneID:

      # tentative 1 ne fonctionne pas
      with open (targetFastaFile) as fastaFile :
        for line in fastaFile :
          if line.startswith(">"):
            continue
          else :
            if line.count("*") != 0 :
              print (geneID)
              print (line)
              line = line.replace("*","X")
              print (line)
              
      
      ''' tentative 2
      with open (targetFastaFile, 'r') as fastaFile :
        contenu = fastaFile.read()
      fastaFile.close

      lignes_contenu = contenu.split("\n")
      for line in lignes_contenu :
        if line.startswith(">"):
          continue
        else :
          if line.count("*") !=0 :
            line = line.replace("*","X")
      contenu = "\n".join(lignes_contenu)

      with open (targetFastaFile, 'w') as fastaFile :
        fastaFile.write(contenu)
      '''


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

  '''
  # Création liste de taxon présents dans le fichier fasta RNAtoCDS 
  currentRNAtoCDSTaxonList = []
  for header in listRNAtoCDSHeader :
    currentTaxonID = getTaxonID(header)
    currentRNAtoCDSTaxonList.append(currentTaxonID)
  #print (currentRNAtoCDSTaxonList)

  # Création liste de taxons présents dans le fichier fasta Protein 
  currentProtTaxonList = []
  for header in listProteinHeader :
    currentTaxonID = getTaxonID(header)
    currentProtTaxonList.append(currentTaxonID)
  #print (currentProtTaxonList)


  # Création liste de taxons absents dans les fichiers fasta RNAtoCDS et Protein
  taxonIDMissingInRNAtoCDSFile = []
  taxonIDMissingInProtFile = []
  for taxonID in orthomamTaxonIDList :
    if taxonID not in currentRNAtoCDSTaxonList :
      taxonIDMissingInRNAtoCDSFile.append(taxonID)
    if taxonID not in currentProtTaxonList :
      taxonIDMissingInProtFile.append(taxonID)
  #print (taxonIDMissingInRNAtoCDSFile)
  '''

  # Comparaisons
  #distance, ratio = levCalclulate(listRNAtoCDSSequence[cptRNAtoCDS],listProteinSequence[cptProt])
  

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
  
  '''
  for taxonID in orthomamTaxonIDList :
    if taxonPresenceInRNAtoCDSFile[taxonID][0] == False:
      print("---\n",taxonPresenceInRNAtoCDSFile,"\n", taxonPresenceInProtFile,"\n---\n")
  '''
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
  
  return (positions)

#------------------------------------------------------------------------------
'''
  Compare RNAtoCDS sequence and AA sequences per gene File.
    Used methods : translateMACSE, makeListOfFastaFile, parseFastaFile 
'''

def checkDownloadedSequences(outputFolder):

  # Translate CDS fasta file (*_NT.fasta) with MACSE (create *_NTtranslate.fasta files) #####
  #translateCDSMACSE(outputFolder)


  # Translate RNAtoCDS fasta file (*_RNAtoCDS.fasta) with MACSE (create *_RNAtoCDStranslate.fasta files) #####
  #translateRNAtoCDSMACSE(outputFolder)

  # Bash primary stats (Homo, Mus, Canis)
  bashCommands (outputFolder)

  # Open GCF.list for create two lists for taxons (Name and ID)
  GCFlistPath = outputFolder+os.path.sep+"GCF"+os.path.sep+"GCF.list"
  orthomamTaxonNameList = []
  orthomamTaxonIDList = []

  with open (GCFlistPath,'r') as GCFlist :
    for line in GCFlist :
      parse = line.strip().split("\t")
      orthomamTaxonNameList.append(parse[0])
      orthomamTaxonIDList.append(parse[-1])
  
  print("\n////////////////////////////////////////////////////////////////////////\n")

  # Open 'taxonForOrthoMaM_v11.list' for create a third list for the classification (based on orthomamTaxonIDList order)
  print ("Récupération des taxons (Noms, ID, Classif) :\n")
  taxonIDtoClassif ={}
  orthomamTaxonClassifList = []

  with open ("taxonForOrthoMaM_v11.list",'r') as classif :
    for line in classif :
      parse = line.strip().split("\t")
      taxonID = parse[1]
      taxonClassif = parse[2]
      taxonIDtoClassif[taxonID]=taxonClassif
    for taxonID in orthomamTaxonIDList :
      orthomamTaxonClassifList = taxonIDtoClassif[taxonID]
  for i in range(0,len(orthomamTaxonNameList)):
    print (orthomamTaxonNameList[i]," / ",orthomamTaxonIDList[i]," / ",orthomamTaxonClassifList)


  # Dictionnary with '$' behind each sequences (AVEC DICO SEQ ET DICO HEADER)
  dicoCDS_WithDollar, dicoCDS_header_WithDollar = parseFastaFile(outputFolder, "NT")
  dicoProt_WithDollar, dicoProt_header_WithDollar = parseFastaFile(outputFolder, "AA")
  dicoRNA_WithDollar, dicoRNA_header_WithDollar = parseFastaFile(outputFolder, "RNA")
  dicoRNAtoCDS_WithDollar, dicoRNAtoCDS_header_WithDollar = parseFastaFile(outputFolder, "RNAtoCDS")
  dicoTranslatedRNAtoCDS_WithDollar, dicoTranslatedRNAtoCDS_header_WithDollar = parseFastaFile(outputFolder, "RNAtoCDStranslate")

  print("\n////////////////////////////////////////////////////////////////////////\n")

  # Stats : differenceMaker() function (line 255)
  
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


  print("Parsing des fichiers et utilisation des dictionnaires (vérifications/stats) :\n")
  # Absence taxons dans les fichier RNAtoCDS et AA :
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
    print ("Taxon ",taxonID," :  absent dans ",len(absenceRNAtoCDS)," fichiers RNAtoCDS et dans ",len(absenceProt)," fichiers Prot")
  print("------------------------------------------------------------------------")

  # Présence taxons dans les fichier RNAtoCDS et AA :
  for taxonID in orthomamTaxonIDList :
    presenceRNAtoCDS = []
    presenceProt = []
    for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      extractPresenceProt = geneIDtoProtFilesStats[geneID][taxonID][0]
      if extractPresenceProt == True :
        presenceRNAtoCDS.append(geneID)
      if extractPresenceProt == True :
        presenceProt.append(geneID)
    print ("Taxon ",taxonID," :  présent dans ",len(presenceRNAtoCDS)," fichiers RNAtoCDS et dans ",len(presenceProt)," fichiers Prot")
  print("------------------------------------------------------------------------")

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
    print ("Taxon ",taxonID," :  a necessité une correction dans ",len(geneIDWithCorrection)," fichiers RNAtoCDS ")
    print ("Taxon ",taxonID," :  n'a pas necessité une correction dans ",len(geneIDWithoutCorrection)," fichiers RNAtoCDS ")
    print("------")
    #print ("Total : ",len(geneIDWithCorrection)+len(geneIDWithoutCorrection))
  print("------------------------------------------------------------------------")

  # Taxons :
  for taxonID in orthomamTaxonIDList :
    geneIDWithMutationOnFirstAminoAcid =[]
    for geneID in dicoTranslatedRNAtoCDS_WithDollar.keys():
      extractPresenceRNAtoCDS = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][0]
      if extractPresenceRNAtoCDS == True :
        extractGeneIDWithDifference = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][2]
        if extractGeneIDWithDifference == False :
          extractgeneIDWithMutationOnFirstAminoAcid = geneIDtoRNAtoCDSFilesStats[geneID][taxonID][3]
          if len(extractgeneIDWithMutationOnFirstAminoAcid) == 1:
            if extractgeneIDWithMutationOnFirstAminoAcid[0] == 1 :
              geneIDWithMutationOnFirstAminoAcid.append(geneID)
    
    print ("Taxon ",taxonID," :  a un changement d'acide aminé en première position dans ",len(geneIDWithMutationOnFirstAminoAcid)," fichiers RNAtoCDS ")
  print("------------------------------------------------------------------------")



  print("\n////////////////////////////////////////////////////////////////////////\n")

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
