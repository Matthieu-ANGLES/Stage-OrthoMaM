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

''' exemple
java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA 
  -out_AA /home/lenovo/Documents/Master/Stage/Projet_Stage/Etape1_AvecVincent_Finale/testDu210517/FASTA/10000_AA_translate.fasta 
  -seq /home/lenovo/Documents/Master/Stage/Projet_Stage/Etape1_AvecVincent_Finale/testDu210517/FASTA/10000_NT.fasta

'''
#------------------------------------------------------------------------------
'''
  MACSE translate RNAtoCDS fasta files
'''

# CHEMINS A REVOIR SELON UTILISATION !!!! (chemin MACSE et répertoires scripts et d'enregistrementdes fastas)

def translateRNAMACSE(outputFolder):
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

''' exemple
java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA 
  -out_AA /home/lenovo/Documents/Master/Stage/Projet_Stage/Etape1_AvecVincent_Finale/testDu210517/FASTA/10000_AA_translate.fasta 
  -seq /home/lenovo/Documents/Master/Stage/Projet_Stage/Etape1_AvecVincent_Finale/testDu210517/FASTA/10000_NT.fasta
'''

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
  Make dictionnary for NT sequence and AA sequences (k = geneFile , v = sequences)
'''

def parseFastaFile (outputFolder, type):

  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,type)
  #print (listFastaFiles)

  seqDictionary = {}
  for fastaFile in listFastaFiles :
    ext = fastaFile.split("_")
    geneID = ext[0]
    pathFile = fastaFilesDirectory+os.path.sep+fastaFile
    with open (pathFile,'r') as file :
      seq = ""
      first = True
      for line in file :
        if line.startswith(">") :
          if first == False :
            if geneID not in seqDictionary.keys():
              # Version dico de listes de seq
              #seqDictionary[geneID]=[seq]
              # Version dico de seq concaténées (seq+"$")
              seqDictionary[geneID]=seq+"$"
            else :
              #seqDictionary[geneID].append(seq)
              seqDictionary[geneID]+=seq+"$"
          seq = ""
          first = False
        else :
          seq += line.strip()
      if first == False :
        #------------------
        # ajout:
        if geneID not in seqDictionary.keys():
          #seqDictionary[geneID]=[seq]
          seqDictionary[geneID]=seq+"$"

        #------------------
        else :
          #seqDictionary[geneID].append(seq)
          seqDictionary[geneID]+=seq+"$"

  return (seqDictionary)

#------------------------------------------------------------------------------
'''
  Compare NT sequence and AA sequences per gene.
    Used methods : translateMACSE, makeListOfFastaFile, parseFastaFile 
    Make lists :
    - geneID with good translation = listGenesOk
    - geneID with bad translation = listGenesNotOk
    - geneID with not all cores taxons = fastaFileWithoutCoreTaxons
'''
############## VOIR POUR CREER DES set() au lieux des listes ?? ##############


def checkDownloadedSequences(outputFolder):

  ##### Translate CDS fasta file (*_NT.fasta) with MACSE (create *_AA_translate.fasta files) #####
  translateCDSMACSE(outputFolder)
  translateRNAMACSE(outputFolder)

  ##### Dictionary of fastafiles Translated (MACSE) #####
  
  #Partie pour dico de listes
  '''
  dicoCDS = parseFastaFile(outputFolder, "NT")
  dicoTranslatedCDS = parseFastaFile(outputFolder, "NTtranslate")
  dicoProt = parseFastaFile(outputFolder, "AA")
  #dicoRNA = parseFastaFile(outputFolder, "RNA")
  #dicoRNAtoCDS = parseFastaFile(outputFolder, "RNAtoCDS")
  dicoTranslatedRNAtoCDS = parseFastaFile(outputFolder, "RNAtoCDStranslate")
  #dicoTranslated = parseFastaFile(outputFolder, "CDSTranslated")
  '''

  # Partie pour dico avec $
  
  # Dictionnary with '$' behind each sequences
  dicoCDS_WithDollar = parseFastaFile(outputFolder, "NT")
  dicoTranslatedCDS_WithDollar = parseFastaFile(outputFolder, "NTtranslate")
  dicoProt_WithDollar = parseFastaFile(outputFolder, "AA")
  dicoRNA_WithDollar = parseFastaFile(outputFolder, "RNA")
  dicoTranslatedRNAtoCDS_WithDollar = parseFastaFile(outputFolder, "RNAtoCDStranslate")
  #print (dicoTranslatedRNAtoCDS_WithDollar)

  # Copy of the dictionaries with removal of the dollars
  dicoCDS = {}
  dicoTranslatedCDS = {}
  dicoProt = {}
  dicoRNA = {}
  dicoTranslatedRNAtoCDS = {}

  for gene in dicoCDS_WithDollar.keys() :
    #print (gene)
    #print (dicoCDS_WithDollar[gene])
    dicoCDS[gene] = dicoCDS_WithDollar[gene].replace("$","")
    dicoTranslatedCDS[gene] = dicoTranslatedCDS_WithDollar[gene].replace("$","")
    dicoProt[gene] = dicoProt_WithDollar[gene].replace("$","")
    
    if gene in dicoTranslatedRNAtoCDS_WithDollar.keys():
      dicoTranslatedRNAtoCDS[gene] = dicoTranslatedRNAtoCDS_WithDollar[gene].replace("$","")
      dicoRNA[gene] = dicoRNA_WithDollar[gene].replace("$","")
  
  
  print("\n////////////////////////////////////////////////////////////////////////\n")

  # Count Good or Bad Translations :
  listGenesOk = []
  listGenesNotOk = []
  listGeneNotInDicoRNAtoCDS = []
  for gene in dicoTranslatedCDS.keys() :
    if gene in dicoTranslatedRNAtoCDS.keys():
      if dicoTranslatedRNAtoCDS[gene] == dicoProt[gene] :
        listGenesOk.append(gene)
      else :
        listGenesNotOk.append(gene)
    else :
      listGeneNotInDicoRNAtoCDS.append(gene)
    
  print ("Genes with good CDS sequences vs proteins references : \n", str(len(listGenesOk)))
  print ("Genes with wrong CDS sequences vs proteins references (or wrong translation) : \n",str(len(listGenesNotOk)))

  print ("Genes not in original CDS fasta file : ",str(len(listGeneNotInDicoRNAtoCDS)),"\n")

  print("------------------------------------------------------------------------")

  # Retrait des codons stop et ajout de X à la place 
  geneCorrectionStopCodonList_RNAtoCDS = []
  ''' un seul changement par seq ?
  for gene in listGenesNotOk :
    if gene in dicoTranslatedRNAtoCDS.keys():
      if dicoTranslatedRNAtoCDS[gene].count("*") != 0 :
        currentSeq = dicoTranslatedRNAtoCDS[gene].replace("*","X")
        if currentSeq == dicoProt[gene] :
          listGenesOk.append(gene)
          listGenesNotOk.remove(gene)
          geneCorrectionStopCodonList_RNAtoCDS.append(gene)
  '''
  seqDistanceEgalOne = []
  seqdistanceEgalTwo = []
  seqDistanceEgalThree = []
  seqDistanceEgalFour = []
  seqDistanceSupEgalFive = []

  seqRatioSup98 = []
  seqRatioSup95Inf98 = []
  seqRatiosup70Inf95 = []
  seqRatioSup60Inf70 = []
  seqRatioInf60 = []
  
  for gene in listGenesNotOk :
    if gene in dicoTranslatedRNAtoCDS.keys():
      if dicoTranslatedRNAtoCDS[gene].count("*") != 0 :
        dicoTranslatedRNAtoCDS[gene] = dicoTranslatedRNAtoCDS[gene].replace("*","X")

      distance, ratio = levCalclulate(dicoTranslatedRNAtoCDS[gene],dicoProt[gene])
      if dicoTranslatedRNAtoCDS[gene] == dicoProt[gene] :
        #listGenesOk.append(gene)
        #listGenesNotOk.remove(gene)
        geneCorrectionStopCodonList_RNAtoCDS.append(gene)

      else : 
        if distance == 1 :
          seqDistanceEgalOne.append(gene)
        if distance == 2 :
          seqdistanceEgalTwo.append(gene)
        if distance == 3 :
          seqDistanceEgalThree.append(gene)
        if distance == 4 :
          seqDistanceEgalFour.append(gene)
        if distance >= 5 :
          seqDistanceSupEgalFive.append(gene)

        if (ratio >= 0.98) :
          seqRatioSup98.append(gene)
        if (ratio >= 0.95) and (ratio < 0.98) :
          seqRatioSup95Inf98.append(gene)
        if (ratio >= 0.70) and (ratio < 0.95) :
          seqRatiosup70Inf95.append(gene)
        if (ratio >= 0.60) and (ratio < 0.70) :
          seqRatioSup60Inf70.append(gene)
        if (ratio < 0.60) :
          seqRatioInf60.append(gene)
  
  for gene in geneCorrectionStopCodonList_RNAtoCDS :
    listGenesOk.append(gene)
    listGenesNotOk.remove(gene)

  print ("After replace ' * ' with ' X ' (in RNAtoCDS too):")
  print ("Genes with good CDS sequences vs proteins references : \n", str(len(listGenesOk)))
  print ("Genes with wrong CDS sequences vs proteins references (or wrong translation) : \n",str(len(listGenesNotOk)))
  print ("Concerned genes : ",str(len(geneCorrectionStopCodonList_RNAtoCDS)),"\n")


  print("------------------------------------------------------------------------")

  
  print ("number of sequences with a distance = 1 : ",str(len(seqDistanceEgalOne)),"\n",seqDistanceEgalOne,"\n")
  print ("number of sequences with a distance = 2 : ",str(len(seqdistanceEgalTwo)),"\n",seqdistanceEgalTwo,"\n")
  print ("number of sequences with a distance = 3 : ",str(len(seqDistanceEgalThree)),"\n",seqDistanceEgalThree,"\n")
  print ("number of sequences with a distance = 4 : ",str(len(seqDistanceEgalFour)),"\n",seqDistanceEgalFour,"\n")
  print ("number of sequences with a distance >= 5 : ",str(len(seqDistanceSupEgalFive)),"\n",seqDistanceSupEgalFive,"\n")

  print("------------------------------------------------------------------------")

  print ("number of sequences with a ratio >= 0.98 : ",str(len(seqRatioSup98)),"\n",seqRatioSup98,"\n")
  print ("number of sequences with a 0.95 > ratio < 0.98 : ",str(len(seqRatioSup95Inf98)),"\n",seqRatioSup95Inf98,"\n")
  print ("number of sequences with a 0.70 > ratio < 0.95 : ",str(len(seqRatiosup70Inf95)),"\n",seqRatiosup70Inf95,"\n")
  print ("number of sequences with a 0.60 > ratio < 0.70 : ",str(len(seqRatioSup60Inf70)),"\n",seqRatioSup60Inf70,"\n")
  print ("number of sequences with a ratio < 0.60 : ",str(len(seqRatioInf60)),"\n",seqRatioInf60,"\n")

  print("------------------------------------------------------------------------")
  

  ##### bash commands #####
  # voir pour amélioration / ajustement chemin ?
    #path = os.getcwd()
    #fastaFilesDirectory = path+os.path.sep+outputFolder+os.path.sep+"FASTA"+os.path.sep

  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"+os.path.sep

  print("\n////////////////////////////////////////////////////////////////////////\n")
  
  '''
  ### Contrôle sur les fichiers CDS (NT)
  '''
  
  print ("Check CDS sequences Files (translated with MACSE) :\n")


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
 
  # Vérif distance seq et absence taxon pour un geneID
  '''
  taxName = ["Canis", "Mus", "Homo"]
  canis, mus, homo = 0

  for tax in taxName :
    for gene in fastaFileWithoutCanis :
      if gene in geneDistanceEgalOne :
        print ()
  '''






  ##### Créations fichiers de sorties #####
  '''
  output = outputFolder+os.path.sep+"output_orthoMaM_marker.list"
  with open (output,"w") as orthoMarker :
    for gene in listGenesOk :
      if gene not in fastaFileWithoutCoreTaxons :
        orthoMarker.write(gene+"\t"+str(dicoCDS_WithDollar[gene].count("$"))+"\t"+dicoCDS[gene]+"\n")
        #orthoMarker.write(gene+"\t"+str(len(dicoCDS[gene]))+"\t"+str(dicoCDS[gene])+"\n")
  print ("The file 'output_orthoMaM_marker.list' was created.")

  output = outputFolder+os.path.sep+"output_badTranslated_marker.list"
  with open (output,"w") as badTransMarker :
    for gene in listGenesNotOk :
      if gene not in fastaFileWithoutCoreTaxons :
        badTransMarker.write(gene+"\t"+str(dicoCDS_WithDollar[gene].count("$"))+"\t"+dicoCDS[gene]+"\n")
        #badTransMarker.write(gene+"\t"+str(len(dicoCDS[gene]))+"\t"+str(dicoCDS[gene])+"\n")
  print ("The file 'output_badTranslated_marker.list' was created.")

  output = outputFolder+os.path.sep+"output_notHomoMusCanis_marker.list"
  with open (output,"w") as exludedMarker :
    for gene in fastaFileWithoutCoreTaxons :
      exludedMarker.write(gene+"\t"+str(dicoCDS_WithDollar[gene].count("$"))+"\t"+dicoCDS[gene]+"\n")
      #exludedMarker.write(gene+"\t"+str(len(dicoCDS[gene]))+"\t"+str(dicoCDS[gene])+"\n")
  print ("The file 'output_notHomoMusCanis_marker.list was' created.","\n")
  '''
  
  # Fichier avec geneID et séquences RNAtoCDS versus Proteine
  output = outputFolder+os.path.sep+"output_RNAtoCDS.stats"
  with open (output,"w") as badTransMarker :

    badTransMarker.write("After replace ' * ' with ' X ' : (Concerned genes : " + str(len(geneCorrectionStopCodonList_RNAtoCDS))+")\n")
    badTransMarker.write("Genes with good RNAtoCDS sequences vs proteins references : \n" + str(len(listGenesOk))+"\n")
    badTransMarker.write("Genes with wrong RNAtoCDS sequences vs proteins references (or wrong translation) : \n" + str(len(listGenesNotOk))+"\n\n")
    
    #badTransMarker.write("# List geneID which required stop codon correction * to X : "+str(len(geneCorrectionStopCodonList_RNAtoCDS))+"\n")
    #badTransMarker.write("# "+str(geneCorrectionStopCodonList_RNAtoCDS)+"\n\n")

    badTransMarker.write("----------------------------------------------------------------------------------------------"+"\n")
    badTransMarker.write("# List geneID without core taxons : "+str(len(fastaFileWithoutCoreTaxons))+"\n\n")

    badTransMarker.write("# List geneID without  Homo Sapiens : "+str(len(fastaFileWithoutHomo))+"\n")
    badTransMarker.write("# "+str(fastaFileWithoutHomo)+"\n\n")

    badTransMarker.write("# List geneID without  Mus musculus : "+str(len(fastaFileWithoutMus))+"\n")
    badTransMarker.write("# "+str(fastaFileWithoutMus)+"\n\n")

    badTransMarker.write("# List geneID without  Canis lupus familiaris : "+str(len(fastaFileWithoutCanis))+"\n")
    badTransMarker.write("# "+str(fastaFileWithoutCanis)+"\n")
    badTransMarker.write("----------------------------------------------------------------------------------------------"+"\n\n")
    
    badTransMarker.write("# List geneID for sequences IN CDS files and NOT IN RNA files : "+str(len(listGeneNotInDicoRNAtoCDS))+"\n")
    badTransMarker.write("# "+str(listGeneNotInDicoRNAtoCDS)+"\n\n")
    badTransMarker.write("----------------------------------------------------------------------------------------------"+"\n\n")
    badTransMarker.write("# List geneID with wrong RNAtoCDS sequence vs proteins reference (or wrong translation) : "+str(len(listGenesNotOk))+"\n")
    badTransMarker.write("# "+str(listGenesNotOk)+"\n\n")


    badTransMarker.write("======================================================================================================"+"\n")
    badTransMarker.write("# RNAtoCDS sequences VS Prot sequences for geneID with translations errors :"+str(len(listGenesNotOk))+"\n")
    badTransMarker.write("======================================================================================================"+"\n")

    for gene in listGenesNotOk :
      #print (gene)
      if gene not in fastaFileWithoutCoreTaxons :
        if gene in dicoTranslatedRNAtoCDS.keys():
          if gene not in geneCorrectionStopCodonList_RNAtoCDS:
            badTransMarker.write("/geneID : "+gene+"\n")
            badTransMarker.write("/RNA sequence :"+"\n"+dicoRNA[gene]+"\n")
            badTransMarker.write("/Translated RNAtoCDS sequence :"+"\n"+dicoTranslatedRNAtoCDS[gene]+"\n")
            badTransMarker.write("/Protein sequence :"+"\n"+dicoProt[gene]+"\n")
            badTransMarker.write("----------------------------------------------------------------------------------------------"+"\n")
  
  print ("The file 'output_RNAtoCDS.stats' was created.")




  # TESTS !
  '''
  for gene in listGenesNotOk :
    print (gene, " :")
    print ("CDS :\n",dicoTranslatedCDS[gene])
    print ("----------------------------------------------")
    if gene not in dicoTranslatedRNAtoCDS.keys():
      print ("RNAtoCDS : No sequence")
    else :
      print ("RNAtoCDS :\n",dicoTranslatedRNAtoCDS[gene],"\n")
    print ("----------------------------------------------")
    print ("Prot :\n",dicoProt[gene],"\n")
  '''

  '''
  for gene in listGenesNotOk :
    print (gene, " :")
    print ("CDS :\n",dicoTranslatedCDS_WithDollar[gene])
    print ("----------------------------------------------")
    if gene not in dicoTranslatedRNAtoCDS_WithDollar.keys():
      print ("RNAtoCDS : No sequence")
    else :
      print ("RNAtoCDS :\n",dicoTranslatedRNAtoCDS_WithDollar[gene],"\n")
    print ("----------------------------------------------")
    print ("Prot :\n",dicoProt_WithDollar[gene],"\n")
  '''

  '''
  for gene in RNAtoCDS_notOk_WithProt :
    print (gene, " :")
    print ("CDS :\n",dicoTranslatedCDS[gene])
    print ("----------------------------------------------")
    if gene not in dicoTranslatedRNAtoCDS.keys():
      print ("RNAtoCDS : No sequence")
    else :
      print ("RNAtoCDS :\n",dicoTranslatedRNAtoCDS[gene],"\n")
    print ("----------------------------------------------")
    print ("Prot :\n",dicoProt[gene],"\n")
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
