#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

import sys, os, re

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
  MACSE translate fasta files
'''

# CHEMINS A REVOIR SELON UTILISATION !!!! (chemin MACSE et répertoires scripts et d'enregistrementdes fastas)

def translateMACSE(outputFolder):
  path = os.getcwd()
  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,"NT")

  for fastaFile in listFastaFiles :
    ext = fastaFile.split("_")
    geneID = ext[0]
    pathFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+fastaFile
    pathOutputFile = path+os.path.sep+fastaFilesDirectory+os.path.sep+geneID+"_AA_translate.fasta"
    cmdMACSE = "java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA -out_AA " + pathOutputFile + " -seq " + pathFile
    os.system(cmdMACSE)

''' exemple
java -jar /home/lenovo/Documents/Logiciels/macse_v2.05.jar -prog translateNT2AA 
  -out_AA /home/lenovo/Documents/Master/Stage/Projet_Stage/Etape1_AvecVincent_Finale/testDu210517/FASTA/10000_AA_translate.fasta 
  -seq /home/lenovo/Documents/Master/Stage/Projet_Stage/Etape1_AvecVincent_Finale/testDu210517/FASTA/10000_NT.fasta

'''

#------------------------------------------------------------------------------
'''
  Make dictionnary for NT sequence and AA sequences (k = geneFile , v = sequences)
'''

def parseFastaFile (outputFolder, type):

  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"
  listFastaFiles = makeListOfFastaFile(fastaFilesDirectory,type)

  seqDictionary = {}
  cptSeqDictionary = {}
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
              # Version dico de seq concaténées
              seqDictionary[geneID]=seq
              cptSeqDictionary[geneID]=1
            else :
              #seqDictionary[geneID].append(seq)
              seqDictionary[geneID]+=seq
              cptSeqDictionary[geneID]+=1
          seq = ""
          first = False
        else :
          seq += line.strip()
      if first == False :
        #seqDictionary[geneID].append(seq)
        seqDictionary[geneID]+=seq
        cptSeqDictionary[geneID]+=1

  return (seqDictionary,cptSeqDictionary)

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

  ##### Translate CDS fasta folder with MACSE #####
  translateMACSE(outputFolder)

  ##### dictionary of CDS standard Translated (MACSE) #####
  
  dicoCDS, dicoCptSeqCDS = parseFastaFile(outputFolder, "translate")
  dicoTranslatedCDS, dicoCptSeqTranslatedCDS = parseFastaFile(outputFolder, "translate")
  dicoProt, dicoCptSeqProt = parseFastaFile(outputFolder, "AA")

  listGenesOk = []
  listGenesNotOk = []
  for gene in dicoTranslatedCDS.keys() :
    if dicoTranslatedCDS[gene] == dicoProt[gene] :
      listGenesOk.append(gene)
    else :
      listGenesNotOk.append(gene)
  print ("Genes with good CDS sequences vs proteins references : \n", str(len(listGenesOk)))
  print ("Genes with wrong CDS sequences vs proteins references (or wrong translation) : \n",str(len(listGenesNotOk)),"\n")
  print ("GeneID with bad translation :")
  print (listGenesNotOk)
  
  ##### bash commands #####
  # voir pour amélioration / ajustement chemin ?
    #path = os.getcwd()
    #fastaFilesDirectory = path+os.path.sep+outputFolder+os.path.sep+"FASTA"+os.path.sep

  fastaFilesDirectory = outputFolder+os.path.sep+"FASTA"+os.path.sep

  print("------------------------------------------------------------------------")
  ### cmd0 et cmd1 ###
  cmd0 = "ls "+fastaFilesDirectory+"*_NT.fasta | wc -l"
  print (cmd0)
  os.system(cmd0)
  for tax in "Canis", "Mus", "Homo" :
    cmd1 = "cat "+fastaFilesDirectory+"*_NT.fasta | grep -c "+tax
    print (cmd1)
    os.system(cmd1)

  print("------------------------------------------------------------------------")
  ### cmd2 et cmd3 ###
  cptSeqWithoutCoreTaxon = 0
  for tax in "Canis", "Mus", "Homo" :
    cmd2 = "grep -c "+tax+" "+fastaFilesDirectory+"*_NT.fasta | grep -c :0 "
    print (cmd2)
    os.system(cmd2)
    cmd3 = "grep -c "+tax+" "+fastaFilesDirectory+"*_NT.fasta| grep :0 > "+"fastaFileWithout_"+tax+".list"
    os.system(cmd3)
    rescmd2 = os.popen(cmd2).readlines()
    cptSeqWithoutCoreTaxon += int(rescmd2[0])
  print ("Total : ", str(cptSeqWithoutCoreTaxon))

  print("------------------------------------------------------------------------")
  ### Récupération informations de cmd3 ###
  fastaFileWithoutCanis = getGeneIDListInBashCommandFile("fastaFileWithout_Canis.list")
  print ("GeneID without Canis ("+str(len(fastaFileWithoutCanis))+" ) : \n",fastaFileWithoutCanis,"\n")
  fastaFileWithoutMus = getGeneIDListInBashCommandFile("fastaFileWithout_Mus.list")
  print ("GeneID without Mus ("+str(len(fastaFileWithoutMus))+" ) : \n",fastaFileWithoutMus,"\n")
  fastaFileWithoutCoreTaxons = fastaFileWithoutCanis + fastaFileWithoutMus
  print ("Total : ",str(len(fastaFileWithoutCoreTaxons)))

  os.remove("fastaFileWithout_Canis.list")
  os.remove("fastaFileWithout_Mus.list")
  os.remove("fastaFileWithout_Homo.list")

  print("------------------------------------------------------------------------")
 
  ##### Créations fichiers de sorties ##### (en cours)
  
  # for gene in listGenesOk :
  #    if gene not in fastaFileWithoutCoreTaxons.keys():
  #       with open ...
  #         gene  dicoCptSeqCDS[gene]   dicoCDS[gene]


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
  Main
'''

if __name__ == "__main__" :
  if len (sys.argv) != 2:
      usage()
      sys.exit()
  else :
    outputFolder = sys.argv[1]  
    checkDownloadedSequences(outputFolder)
