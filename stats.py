#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

from statistics import *
#from decimal import Decimal
import numpy
import sys


#with open ("res.canis", 'r') as canis, open ("res.mus", 'r') as mus, open ("res.choloepus", 'r') as choloepus, open ("res.loxodonta", 'r') as loxodonta :

def usage() :
    print ("Ajouter fichier pour analyse (ex: res.canis")
    # améliorer

taxon = sys.argv[1]

if len(sys.argv) != 2 :
    usage()

def cdsStats(taxon):
    fileName = taxon+"_nbCDS.stat"
    liste_cpt_CDS = []
    with open (fileName, 'r') as statFile :
        for line in statFile :
            parse = line.strip().split(":")
            val = int(parse[1])
            liste_cpt_CDS.append(val)
        print ("CDS par gène pour",taxon," : ")
        print ("   Moyenne : ",(numpy.round(mean(liste_cpt_CDS),2)),"| Max : ",max(liste_cpt_CDS),"| Min :",min(liste_cpt_CDS),"| Ecart type : ",numpy.round(pstdev(liste_cpt_CDS),2))


cdsStats(taxon)