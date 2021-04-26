#!/usr/bin/env python3
#_*_ coding:utf-8 _*_

def verifSymbol():
    with open ("maListeOtho.tsv", 'r') as Ortho, open ("gene3.list", 'r') as geneList :
        listeGenesSymbols = []

        # Transo fichier gene en liste de symbole de CDS
        for G in geneList :
            parseG = G.split("\t")
            if parseG[0] == "CDS" :
                listeGenesSymbols.append(parseG[1].upper())

        # Transfo fichier ortho en listes
        nb_ligne = -1
        listeOrtho_taxon = []
        listeOrtho_symbol = []
        listeOrtho_id = []

        for O in Ortho :
            nb_ligne +=1
            parseO = O.split("\t")
            listeOrtho_taxon.append(parseO[0].upper())
            listeOrtho_symbol.append(parseO[1].upper())
            listeOrtho_id.append(parseO[2])


        listeSymbolesDifferentsDeLaListeInitiale = []
        listeTaxons_SymbolesDifferentsDeLaListeInitiale = []

        listeSymboleDatasets_HomoSapiens = []
        listeIDDatasets_HomoSapiens = []
        listeTaxonDatasets_HomoSapiens = []

        listeSymboleDatasets_MusMusculus = []
        listeIDDatasets_MusMusculus = []

        listeSymboleDatasets_CanisLupus = []
        listeIDDatasets_CanisLupus = []

        # Liste des gènes recherchés avec Datasets
        for i in range (0,nb_ligne) :
            if (listeOrtho_symbol[i] in listeGenesSymbols):
                #if (listeOrtho_taxon[i] == "Homo sapiens") :
                if listeOrtho_taxon[i].startswith("H"):
                    listeSymboleDatasets_HomoSapiens.append(listeOrtho_symbol[i])
                    listeIDDatasets_HomoSapiens.append(listeOrtho_id[i])

                #if (listeOrtho_taxon[i] == "Mus musculus") :
                if listeOrtho_taxon[i].startswith("M"):
                    listeSymboleDatasets_MusMusculus.append(listeOrtho_symbol[i])
                    listeIDDatasets_MusMusculus.append(listeOrtho_id[i])

                #if (listeOrtho_taxon[i] == "Canis lupus familiaris") :
                if listeOrtho_taxon[i].startswith("C"):
                    listeSymboleDatasets_CanisLupus.append(listeOrtho_symbol[i])
                    listeIDDatasets_CanisLupus.append(listeOrtho_id[i])
            else :
                listeSymbolesDifferentsDeLaListeInitiale.append(listeOrtho_symbol[i])
                listeTaxons_SymbolesDifferentsDeLaListeInitiale.append(listeOrtho_taxon[i])


       # Liste des gènes NON-recherchés avec Datasets (si pb de symbole non connu par NCBI ou pb de fichier zip)
        liste_CDS_rechercheDatasets = []
        listeGeneNonRechercheParDatasets = []
        
        for gene in listeGenesSymbols :
            #print (gene)
            if gene in listeOrtho_symbol :
                liste_CDS_rechercheDatasets.append(gene)
            else :
                #print ("non trouvé")
                listeGeneNonRechercheParDatasets.append(gene)


        # Listes des gènes en communs pour :

            # Homo (seul)
        listeSymbol_Homo = []
            # Homo_Mus_Canis
        listeSymbol_Homo_Mus_Canis = []
            # Homo_Mus (non Canis)
        listeSymbol_Homo_Mus = []
            # Homo_Canis (non Mus)
        listeSymbol_Homo_Canis = []

        for gene in listeSymboleDatasets_HomoSapiens :
            if (gene not in listeSymboleDatasets_MusMusculus) and (gene not in listeSymboleDatasets_CanisLupus) :
                listeSymbol_Homo.append(gene)            
            if (gene in listeSymboleDatasets_MusMusculus) and (gene in listeSymboleDatasets_CanisLupus) :
                listeSymbol_Homo_Mus_Canis.append(gene)
            if (gene in listeSymboleDatasets_MusMusculus) and (gene not in listeSymboleDatasets_CanisLupus):
                listeSymbol_Homo_Mus.append(gene)
            if (gene in listeSymboleDatasets_CanisLupus) and (gene not in listeSymboleDatasets_MusMusculus) :
                listeSymbol_Homo_Canis.append(gene)

            # Canis (seul)
        listeSymbol_Canis = []
            # Mus_Canis (non Homo)
        listeSymbol_Mus_Canis = []

        for gene in listeSymboleDatasets_CanisLupus :
            if (gene in listeSymboleDatasets_MusMusculus) and (gene not in listeSymboleDatasets_HomoSapiens) :
                listeSymbol_Mus_Canis.append(gene)
            if (gene not in listeSymboleDatasets_MusMusculus) and (gene not in listeSymboleDatasets_HomoSapiens) :
                listeSymbol_Canis.append(gene)

            # Mus (seul)
        listeSymbol_Mus = []

        for gene in listeSymboleDatasets_MusMusculus :
            if (gene not in listeSymboleDatasets_CanisLupus) and (gene not in listeSymboleDatasets_HomoSapiens) :
                listeSymbol_Mus.append(gene)


        # Récap
        nb_CDS_genome = len(listeGenesSymbols)
        nb_CDS_rechercheDatasets = len(liste_CDS_rechercheDatasets)
        nb_cds_NONrechercheDatasets = len(listeGeneNonRechercheParDatasets)
        
        print("Nombre de gènes (CDS) du génome humain (pour requête): ",nb_CDS_genome)
        print("Nombre de gènes recherchés par datasets : ",nb_CDS_rechercheDatasets)
        print("Nombre de gènes NON recherchés par datasets : ",nb_cds_NONrechercheDatasets)

        print("\n")
        print("Nombre de gène pour Homo Sapiens : ",len(listeSymboleDatasets_HomoSapiens))
        print("Nombre de gène pour Mus Musculus : ",len(listeSymboleDatasets_MusMusculus))
        print("Nombre de gène pour Canis lupus familiaris : ",len(listeSymboleDatasets_CanisLupus))

        print("\n")
        print("Nombre de symboles différents de la liste initiale : ", len(listeSymbolesDifferentsDeLaListeInitiale))
        #print(listeSymbolesDifferentsDeLaListeInitiale)
        #print(listeTaxons_SymbolesDifferentsDeLaListeInitiale)

        print("\n")
        print("Valeurs pour Diagramme de Venn :")
        # Valeurs à partir de la liste des gènes de Homo
        print("Nombre de symboles Homo Sapiens seul (non Mus et non Canis) : ",len(listeSymbol_Homo))
        print("Nombre de symboles Mus Musculus seul (non Homo et non Canis) : ",len(listeSymbol_Mus))
        print("Nombre de symboles Canis seul (non Mus et non Homo) : ",len(listeSymbol_Canis))

        print("Nombre de symboles communs Homo, Mus (non Canis) : ",len(listeSymbol_Homo_Mus))
        print("Nombre de symboles communs Homo, Canis (non Mus) : ",len(listeSymbol_Homo_Canis))
        print("Nombre de symboles communs Mus, Canis (non Homo) : ",len(listeSymbol_Mus_Canis))

        print("Nombre de symboles communs Homo Mus et Canis : ",len(listeSymbol_Homo_Mus_Canis))



################## MAIN ##################


#listesHumanGenomeNCBI()
verifSymbol()