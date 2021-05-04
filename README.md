# Stage-OrthoMaM
_____________________________________________________________________________________________________________________________________________________________
Etape1

addGeneID.py :
Récupère et génère une liste de gene-id à partir d’une liste de taxons.
En cours d'amélioration. 

addGenomes.py :
Récupère et ajoute dans un dossier les séquences fasta (en fonction de la liste de gene-id de la fonction addGeneID.py)
En cours d'amélioration. 

orthologParse.py :
Récupère des listes de gene_id et des orthologues 1:1 par taxons cherchés (avec leur gene_id = other_gene_id dans le fichier d’entrée gene-ortholog.gz)

maListeOrtho.py : Test stratégie DATASETS
Nécessite dans le répertoire courant :
- datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/quickstarts/command-line-tools/)
- gene3.list (https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/ puis commande "tail -n+2 GCF_000001405.39_GRCh38.p13_feature_table.txt | cut -f1,15,16,20 |sort -u > gene3.list")
_____________________________________________________________________________________________________________________________________________________________
