# Stage-OrthoMaM

_____________________________________________________________________________________________________________________________________________________________
Etape1

Main file :

	buildOrtholog.py
 
Secondary files (imported) :
  
    selectHumanGeneID.py
    
    downloadRelevantGCF.py
    
    addGenomes.py
    
    checkSequences.py

How to run it :

./buildOrtholog.py assembly_summary.txt core_species.list gene_orthologs orthologFasta

	summaryFile = sys.argv[1]	        #  assembly_summary.txt (to be downloded from the NCBI, see below)
	
	coreTaxonList = sys.argv[2]	      	#  core_species.list (list of the core taxa)
	
	orthologFile = sys.argv[3]		#  gene_orthologs (to be downloded from the NCBI, see below)
	
	orthologFasta = sys.argv[4]	  	#  working directory

_____________________________________________________________________________________________________________________________________________________________
Fonction usage du module buildOrthologs :
    
    This program build ortholog fasta files of orthologous genes using the human gene identifier as a cross reference and three core taxa.
    

    - assembly_summary.tsv : is supposed to have 22 fields including :
        assembly_accession  refseq_category taxid   organisme_name  ftp_path (to download the GCF files)
        (following NCBI convention: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/ assembly_summary.txt)
        WARNING : REMEMBER TO DELETE THE HYBRID TAXON (30522) FROM THE FILE (Bos indicus x Bos taurus)

    - core_species.list : should contain one taxon id per line, for a core made of Homo Sapiens, Mus musculus and Canis Lupus
    familiaris it will be (human should always be first):
        9606        (Homo sapiens)
        10090       (Mus musculus)
        9615        (Canis Lupus familiaris)
	
    - gene_orthologs.tsv : is supposed to have 5 fields and to contain only 1:1 ortholgs:
        tax_id GeneID  relationship    Other_tax_id    Other_GeneID 
        (following NCBI convention: http://ftp.ncbi.nlm.nih.gov/gene/DATA/ gene-ortholog.gz )



    Usefull links : 
      - https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice
      - https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.ncbi_s_annotation_displayed_on

_____________________________________________________________________________________________________________________________________________________________
