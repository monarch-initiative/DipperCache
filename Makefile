# should run as
# make  --no-builtin-variables --no-builtin-rules
# since this is not performing any traditional build tasks

DIPPER = $$HOME/Projects/Monarch/dipper
WGET = /usr/bin/wget --timestamping
FULLPTH := --force-directories --no-host-directories
RM := rm --force --recursive --verbose

.PHONY: cruft recent clean dipper/scripts

all:  animalqtldb bgee clinvar ctd flybase genereviews go gwascatalog \
	hgnc hpoa impc kegg  mmrrc  monochrom mpd ncbigene omia  orphanet \
	panther reactome rgd sgd string  wormbase zfin zfinslim

# ommited for cause
	# coriell ensembl eom mgi monarch mychem mychem omim ucscbands

# newish updates
recent:
	find ./*/ -mtime -5 -ls

cruft:
	make -n | sed 's|^mkdir |\nmkdir |g' > $@

##########################################
# animalqtldb
AQTLDL = https://www.animalgenome.org/QTLdb
CDAQTL = cd animalqtldb ;

AQTLGI = gene_info.gz \
		Bos_taurus.gene_info.gz \
		Sus_scrofa.gene_info.gz \
		Gallus_gallus.gene_info.gz \
		Ovis_aries.gene_info.gz \
		Oncorhynchus_mykiss.gene_info.gz \
		Equus_caballus.gene_info.gz

AQTLTMP = QTL_Btau_4.6.gff.txt.gz \
		QTL_EquCab2.0.gff.txt.gz \
		QTL_GG_4.0.gff.txt.gz \
		QTL_OAR_3.1.gff.txt.gz

AQTLVER = pig_QTLdata.txt \
		sheep_QTLdata.txt \
		cattle_QTLdata.txt \
		chicken_QTLdata.txt \
		horse_QTLdata.txt \
		rainbow_trout_QTLdata.txt

animalqtldb: ncbigene animalqtldb/ \
		$(foreach spc, $(AQTLGI), animalqtldb/$(spc)) \
		$(foreach spc, $(AQTLTMP), animalqtldb/$(spc)) \
		$(foreach spc, $(AQTLVER), animalqtldb/$(spc))

animalqtldb/: ; mkdir $@

# AQTL_TMP_FILES
$(foreach spc, $(AQTLTMP), animalqtldb/$(spc)):
	$(CDAQTL) $(WGET) $(AQTLDL)/tmp/$(notdir $@)

# AQTL_VER_FILES
$(foreach spc, $(AQTLVER), animalqtldb/$(spc)):
	$(CDAQTL) $(WGET) $(AQTLDL)/export/KSUI8GFHOT6/$(notdir $@)

# GENEINFO_FILES
# these are all created under ncbigene first then linked here
# so the distinction of the locally generated ones becomes moot

$(foreach spc, $(AQTLGI), animalqtldb/$(spc)): $(foreach spc, $(AQTLGI),ncbigene/$(spc))
	unlink $@; $(CDAQTL) ln -s ../ncbigene/$(notdir $@) $(notdir $@)

animalqtldb_clean: ;  $(RM) animalqtldb/*
##########################################
CDBGE = cd bgee/ ;
bgee: dipper bgee/ \
		bgee/bgee.sqlite3.gz

bgee/: ; mkdir $@
bgee/bgee.sqlite3.gz: bgee/bgee_sqlite3.sql
	$(CDBG) gzip --force bgee.sqlite3  # appends .gz,  caused reindexing

bgee/bgee.sqlite3: bgee/bgee_sqlite3.sql
	$(CDBGE) /usr/bin/sqlite3 -mmap 3G bgee.sqlite3 < bgee_sqlite3.sql

bgee/bgee_sqlite3.sql:  bgee/sql_lite_dump.sql
	$(CDBGE) ../dipper/scripts/mysql2sqlite $? > $@ ;\
	echo -e "\nvacuum;analyze;" >> $@

bgee/sql_lite_dump.sql:  bgee/sql_lite_dump.tar.gz
	$(CDBGE) /bin/tar -xzf sql_lite_dump.tar.gz $(notdir $@)

bgee/sql_lite_dump.tar.gz:
	$(CDBGE) $(WGET) ftp://ftp.bgee.org/current/sql_lite_dump.tar.gz

bgee_clean: ; $(RM) bgee/*
########################################
BGDL = https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release
CDBOG = cd biogrid ;
biogrid: biogrid/ \
	biogrid/BIOGRID-ALL-LATEST.mitab.zip \
	biogrid/BIOGRID-IDENTIFIERS-LATEST.tab.zip

biogrid/: ; mkdir $@

biogrid/BIOGRID-ALL-LATEST.mitab.zip:
	$(CDBOG) $(WGET) $(BGDL)/BIOGRID-ALL-LATEST.mitab.zip
biogrid/BIOGRID-IDENTIFIERS-LATEST.tab.zip:
	$(CDBOG) $(WGET) $(BGDL)/BIOGRID-IDENTIFIERS-LATEST.tab.zip

# TODO get rid of obsuscating  name changes in py (then delete here too)
biogrid/identifiers.tab.zip:  biogrid/BIOGRID-IDENTIFIERS-LATEST.tab.zip
	unlink $@; $(CDBOG) \
	ln -s BIOGRID-IDENTIFIERS-LATEST.tab.zip identifiers.tab.zip
biogrid/interactions.mitab.zip: biogrid/BIOGRID-ALL-LATEST.mitab.zip
	unlink $@; $(CDBOG) \
	ln -s BIOGRID-ALL-LATEST.mitab.zip interactions.mitab.zip

biogrid_clean: ;  $(RM) biogrid/*
##########################################
CVFTP = ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar
clinvar: clinvar/ \
	clinvar/ClinVarFullRelease_00-latest.xml.gz \
	clinvar/gene_condition_source_id

clinvar/: ; mkdir $@

clinvar/ClinVarFullRelease_00-latest.xml.gz:
	cd clinvar; $(WGET) $(CVFTP)/xml/ClinVarFullRelease_00-latest.xml.gz
clinvar/gene_condition_source_id:
	cd clinvar; $(WGET) $(CVFTP)/gene_condition_source_id

clinvar_clean: ; $(RM) clinvar/*
##########################################
# public licence issues
#coriell: coriell/ ; mkdir $@
##########################################
ctd: ctd/ ctd/CTD_chemicals_diseases.tsv.gz

# (depreciated) CTD_genes_pathways.tsv.gz CTD_genes_diseases.tsv.gz
# there is an API not called here
# http://ctdbase.org/tools/batchQuery.go?q

ctd/: ; mkdir $@

ctd/CTD_chemicals_diseases.tsv.gz:
	cd ctd; $(WGET) http://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz

ctd_clean: ; $(RM) ctd/*
##########################################
# included here for dipper/scripts/
dipper: dipper/ dipper/scripts
dipper/: ; git clone https://github.com/monarch-initiative/dipper.git

dipper/scripts:
	cd dipper; git pull

dipper_clean: ; $(RM) dipper
##########################################
# used in more than one ingest
OBO = http://purl.obolibrary.org/obo
eco: eco/ \
	eco/gaf-eco-mapping.txt

eco/: ; mkdir $@

eco/gaf-eco-mapping.txt:
	cd eco/; $(WGET) $(OBO)/eco/gaf-eco-mapping.txt

eco_clean: ; $(RM) eco/*
##########################################
# Lots to reconsider here ...
#ENSURL = https://uswest.ensembl.org
# alternativly slower www.ensembl.org

# they provide per species turtle ensembl to external id mappings
# ftp.ensembl.org/pub/current_rdf/  species
#------------------------------------------------------------------------
# (cheap) favorite species in dipper can be found in the QC Readmes on M4
# (cheap) species available in ensembl
# curl ftp://ftp.ensembl.org/pub/current_rdf/ > ensembl_species_name
# for ds in $(cat dipper_species_name); do grep -E " $ds$" ensembl_species_name ; done | cut -c57- > names_to_fetch

# ordered by "popularity"
ENSSPC = mus_musculus homo_sapiens drosophila_melanogaster bos_taurus danio_rerio
#	caenorhabditis_elegans sus_scrofa rattus_norvegicus gallus_gallus \
#	canis_lupus_familiaris pan_troglodytes macaca_mulatta monodelphis_domestica \
#	equus_caballus felis_catus ornithorhynchus_anatinus	arabidopsis_thaliana \
#	xenopus_(silurana)_tropicalis anolis_carolinensis takifugu_rubripe \
#	saccharomyces_cerevisiae_s288c schizosaccharomyces_pombe dictyostelium_discoideum \
#	ovis_aries escherichia_coli oncorhynchus_mykiss mus_setulosus mus \
#	mus_musculus_molossinus	mus_spretus	mus_musculus_musculus mus_musculus_castaneus \
#	sus_scrofa_domestica mus_musculus_musculus_x_m._m._domesticus capra_hircus \
#	macaca_nemestrina mus_musculus_bactrianus
ENSRDF := ftp://ftp.ensembl.org/pub/current_rdf
ENSRDF_TARGET := $(foreach species, $(ENSSPC), ensrdf/$(species).ttl.gz)
ensrdf:  ensrdf/ \
		$(ENSRDF_TARGET)

ensrdf/: ; mkdir $@
$(ENSRDF_TARGET):
	cd ensrdf; $(WGET) $(ENSRDF)/$(subst ensrdf/,,$(subst .ttl.gz,,$@))/$(subst ensrdf/,,$@)


# TODO Consider transplanting dipper python api fetch code to here?
#ENSMRT = $(ENSURL)/biomart/martservice?query=
#ensembl: ensembl/ \
#	ensembl/ensembl_9606.txt \
#	ensembl/ensembl_7955.txt \
#	ensembl/ensembl_10090.txt \
#	ensembl/ensembl_28377.txt \
#	ensembl/ensembl_3702.txt \
#	ensembl/ensembl_9913.txt \
#	ensembl/ensembl_6239.txt \
#	ensembl/ensembl_9615.txt \
#	ensembl/ensembl_9031.txt \
#	ensembl/ensembl_44689.txt \
#	ensembl/ensembl_7227.txt \
#	ensembl/ensembl_9796.txt \
#	ensembl/ensembl_9544.txt \
#	ensembl/ensembl_13616.txt \
#	ensembl/ensembl_9258.txt \
#	ensembl/ensembl_9823.txt \
#	ensembl/ensembl_10116.txt \
#	ensembl/ensembl_4896.txt \
#	ensembl/ensembl_31033.txt \
#	ensembl/ensembl_8364.txt \
#	ensembl/ensembl_4932.txt
#ensembl/: ; mkdir $@
#ensembl/ensembl_9606.txt:

#ensembl_clean: ;  $(RM) ensembl
##########################################
#eom: eom/ ; mkdir $@
# a dead end fron the 'disco' days?
##########################################
FLYFTP = ftp://ftp.flybase.net
FLYPRE = releases/current/precomputed_files
CDFLY = cd flybase/ ;

flybase: flybase/ \
		 flybase/md5sum.txt \
		 flybase/disease_model_annotations.tsv.gz \
		 flybase/species.ab.gz \
		 flybase/fbal_to_fbgn_fb.tsv.gz \
		 flybase/fbrf_pmid_pmcid_doi_fb.tsv.gz

flybase/: ; mkdir $@

flybase/md5sum.txt:
	$(CDFLY) $(WGET) $(FLYFTP)/$(FLYPRE)/md5sum.txt

flybase/disease_model_annotations.tsv.gz: flybase/md5sum.txt
	$(CDFLY)  \
	fname=$$(fgrep "/human_disease/disease_model_annotations" md5sum.txt| cut -f2- -d'/') ; \
	$(WGET) $(FULLPTH) $(FLYFTP)/$(FLYPRE)/$$fname ; \
	unlink disease_model_annotations.tsv.gz; \
	ln -s $(FLYPRE)/$$fname  disease_model_annotations.tsv.gz

flybase/species.ab.gz: flybase/md5sum.txt
	$(CDFLY) $(WGET) $(FLYFTP)/$(FLYPRE)/species/species.ab.gz

flybase/fbal_to_fbgn_fb.tsv.gz: flybase/md5sum.txt
	$(CDFLY)  \
	fname=$$(fgrep "alleles/fbal_to_fbgn_fb" md5sum.txt| cut -f2- -d'/') ; \
	$(WGET) $(FULLPTH) $(FLYFTP)/$(FLYPRE)/$$fname ; \
	unlink$(notdir $@) ; ln -s $(FLYPRE)/$$fname $(notdir $@)

flybase/fbrf_pmid_pmcid_doi_fb.tsv.gz: flybase/md5sum.txt
	$(CDFLY) \
	fname=$$(fgrep "references/fbrf_pmid_pmcid_doi_fb" md5sum.txt| cut -f2- -d'/') ; \
	$(WGET) $(FULLPTH) $(FLYFTP)/$(FLYPRE)/$$fname ; \
	unlink $(notdir $@); ln -s $(FLYPRE)/$$fname $(notdir $@)

flybase_clean: ;  $(RM) flybase/*

##########################################
GRDL = http://ftp.ncbi.nih.gov/pub/GeneReviews
genereviews: genereviews/ \
		genereviews/NBKid_shortname_OMIM.txt \
		genereviews/GRtitle_shortname_NBKid.txt

genereviews/: ; mkdir $@
genereviews/NBKid_shortname_OMIM.txt:
	cd genereviews; $(WGET) $(GRDL)/NBKid_shortname_OMIM.txt
genereviews/GRtitle_shortname_NBKid.txt:
	cd genereviews; $(WGET) $(GRDL)/GRtitle_shortname_NBKid.txt

genereviews_clean: ;  $(RM) genereviews/*

##########################################
GOADL := http://current.geneontology.org
FTPEBI := ftp://ftp.uniprot.org/pub/databases
UPCRKB := uniprot/current_release/knowledgebase

GOASPC = fb \
		zfin \
		mgi \
		rgd \
		wb \
		goa_pig \
		goa_chicken \
		goa_human \
		goa_cow \
		goa_dog \
		sgd \
		pombase \
		dictybase \
		aspgd

go: eco go/ \
	$(foreach species, $(GOASPC), go/$(species).gaf.gz) \
	go/go-refs.json \
	go/GO.references \
	go/idmapping_selected.tab.gz \
	go/gaf-eco-mapping.txt

go/: ; mkdir $@

$(foreach species, $(GOASPC), go/$(species).gaf.gz):
	cd go/; $(WGET) $(GOADL)/annotations/$(notdir $@)

go/go-refs.json:
	cd go/; $(WGET) http://current.geneontology.org/metadata/go-refs.json
go/GO.references:  # TODO depreicated in favor of go-refs.json
	cd go/; $(WGET) http://www.geneontology.org/doc/GO.references

go/idmapping_selected.tab.gz:  # expensive
	cd go/; $(WGET) $(FTPEBI)/$(UPCRKB)/idmapping/idmapping_selected.tab.gz
go/gaf-eco-mapping.txt: eco/gaf-eco-mapping.txt
	unlink $@;cd go/; ln -s ../$< gaf-eco-mapping.txt

go_clean: ; $(RM) go/*
##########################################
GWASFTP = ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest
GWASFILE = gwas-catalog-associations_ontology-annotated.tsv
gwascatalog: gwascatalog/ gwascatalog/$(GWASFILE)

gwascatalog/: ; mkdir $@
gwascatalog/$(GWASFILE):
	cd gwascatalog; $(WGET) $(GWASFTP)/$(GWASFILE)
# SO & MONDO ontologies need their own cache
gwascatalog_clean:  ;  $(RM) gwascatalog/*
##########################################
EBIFTP := ftp://ftp.ebi.ac.uk
EBIPTH := pub/databases/genenames/new/tsv

hgnc: hgnc/ \
		hgnc/hgnc_complete_set.txt

hgnc/: ; mkdir $@
hgnc/hgnc_complete_set.txt:
	cd hgnc/; $(WGET) $(EBIFTP)/$(EBIPTH)/hgnc_complete_set.txt

# Warning: File 'hgnc/hgnc_complete_set.txt' has modification time 18212 s in the future

hgnc_clean:  ;  $(RM) hgnc/*
##########################################
# fragile
PNR = http://compbio.charite.de/jenkins/job/hpo.annotations.current
HPOADL2 = $(PNR)/lastSuccessfulBuild/artifact/misc_2018
hpoa: hpoa/ \
	hpoa/phenotype.hpoa

hpoa/: ; mkdir $@
hpoa/phenotype.hpoa:
	cd hgnc; $(WGET) $(HPOADL2)/phenotype.hpoa

hpoa_clean:  ; $(RM) hpoa/*
##########################################
IMPCDL = ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv
impc: impc/  impc/checksum.md5 \
	  impc/ALL_genotype_phenotype.csv.gz

impc/: ; mkdir $@
impc/checksum.md5:
	cd impc; $(WGET) $(IMPCDL)/checksum.md5
impc/ALL_genotype_phenotype.csv.gz: impc/checksum.md5
	cd impc; $(WGET) $(IMPCDL)/ALL_genotype_phenotype.csv.gz

impc_clean:  ;  $(RM) impc/*
##########################################
#/list
KEGGG = http://rest.genome.jp
# link
KEGGK = http://rest.kegg.jp

kegg: kegg/ \
	kegg/disease \
	kegg/pathway \
	kegg/hsa_genes \
	kegg/orthology \
	kegg/disease_gene \
	kegg/omim \
	kegg/omim2gene \
	kegg/human_gene2pathway \
	kegg/hsa_orthologs \
    kegg/mmu \
    kegg/rno \
    kegg/dme \
    kegg/dre \
    kegg/cel \
    kegg/pubmed \
    kegg/ds \
    kegg/ko

kegg/: 	; mkdir $@
# note choosing native name except when there is a conflict
# conflicts interfer with --timestamping
kegg/disease:
	cd kegg; $(WGET) $(KEGGG)/list/disease
kegg/pathway:
	cd kegg; $(WGET) $(KEGGG)/list/pathway
kegg/hsa_genes:
	cd kegg; $(WGET) $(KEGGG)/list/hsa -O hsa_genes
kegg/orthology:
	cd kegg; $(WGET) $(KEGGG)/list/orthology
kegg/disease_gene:
	cd kegg; $(WGET) $(KEGGK)/link/disease/hsa -O disease_gene
kegg/omim:
	cd kegg; $(WGET) $(KEGGG)/link/disease/omim
kegg/omim2gene:
	cd kegg; $(WGET) $(KEGGG)/link/omim/hsa -O omim2gene
kegg/ncbi:
	cd kegg; $(WGET) $(KEGGG)/conv/ncbi-geneid/hsa -O ncbi
kegg/human_gene2pathway:
	cd kegg; $(WGET) $(KEGGK)/link/pathway/hsa -O human_gene2pathway
kegg/hsa_orthologs:
	cd kegg; $(WGET) $(KEGGK)/link/orthology/hsa -O hsa_orthologs
kegg/mmu:
	cd kegg; $(WGET) $(KEGGK)/link/orthology/mmu
kegg/rno:
	cd kegg; $(WGET) $(KEGGK)/link/orthology/rno
kegg/dme:
	cd kegg; $(WGET) $(KEGGK)/link/orthology/dme
kegg/dre:
	cd kegg; $(WGET) $(KEGGK)/link/orthology/dre
kegg/cel:
	cd kegg; $(WGET) $(KEGGK)/link/orthology/cel
kegg/pubmed:
	cd kegg; $(WGET) $(KEGGK)/link/pathway/pubmed
kegg/ds:
	cd kegg; $(WGET) $(KEGGK)/link/pathway/ds
kegg/ko:
	cd kegg; $(WGET) $(KEGGK)/link/pathway/ko

kegg_clean:  ;  $(RM) kegg/*
##########################################
# pulls via sql queries
#mgi: mgi/
#mgi/; mkdir $@
##########################################
mmrrc: mmrrc/ mmrrc/mmrrc_catalog_data.csv
mmrrc/: ; mkdir $@
mmrrc/mmrrc_catalog_data.csv:
	cd mmrrc; $(WGET) https://www.mmrrc.org/about/mmrrc_catalog_data.csv

mmrrc_clean:  ;  $(RM) mmrrc/*
##########################################
# private -- why?
#monarch: monarch/
#monarch/: ; mkdir $@
#monarch_clean:  ;  $(RM) $(substr _clean,/,$@)
##########################################
# TODO reconsider layout make ucscbands prime
MCDL = http://hgdownload.cse.ucsc.edu/goldenPath
CDMC = cd monochrom/ ;
CBI = cytoBandIdeo.txt.gz
monochrom: monochrom/ \
	monochrom/9606cytoBand.txt.gz \
	monochrom/10090cytoBand.txt.gz \
	monochrom/7955cytoBand.txt.gz \
	monochrom/10116cytoBand.txt.gz \
	monochrom/bosTau7cytoBand.txt.gz \
	monochrom/galGal4cytoBand.txt.gz \
	monochrom/susScr3cytoBand.txt.gz \
	monochrom/oviAri3cytoBand.txt.gz \
	monochrom/equCab2cytoBand.txt.gz

monochrom/: ; mkdir $@
# accomadate existing ingest given names
monochrom/9606cytoBand.txt.gz: monochrom/hg19/ monochrom/hg19/cytoBand.txt.gz
	$(CDMC) unlink 9606cytoBand.txt.gz ;\
	ln -s hg19/cytoBand.txt.gz 9606cytoBand.txt.gz
monochrom/hg19/: ; mkdir $@
monochrom/hg19/cytoBand.txt.gz:
	cd monochrom/hg19; $(WGET) $(MCDL)/hg19/database/cytoBand.txt.gz

monochrom/10090cytoBand.txt.gz: monochrom/mm10/ monochrom/mm10/$(CBI)
	$(CDMC) unlink 10090cytoBand.txt.gz ;\
	ln -s mm10/$(CBI) 10090cytoBand.txt.gz  # note dropping 'Ideo' (to fix)
monochrom/mm10/: ; mkdir $@
monochrom/mm10/$(CBI):
	cd monochrom/mm10/ ; $(WGET) $(MCDL)/mm10/database/$(CBI)

monochrom/7955cytoBand.txt.gz: monochrom/danRer10/ monochrom/danRer10/$(CBI)
	$(CDMC) unlink 7955cytoBand.txt.gz ; \
	ln -s  danRer10/$(CBI) 7955cytoBand.txt.gz
monochrom/danRer10/: ; mkdir $@
monochrom/danRer10/$(CBI):
	cd monochrom/danRer10/ ;  $(WGET) $(MCDL)/danRer10/database/$(CBI)

monochrom/10116cytoBand.txt.gz: monochrom/rn6/  monochrom/rn6/$(CBI)
	$(CDMC) unlink 10116cytoBand.txt.gz ; \
	ln -s rn6/$(CBI) 10116cytoBand.txt.gz
monochrom/rn6/: ; mkdir $@
monochrom/rn6/$(CBI):
	cd monochrom/rn6/; $(WGET) $(MCDL)/rn6/database/$(CBI)

monochrom/bosTau7cytoBand.txt.gz: monochrom/bosTau7/ monochrom/bosTau7/$(CBI)
	$(CDMC) unlink bosTau7cytoBand.txt.gz; \
	ln -s bosTau7/$(CBI) bosTau7cytoBand.txt.gz
monochrom/bosTau7/: ; mkdir $@

monochrom/bosTau7/$(CBI):
	cd monochrom/bosTau7/; $(WGET) $(MCDL)/bosTau7/database/$(CBI)

monochrom/galGal4cytoBand.txt.gz: monochrom/galGal4/ monochrom/galGal4/$(CBI)

	$(CDMC) unlink galGal4cytoBand.txt.gz; \
	ln -s galGal4/$(CBI) galGal4cytoBand.txt.gz
monochrom/galGal4/: ; mkdir $@
monochrom/galGal4/cytoBandIdeo.txt.gz:
	cd monochrom/galGal4/; $(WGET) $(MCDL)/galGal4/database/$(CBI)

monochrom/susScr3cytoBand.txt.gz: monochrom/susScr3/ monochrom/susScr3/$(CBI)
	$(CDMC) unlink susScr3cytoBand.txt.gz ; \
	ln -s susScr3/$(CBI) susScr3cytoBand.txt.gz
monochrom/susScr3/: ; mkdir $@
monochrom/susScr3/cytoBandIdeo.txt.gz:
	cd monochrom/susScr3/; $(WGET) $(MCDL)/susScr3/database/$(CBI)

monochrom/oviAri3cytoBand.txt.gz: monochrom/oviAri3/ monochrom/oviAri3/$(CBI)
	$(CDMC) unlink oviAri3cytoBand.txt.gz ; \
	ln -s oviAri3/$(CBI) oviAri3cytoBand.txt.gz
monochrom/oviAri3/: ; mkdir $@
monochrom/oviAri3/cytoBandIdeo.txt.gz:
	cd monochrom/oviAri3/; $(WGET) $(MCDL)/oviAri3/database/$(CBI)

monochrom/equCab2cytoBand.txt.gz: monochrom/equCab2/ monochrom/equCab2/$(CBI)
	$(CDMC) unlink equCab2cytoBand.txt.gz ; \
	ln -s equCab2$(CBI) equCab2cytoBand.txt.gz
monochrom/equCab2/: ; mkdir $@
monochrom/equCab2/cytoBandIdeo.txt.gz:
	cd monochrom/equCab2/; $(WGET) $(MCDL)/equCab2/database/$(CBI)

monochrom_clean: ;  $(RM) monochrom/*
##########################################
MPDDL = https://phenomedoc.jax.org/MPD_downloads
mpd: mpd/ \
	mpd/ontology_mappings.csv \
	mpd/straininfo.csv \
	mpd/measurements.csv \
# 	mpd_datasets_metadata.xml.gz  # TODO should this remain unused?

mpd/: ; mkdir $@
mpd/ontology_mappings.csv:
	cd mpd; $(WGET) $(MPDDL)/ontology_mappings.csv
mpd/straininfo.csv:
	cd mpd; $(WGET) $(MPDDL)/straininfo.csv
mpd/measurements.csv:
	cd mpd; $(WGET) $(MPDDL)/measurements.csv
mpd/strainmeans.csv.gz:
	cd mpd; $(WGET) $(MPDDL)/strainmeans.csv.gz

mpd_clean: ;  $(RM) mpd/*
##########################################
#mychem: mychem/ ; mkdir $@
# api driven
##########################################
NCBIFTP = ftp://ftp.ncbi.nih.gov/gene/DATA
ncbigene: ncbigene/ \
	ncbigene/gene_info.gz \
	ncbigene/gene_history.gz \
	ncbigene/gene2pubmed.gz \
	ncbigene/gene_group.gz \
	ncbigene/Gallus_gallus.gene_info.gz \
	ncbigene/Sus_scrofa.gene_info.gz \
	ncbigene/Bos_taurus.gene_info.gz \
	ncbigene/Equus_caballus.gene_info.gz \
	ncbigene/Ovis_aries.gene_info.gz \
	ncbigene/Oncorhynchus_mykiss.gene_info.gz

ncbigene/: ; mkdir $@

ncbigene/gene_info.gz:
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene_info.gz
	@ # test if in the future (GMT?)
	@ # touch -r "gene_info.gz" -d '-8 hour' "gene_info.gz"
ncbigene/gene_history.gz:
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene_history.gz
ncbigene/gene2pubmed.gz:
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene2pubmed.gz
ncbigene/gene_group.gz:
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene_group.gz

GENEINFO = $(NCBIFTP)/GENE_INFO

ncbigene/Gallus_gallus.gene_info.gz:
	cd ncbigene/ ; $(WGET) $(GENEINFO)/Non-mammalian_vertebrates/Gallus_gallus.gene_info.gz
ncbigene/Sus_scrofa.gene_info.gz:
	cd ncbigene/ ; $(WGET) $(GENEINFO)/Mammalia/Sus_scrofa.gene_info.gz
ncbigene/Bos_taurus.gene_info.gz:
	cd ncbigene/ ; $(WGET) $(GENEINFO)/Mammalia/Bos_taurus.gene_info.gz

# GENEINFO_LOCAL_FILES	(maybe partion out more if they could help with other ingests)
ncbigene/Equus_caballus.gene_info.gz: ncbigene/gene_info.gz
	cd ncbigene/ ; /bin/zgrep "^9796[^0-9]" gene_info.gz | /bin/gzip > Equus_caballus.gene_info.gz
ncbigene/Ovis_aries.gene_info.gz: ncbigene/gene_info.gz
	cd ncbigene/ ; /bin/zgrep "^9940[^0-9]" gene_info.gz | /bin/gzip > Ovis_aries.gene_info.gz
ncbigene/Oncorhynchus_mykiss.gene_info.gz: ncbigene/gene_info.gz
	cd ncbigene/ ; /bin/zgrep "^8022[^0-9]" gene_info.gz | /bin/gzip > Oncorhynchus_mykiss.gene_info.gz

ncbigene_clean: ;  $(RM) ncbigene/*
##########################################
#
omia: omia/ \
	omia/omia.xml.gz
#	causal_mutations.tab   # TODO not in use

omia/: ; mkdir $@
omia/omia.xml.gz:
	cd omia; $(WGET) http://compldb.angis.org.au/dumps/omia.xml.gz
	# http://omia.angis.org.au/dumps/omia.xml.gz  # broken alt

omia_clean: ; $(RM) omia/*
##########################################
#omim: omim/ ; mkdir $@
# private api & keys can't go here
##########################################
orphanet: orphanet/ \
		orphanet/en_product6.xml

orphanet/: ; mkdir $@
orphanet/en_product6.xml:
	cd orphanet/; $(WGET) http://www.orphadata.org/data/xml/en_product6.xml

orphanet_clean: ; $(RM) orphanet/*
##########################################
PNTHDL = ftp://ftp.pantherdb.org/ortholog/current_release
panther: panther/ \
	panther/RefGenomeOrthologs.tar.gz \
	panther/Orthologs_HCOP.tar.gz

panther/: ; mkdir $@
panther/RefGenomeOrthologs.tar.gz:
	cd panther; $(WGET) $(PNTHDL)/RefGenomeOrthologs.tar.gz
panther/Orthologs_HCOP.tar.gz:
	cd panther; $(WGET) $(PNTHDL)/Orthologs_HCOP.tar.gz

panther_clean: ; $(RM) panther/*
##########################################
RCTDL = https://www.reactome.org/download/current
reactome: reactome/ \
		reactome/Ensembl2Reactome.txt \
		reactome/ChEBI2Reactome.txt \
		reactome/gaf-eco-mapping.txt

reactome/: ; mkdir $@
reactome/Ensembl2Reactome.txt:
	cd reactome; $(WGET) $(RCTDL)/Ensembl2Reactome.txt
reactome/ChEBI2Reactome.txt:
	cd reactome; $(WGET) $(RCTDL)/ChEBI2Reactome.txt

reactome/gaf-eco-mapping.txt:  eco/gaf-eco-mapping.txt
	cd reactome; unlink $(notdir $@); ln -s ../$< $(notdir $@)

reactome_clean: ; $(RM) rectome/*
##########################################

RGDFTP = ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology
rgd: rgd/ \
	rgd/rattus_genes_mp

rgd/: ; mkdir $@
rgd/rattus_genes_mp:
	cd rgd/; $(WGET) $(RGDFTP)/rattus_genes_mp

rgd_clean: ; $(RM) rdg/*
##########################################
SGDDL = https://downloads.yeastgenome.org/curation/literature
sgd: sdg/ sdg/phenotype_data.tab

sdg/: ; mkdir $@
sdg/phenotype_data.tab:
	cd sdg; $(WGET) $(SGDDL)/phenotype_data.tab
sgd_clean: ; $(RM) sdg/*
##########################################
STRING = https://string-db.org
# redirects on dl to a static file server
STRSTA = https://stringdb-static.org

STRDL = $(STRSTA)/download
STRMAP = $(STRING)/mapping_files/entrez

# not hard coding these would be better ...
STRVER = 11.0
STRYR = 2018

SRTPTH = protein.links.detailed.v$(STRVER)

STRTAX = 9606 10090 7955 7227 6239 4932
STRSPC = celegans fly human mouse yeast zebrafish

string: string/ \
		string/version \
		$(foreach txid, $(STRTAX), string/$(txid).$(SRTPTH).txt.gz) \
		$(foreach species, $(STRSPC), string/$(species).entrez_2_string.$(STRYR).tsv.gz)

string/: ; mkdir $@

string/version:
	cd string ; $(WGET) $(STRING)/api/tsv-no-header/version  ; \
	if [ "$(STRVER)" != "$$(cut -f 1 version)" ] ;then echo "NEW VERSION of STRING!" ; \
	else  echo "same version of STRING" ; fi

$(foreach txid, $(STRTAX), string/$(txid).$(SRTPTH).txt.gz):
	cd string; $(WGET) $(STRDL)/$(SRTPTH)/$(notdir $@)
	#cd string; $(WGET) $(STRDL)/$(SRTPTH)/$(subst string/,,$@)

$(foreach species, $(STRSPC), string/$(species).entrez_2_string.$(STRYR).tsv.gz):
	cd string; $(WGET) $(STRMAP)/$(notdir $@)

string_clean: ;  $(RM) string/*
##########################################
# this data is shared  with Monochrom.
# need to choose one (this one) and make the other symlinks
#
#HGGP = http://hgdownload.cse.ucsc.edu/goldenPath
#ucscbands:  ucscbands/ \
#	hg19_cytoBand.txt.gz \
#	mm10_cytoBandIdeo.txt.gz \
#	danRer11_cytoBandIdeo.txt.gz \
#
#
#
#ucscbands/: ; mkdir $@
#
#HGGP + '/hg19/database/cytoBand.txt.gz
#HGGP + '/mm10/database/cytoBandIdeo.txt.gz
#HGGP + '/danRer11/database/cytoBandIdeo.txt.gz
#
#
#ucscbands_clean: ; $(RM) ucscbands/*
##########################################
#udp: udp/ ; mkdir $@
#  private access
##########################################
WBFTP = ftp://ftp.wormbase.org
#WBDEV = /pub/wormbase/releases//current-development-release  unused?
WBPROD = pub/wormbase/releases/current-production-release
WBSPC = c_elegans/PRJNA13758
CDWB = cd wormbase/ ;
WSNUM = sed -n '1 s|.*\.\(WS[0-9]\{3\}\)\..*|\1|p' CHECKSUMS
wormbase: wormbase/ \
		wormbase/CHECKSUMS \
		wormbase/letter \
		wormbase/c_elegans.PRJNA13758.geneIDs.txt.gz  \
		wormbase/c_elegans.PRJNA13758.annotations.gff3.gz \
		wormbase/c_elegans.PRJNA13758.xrefs.txt.gz \
		wormbase/phenotype_association.wb \
		wormbase/rnai_phenotypes.wb \
		wormbase/disease_association.wb \
		wormbase/pub_xrefs.txt

wormbase/: ; mkdir $@
wormbase/CHECKSUMS:
	$(CDWB) $(WGET) $(WBFTP)/$(WBPROD)/$(notdir $@)
wormbase/letter:  wormbase/CHECKSUMS
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)); \
	$(WGET) $(WBFTP)/$(WBPROD)/letter_$$wbnum ;\
	ln -s /letter_$$wbnum $(notdir $@)
wormbase/c_elegans.PRJNA13758.geneIDs.txt.gz:  wormbase/CHECKSUMS
	#species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS273.geneIDs.txt.gz
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/species/$(WBSPC)/annotation/$(WBSPC).$$wbnum.geneIDs.txt.gz;\
	ln -s $(WBPROD)/species/$(WBSPC)/annotation/$(WBSPC).$$wbnum.geneIDs.txt.gz $(notdir $@)
wormbase/c_elegans.PRJNA13758.annotations.gff3.gz: wormbase/CHECKSUMS
	# species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS273.annotations.gff3.gz
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/species/$(WBSPC)/$(WBSPC).$$wbnum.annotations.gff3.gz;\
	ln -s $(WBPROD)/species/$(WBSPC)/$(WBSPC).$$wbnum.annotations.gff3.gz $(notdir $@)
wormbase/c_elegans.PRJNA13758.xrefs.txt.gz: wormbase/CHECKSUMS
	# species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS273.xrefs.txt.gz
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/species/$(WBSPC)/annotation/$(WBSPC).$$wsnum.xrefs.txt.gz;\
	ln -s $(WBPROD)/species/$(WBSPC)/annotation/$(WBSPC).$$wsnum.xrefs.txt.gz $(notdir $@)
wormbase/phenotype_association.wb: wormbase/CHECKSUMS
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/ONTOLOGY/phenotype_association.$$wsnum.wb;\
	ln -s $(WBPROD)/ONTOLOGY/phenotype_association.$$wsnum.wb $(notdir $@)
wormbase/rnai_phenotypes.wb: wormbase/CHECKSUMS
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/ONTOLOGY/rnai_phenotypes.$$wsnum.wb;\
	ln -s $(WBPROD)/ONTOLOGY/rnai_phenotypes.$$wsnum.wb $(notdir $@)
wormbase/disease_association.wb: wormbase/CHECKSUMS
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET $(FULLPTH) $(WBFTP)/$(WBPROD)/ONTOLOGY/rnai_phenotypes.$$wsnum.wb;\
	ln -s $(WBPROD)/ONTOLOGY/rnai_phenotypes.$$wsnum.wb $(notdir $@)
# api call so no date or file version
wormbase/pub_xrefs.txt:
	$(CDWB) $(WGET) -O $(notdir $@) \
	http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref

wormbase_clean: ;  $(RM) wormbase/*
##########################################
ZFDL = https://zfin.org/downloads
# not        https://purl.obolibrary.org/obo/zp/src/curation
ZPCURATION = http://purl.obolibrary.org/obo/zp
ZFFILES = \
	genotype_features.txt \
	phenotype_fish.txt \
	zfinpubs.txt \
	Morpholinos.txt \
	pheno_environment_fish.txt \
	stage_ontology.txt \
	mappings.txt \
	genotype_backgrounds.txt \
	genbank.txt \
	uniprot.txt \
	gene.txt \
	human_orthos.txt \
	features.txt \
	features-affected-genes.txt \
	gene_marker_relationship.txt \
	CRISPR.txt \
	TALEN.txt \
	pub_to_pubmed_id_translation.txt \
	E_zfin_gene_alias.gff3 \
	fish_model_disease.txt \
	fish_components_fish.txt \
	phenoGeneCleanData_fish.txt

zfin: zfin/ \
		$(addprefix zfin/, $(ZFFILES)) \
		zfin/id_map_zfin.tsv

zfin/: ; mkdir $@

$(addprefix zfin/, $(ZFFILES)):
	cd zfin/; $(WGET) $(ZFDL)/$(notdir $@)

zfin/id_map_zfin.tsv:
	cd zfin/; $(WGET) $(ZPCURATION)/$(notdir $@)

zfin_clean: ;  $(RM) zfin/*
##########################################
zfinslim:  zfin \
		zfin/phenoGeneCleanData_fish.txt \
		zfin/id_map_zfin.tsv

# zfinslim/: ; mkdir $@  # just stick it in zdin

zfinslim/phenoGeneCleanData_fish.txt: zfin/phenoGeneCleanData_fish.txt
zfinslim/zp-mapping-2019.txt: zfin/id_map_zfin.tsv

zfinslim_clean: ; $(RM) zfin/id_map_zfin.tsv zfin/phenoGeneCleanData_fish.txt
##########################################


clean: animalqtldb_clean bgee_clean clinvar_clean  ctd_clean eco_clean\
	flybase_clean genereviews_clean go_clean gwascatalog_clean hgnc_clean \
	hpoa_clean impc_clean kegg_clean  mmrrc_clean  monochrom_clean \
	mpd_clean ncbigene_clean omia_clean  orphanet_clean panther_clean reactome_clean \
	rgd_clean sgd_clean string_clean  wormbase_clean zfin_clean zfinslim_clean
	# ensembl_clean  eom_clean # mychem_clean # mychem_clean # omim_clean # udp_clean
	# mgi_clean monarch_clean ucscbands_clean
