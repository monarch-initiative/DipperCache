
# Since this is not performing any traditional C build tasks avoid
# unexpected behavior with accidental collisions with preset variables/rules
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

WGET = /usr/bin/wget --timestamping --no-verbose
# addtional wget arguments
FULLPTH := --force-directories --no-host-directories

RM := rm --force --recursive --verbose

CLEAN = rm --force --recursive --verbose $(dir $@)/*

# For the first dependency to appear as the target
# note $(realpath) does not seem to be available for
# || ($(realpath $@) !=  $(realpath $<)) ] ; then ...
#SYMLINK = if [ ! -L "$@" ] ; then cd  $(dir $@); unlink "$(notdir $@)"; ln -sf "../$<" "$(notdir $@)"; fi
SYMLINK = ln --force --no-dereference --symbolic $< $@; fi

# when a remote server does not set last-modified headers we can only test
# if a new file is the same or different from the one we already heve.
# Last-modified header missing workaround
COPYCHANGED = if [ "$$(md5sum $<|cut -c 1-32)" != "$$(md5sum $@|cut -c 1-32)" ] ; then cp -fp $< $@ ; fi

# May be used in more than one ingest
OBO = http://purl.obolibrary.org/obo
GITRAW = https://raw.githubusercontent.com
GITMON = $(GITRAW)/monarch-initiative


.PHONY: help cruft recent clean dipper/scripts tree

SOURCES = animalqtldb \
	biogrid \
	bgee \
	clinvar \
	ctd \
	eco \
	flybase \
	genereviews \
	go \
	gwascatalog \
	hgnc \
	hpoa \
	impc \
	kegg  \
	mmrrc  \
	monochrom \
	mpd \
	ncbigene \
	omia  \
	orphanet \
	panther \
	reactome \
	rgd \
	sgd \
	string  \
	wormbase \
	zfin \
	zfinslim

all:  $(SOURCES) dipper

# unsatisifiable target used as a depencency
# to force recipies within independent targets
FORCE:

# ommited for cause
# coriell ensembl eom mgi monarch mychem mychem omim ucscbands

# Borrowed from https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
help: ## Display this help section
	@awk 'BEGIN {FS = ":.*?## "}\
		/^[a-zA-Z0-9_-]+:.*?## / {\
			printf "\033[36m%-38s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
#    .DEFAULT_GOAL := help

recent:  ## See what has been renewed in the last several days
	find ./*/ -mtime -5 -ls

cruft:  ## Preview which commands would be executed now
	make -n | sed 's|^mkdir |\nmkdir |g' > $@

tree:  ## Diffable metadata snapshot of results
	/usr/bin/tree --sort=name --timefmt "%Y%0m%0d %0T" -D -si -f -I dipper -n -o dipper_cache.tree ;\
	# git diff $@; git add $@ ;git commit -m "generated" $@  # not till in repo

# report:  ## A set of metrics to help track and access the health of the cache


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
		QTL_OAR_3.1.gff.txt.gz \
		QTL_SS_10.2.gff.txt.gz

AQTLVER = pig_QTLdata.txt \
		sheep_QTLdata.txt \
		cattle_QTLdata.txt \
		chicken_QTLdata.txt \
		horse_QTLdata.txt \
		rainbow_trout_QTLdata.txt

animalqtldb: ncbigene animalqtldb/ ## \
		$(foreach spc, $(AQTLGI), animalqtldb/$(spc)) \
		$(foreach spc, $(AQTLTMP), animalqtldb/$(spc)) \
		$(foreach spc, $(AQTLVER), animalqtldb/$(spc)) \
		animalqtldb/trait_mappings.csv

animalqtldb/: ; mkdir $@

animalqtldb/trait_mappings.csv: FORCE
	$(CDAQTL) $(WGET) $(AQTLDL)/export/$(notdir $@)

# AQTL_TMP_FILES
$(foreach spc, $(AQTLTMP), animalqtldb/$(spc)): FORCE
	$(CDAQTL) $(WGET) $(AQTLDL)/tmp/$(notdir $@)

# AQTL_VER_FILES
$(foreach spc, $(AQTLVER), animalqtldb/$(spc)): FORCE
	$(CDAQTL) $(WGET) $(AQTLDL)/export/KSUI8GFHOT6/$(notdir $@)

# GENEINFO_FILES
# these are all created under ncbigene first then linked here
# so the distinction of the locally generated ones becomes moot

$(foreach spc, $(AQTLGI), animalqtldb/$(spc)): $(foreach spc, $(AQTLGI),ncbigene/$(spc))
	# unlink $@; $(CDAQTL) ln -s ../ncbigene/$(notdir $@) $(notdir $@)
	$(SYMLINK)

animalqtldb_clean: ;  $(CLEAN)  # $(RM) animalqtldb/*
##########################################
CDBGE = cd bgee/ ;
bgee: dipper bgee/ \
		bgee/bgee.sqlite3.gz

bgee/: ; mkdir $@
bgee/bgee.sqlite3.gz: bgee/bgee_sqlite3.sql
	gzip --stdout $? > $@

bgee/bgee.sqlite3: bgee/bgee_sqlite3.sql
	$(CDBGE) /usr/bin/sqlite3 -mmap 3G bgee.sqlite3 < bgee_sqlite3.sql

bgee/bgee_sqlite3.sql:  bgee/sql_lite_dump.sql
	./dipper/scripts/mysql2sqlite $? > $@ ;\
	echo -e "\nvacuum;analyze;" >> $@

bgee/sql_lite_dump.sql:  bgee/sql_lite_dump.tar.gz
	$(CDBGE) /bin/tar -xzf sql_lite_dump.tar.gz $(notdir $@)

bgee/sql_lite_dump.tar.gz: FORCE
	$(CDBGE) $(WGET) ftp://ftp.bgee.org/current/sql_lite_dump.tar.gz

bgee_clean: ; $(RM) bgee/*
########################################
BGFP = Download/BioGRID/Latest-Release
BGDL = https://downloads.thebiogrid.org/$(BGFP)

CDBOG = cd biogrid ;
biogrid: biogrid/ \
	biogrid/$(BGFP)/BIOGRID-ALL-LATEST.mitab.zip \
	biogrid/BIOGRID-ALL-LATEST.mitab.zip \
	biogrid/$(BGFP)/BIOGRID-IDENTIFIERS-LATEST.tab.zip \
	biogrid/BIOGRID-IDENTIFIERS-LATEST.tab.zip

biogrid/: ; mkdir $@

biogrid/$(BGFP)/BIOGRID-ALL-LATEST.mitab.zip: FORCE
	$(CDBOG) $(WGET) $(FULLPTH) $(BGDL)/BIOGRID-ALL-LATEST.mitab.zip
biogrid/BIOGRID-ALL-LATEST.mitab.zip: biogrid/$(BGFP)/BIOGRID-ALL-LATEST.mitab.zip
	# Last-modified header missing workaround
	$(COPYCHANGED)
biogrid/$(BGFP)/BIOGRID-IDENTIFIERS-LATEST.tab.zip: FORCE
	$(CDBOG) $(WGET) $(FULLPTH) $(BGDL)/BIOGRID-IDENTIFIERS-LATEST.tab.zip
biogrid/BIOGRID-IDENTIFIERS-LATEST.tab.zip: biogrid/$(BGFP)/BIOGRID-IDENTIFIERS-LATEST.tab.zip
	# Last-modified header missing workaround
	$(COPYCHANGED)

biogrid_clean: ;  $(RM) biogrid/*
##########################################
CVFTP = ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar
clinvar: clinvar/ \
	clinvar/ClinVarFullRelease_00-latest.xml.gz \
	clinvar/gene_condition_source_id

clinvar/: ; mkdir $@

clinvar/ClinVarFullRelease_00-latest.xml.gz: FORCE
	cd clinvar; $(WGET) $(CVFTP)/xml/ClinVarFullRelease_00-latest.xml.gz
clinvar/gene_condition_source_id: FORCE
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

ctd/CTD_chemicals_diseases.tsv.gz: FORCE
	cd ctd; $(WGET) http://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz

ctd_clean: ; $(RM) ctd/*
##########################################
# included here for dipper/scripts/
dipper: dipper/ dipper/scripts
dipper/: ; git clone https://github.com/monarch-initiative/dipper.git

dipper/scripts: FORCE
	cd dipper; git pull

dipper_clean: ; $(RM) dipper
##########################################
eco: eco/ \
	eco/gaf-eco-mapping.yaml

eco/: ; mkdir $@

eco/gaf-eco-mapping.txt: eco/
	cd eco/; $(WGET) $(OBO)/eco/gaf-eco-mapping.txt

eco/gaf-eco-mapping.yaml: eco/gaf-eco-mapping.txt
	awk -F'\t' 'BEGIN{print "---"} \
	/^[^#]/ && $$2=="Default"{print "\"" $$1 "\": \"" $$3  "\""} \
	/^[^#]/ && $$2!="Default"{print "\"" $$1 "-" $$2 "\": \"" $$3 "\""}' \
	$< | sort  -k2,2 -t' '> $@

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

# ordered by "popularity" == mention in dipper output
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
$(ENSRDF_TARGET): FORCE
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
#  hp-to-eom-mapping.tsv
# a dead end from the 'disco' days?
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

flybase/md5sum.txt: FORCE
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
	unlink $(notdir $@) ; ln -s $(FLYPRE)/$$fname $(notdir $@)

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
genereviews/NBKid_shortname_OMIM.txt: FORCE
	cd genereviews; $(WGET) $(GRDL)/NBKid_shortname_OMIM.txt
genereviews/GRtitle_shortname_NBKid.txt: FORCE
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
	go/idmapping_selected.tab.gz \
	go/gaf-eco-mapping.txt \
	go/gaf-eco-mapping.yaml

go/: ; mkdir $@

$(foreach species, $(GOASPC), go/$(species).gaf.gz): FORCE
	cd go/; $(WGET) $(GOADL)/annotations/$(notdir $@)
go/go-refs.json: FORCE
	cd go/; $(WGET) http://current.geneontology.org/metadata/go-refs.json
go/idmapping_selected.tab.gz:  FORCE  # expensive
	cd go/; $(WGET) $(FTPEBI)/$(UPCRKB)/idmapping/idmapping_selected.tab.gz
go/gaf-eco-mapping.txt: eco/gaf-eco-mapping.txt
	# unlink $@;cd go/; ln -s ../$< $(notdir $@)
	$(SYMLINK)
go/gaf-eco-mapping.yaml: eco/gaf-eco-mapping.yaml
	# unlink $@;cd go/; ln -s ../$< $(notdir $@)
	$(SYMLINK)

go_clean: ; $(RM) go/*
##########################################
GWASFTP = ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest
GWASFILE = gwas-catalog-associations_ontology-annotated.tsv
gwascatalog: gwascatalog/ gwascatalog/$(GWASFILE)

gwascatalog/: ; mkdir $@
gwascatalog/$(GWASFILE): FORCE
	cd gwascatalog; $(WGET) $(GWASFTP)/$(GWASFILE)
# SO & MONDO ontologies need their own ontology cache
gwascatalog_clean:  ;  $(RM) gwascatalog/*
##########################################
EBIFTP := ftp://ftp.ebi.ac.uk
EBIPTH := pub/databases/genenames/new/tsv

hgnc: hgnc/ \
		hgnc/hgnc_complete_set.txt

hgnc/: ; mkdir $@
hgnc/hgnc_complete_set.txt: FORCE
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
hpoa/phenotype.hpoa: FORCE
	cd hgnc; $(WGET) $(HPOADL2)/phenotype.hpoa

hpoa_clean:  ; $(RM) hpoa/*
##########################################
IMPCDL = ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv
impc: impc/  impc/checksum.md5 \
	  impc/ALL_genotype_phenotype.csv.gz

impc/: ; mkdir $@
impc/checksum.md5: FORCE
	cd impc; $(WGET) $(IMPCDL)/checksum.md5
impc/ALL_genotype_phenotype.csv.gz: impc/checksum.md5
	cd impc; $(WGET) $(IMPCDL)/ALL_genotype_phenotype.csv.gz

impc_clean:  ;  $(RM) impc/*
##########################################
#/list
KEGGG = http://rest.genome.jp
# link
KEGGK = http://rest.kegg.jp

KGFP = list/disease \
    list/pathway \
    list/orthology \
    link/disease/omim \
    link/omim/hsa \
    list/hsa

KKFP = link/orthology/mmu \
    link/orthology/rno \
    link/orthology/dme \
    link/orthology/dre \
    link/orthology/cel \
    link/pathway/pubmed \
    link/pathway/ds \
    link/pathway/ko  \
    link/orthology/hsa \
    link/pathway/hsa \
    link/disease/hsa

kegg: kegg/ \
	kegg/list/disease \
	kegg/list/pathway \
	kegg/list/orthology \
	kegg/link/disease/omim \
	kegg/link/omim/hsa \
	kegg/list/hsa \
	kegg/link/orthology/mmu \
	kegg/link/orthology/rno \
	kegg/link/orthology/dme \
	kegg/link/orthology/dre \
	kegg/link/orthology/cel \
	kegg/link/pathway/pubmed \
	kegg/link/pathway/ds \
	kegg/link/pathway/ko \
	kegg/link/orthology/hsa \
	kegg/link/pathway/hsa \
	kegg/link/disease/hsa \
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

#kg_orig: $(foreach file, $(KGFP), kegg/$(file))
#	for fp in $^; do $(WGET) $(FULLPTH) $(KEGGG)/$$fp; done
#kk_orig: $(foreach file, $(KKFP), kegg/$(file))
#	for fp in $^; do $(WGET) $(FULLPTH) $(KEGGK)/$$fp; done

kegg/list/disease :  FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGG)/$(subst kegg/,,$@)
kegg/list/pathway : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGG)/$(subst kegg/,,$@)
kegg/list/orthology : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGG)/$(subst kegg/,,$@)
kegg/link/disease/omim : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGG)/$(subst kegg/,,$@)
kegg/link/omim/hsa : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGG)/$(subst kegg/,,$@)
kegg/list/hsa : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGG)/$(subst kegg/,,$@)
kegg/link/orthology/mmu : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/orthology/rno : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/orthology/dme : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/orthology/dre : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/orthology/cel : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/pathway/pubmed: FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/pathway/ds : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/pathway/ko : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/orthology/hsa : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/pathway/hsa : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)
kegg/link/disease/hsa : FORCE
	cd kegg; $(WGET) $(FULLPTH) $(KEGGK)/$(subst kegg/,,$@)

# note choosing native name except when there is a conflict  (hsa)
# conflicts would interfer with --timestamping
#   ... if their server supported Last-modified headers

kegg/disease:  kegg/list/disease
	$(COPYCHANGED)
kegg/pathway: kegg/list/pathway
	$(COPYCHANGED)
kegg/hsa_genes: kegg/list/hsa
	$(COPYCHANGED)
kegg/orthology: kegg/list/orthology
	$(COPYCHANGED)
kegg/disease_gene: kegg/link/disease/hsa
	$(COPYCHANGED)
kegg/omim: kegg/link/disease/omim
	$(COPYCHANGED)
kegg/omim2gene: kegg/link/omim/hsa
	$(COPYCHANGED)
kegg/ncbi: kegg/conv/ncbi-geneid/hsa
	$(COPYCHANGED)
kegg/human_gene2pathway: kegg/link/pathway/hsa
	$(COPYCHANGED)
kegg/hsa_orthologs: kegg/link/orthology/hsa
	$(COPYCHANGED)
kegg/mmu: kegg/link/orthology/mmu
	$(COPYCHANGED)
kegg/rno: kegg/link/orthology/rno
	$(COPYCHANGED)
kegg/dme: kegg/link/orthology/dme
	$(COPYCHANGED)
kegg/dre: kegg/link/orthology/dre
	$(COPYCHANGED)
kegg/cel: kegg/link/orthology/cel
	$(COPYCHANGED)
kegg/pubmed: kegg/link/pathway/pubmed
	$(COPYCHANGED)
kegg/ds: kegg/link/pathway/ds
	$(COPYCHANGED)
kegg/ko: kegg/link/pathway/ko
	$(COPYCHANGED)

kegg_clean:  ;  $(RM) kegg/*
##########################################
# pulls via sql queries
#mgi: mgi/
#mgi/; mkdir $@
##########################################
mmrrc: mmrrc/ mmrrc/mmrrc_catalog_data.csv
mmrrc/: ; mkdir $@
mmrrc/mmrrc_catalog_data.csv: FORCE
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
	$(CDMC) unlink 9606cytoBand.txt.gz ; ln -s hg19/cytoBand.txt.gz 9606cytoBand.txt.gz
monochrom/hg19/: ; mkdir $@
monochrom/hg19/cytoBand.txt.gz: FORCE
	cd monochrom/hg19; $(WGET) $(MCDL)/hg19/database/cytoBand.txt.gz

monochrom/10090cytoBand.txt.gz: monochrom/mm10/$(CBI)  monochrom/mm10/
	$(CDMC) unlink 10090cytoBand.txt.gz ;\
	ln -s mm10/$(CBI) 10090cytoBand.txt.gz  # note dropping 'Ideo' (to fix)
monochrom/mm10/: ; mkdir $@
monochrom/mm10/$(CBI): FORCE
	cd monochrom/mm10/ ; $(WGET) $(MCDL)/mm10/database/$(CBI)

monochrom/7955cytoBand.txt.gz: monochrom/danRer10/$(CBI) monochrom/danRer10/
	$(CDMC) unlink 7955cytoBand.txt.gz ; ln -s  danRer10/$(CBI) 7955cytoBand.txt.gz
monochrom/danRer10/: ; mkdir $@
monochrom/danRer10/$(CBI): FORCE
	cd monochrom/danRer10/ ;  $(WGET) $(MCDL)/danRer10/database/$(CBI)
monochrom/10116cytoBand.txt.gz: monochrom/rn6/$(CBI)  monochrom/rn6/
	$(CDMC) unlink 10116cytoBand.txt.gz ; ln -s rn6/$(CBI) 10116cytoBand.txt.gz
monochrom/rn6/: ; mkdir $@
monochrom/rn6/$(CBI): FORCE
	cd monochrom/rn6/; $(WGET) $(MCDL)/rn6/database/$(CBI)

monochrom/bosTau7cytoBand.txt.gz: monochrom/bosTau7/$(CBI) monochrom/bosTau7/
	$(CDMC) unlink bosTau7cytoBand.txt.gz; ln -s bosTau7/$(CBI) bosTau7cytoBand.txt.gz
monochrom/bosTau7/: ; mkdir $@

monochrom/bosTau7/$(CBI): FORCE
	cd monochrom/bosTau7/; $(WGET) $(MCDL)/bosTau7/database/$(CBI)

monochrom/galGal4cytoBand.txt.gz: monochrom/galGal4/$(CBI) monochrom/galGal4/
	$(CDMC) unlink galGal4cytoBand.txt.gz; ln -s galGal4/$(CBI) galGal4cytoBand.txt.gz
monochrom/galGal4/: ; mkdir $@
monochrom/galGal4/cytoBandIdeo.txt.gz: FORCE
	cd monochrom/galGal4/; $(WGET) $(MCDL)/galGal4/database/$(CBI)

monochrom/susScr3cytoBand.txt.gz: monochrom/susScr3/$(CBI) monochrom/susScr3/
	$(CDMC) unlink susScr3cytoBand.txt.gz ; ln -s susScr3/$(CBI) susScr3cytoBand.txt.gz
monochrom/susScr3/: ; mkdir $@
monochrom/susScr3/cytoBandIdeo.txt.gz: FORCE
	cd monochrom/susScr3/; $(WGET) $(MCDL)/susScr3/database/$(CBI)

monochrom/oviAri3cytoBand.txt.gz: monochrom/oviAri3/$(CBI) monochrom/oviAri3/
	$(CDMC) unlink oviAri3cytoBand.txt.gz ; ln -s oviAri3/$(CBI) oviAri3cytoBand.txt.gz
monochrom/oviAri3/: ; mkdir $@
monochrom/oviAri3/cytoBandIdeo.txt.gz: FORCE
	cd monochrom/oviAri3/; $(WGET) $(MCDL)/oviAri3/database/$(CBI)

monochrom/equCab2cytoBand.txt.gz: monochrom/equCab2/$(CBI) monochrom/equCab2/
	$(CDMC) unlink equCab2cytoBand.txt.gz ; ln -s equCab2/$(CBI) equCab2cytoBand.txt.gz
monochrom/equCab2/: ; mkdir $@
monochrom/equCab2/cytoBandIdeo.txt.gz: FORCE
	cd monochrom/equCab2/; $(WGET) $(MCDL)/equCab2/database/$(CBI)

monochrom_clean: ;  $(RM) monochrom/*
##########################################
MPDDL = https://phenomedoc.jax.org/MPD_downloads
mpd: mpd/ \
	mpd/ontology_mappings.csv \
	mpd/straininfo.csv \
	mpd/measurements.csv \
	mpd/strainmeans.csv.gz
# 	mpd_datasets_metadata.xml.gz  # TODO should this remain unused?

mpd/: ; mkdir $@
mpd/ontology_mappings.csv: FORCE
	cd mpd; $(WGET) $(MPDDL)/ontology_mappings.csv
mpd/straininfo.csv: FORCE
	cd mpd; $(WGET) $(MPDDL)/straininfo.csv
mpd/measurements.csv: FORCE
	cd mpd; $(WGET) $(MPDDL)/measurements.csv
mpd/strainmeans.csv.gz: FORCE
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

ncbigene/gene_info.gz: FORCE
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene_info.gz
	@ # test if in the future (GMT?)
	@ # touch -r "gene_info.gz" -d '-8 hour' "gene_info.gz"
ncbigene/gene_history.gz: FORCE
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene_history.gz
ncbigene/gene2pubmed.gz: FORCE
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene2pubmed.gz
ncbigene/gene_group.gz: FORCE
	cd ncbigene/ ; $(WGET) $(NCBIFTP)/gene_group.gz

GENEINFO = $(NCBIFTP)/GENE_INFO

ncbigene/Gallus_gallus.gene_info.gz: FORCE
	cd ncbigene/ ; $(WGET) $(GENEINFO)/Non-mammalian_vertebrates/Gallus_gallus.gene_info.gz
ncbigene/Sus_scrofa.gene_info.gz: FORCE
	cd ncbigene/ ; $(WGET) $(GENEINFO)/Mammalia/Sus_scrofa.gene_info.gz
ncbigene/Bos_taurus.gene_info.gz: FORCE
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
omia/omia.xml.gz: FORCE
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
orphanet/en_product6.xml: FORCE
	cd orphanet/; $(WGET) http://www.orphadata.org/data/xml/en_product6.xml

orphanet_clean: ; $(RM) orphanet/*

##########################################

owl: owl/ \
	owl/geno.owl \
	owl/faldo.ttl \
	owl/eco.owl \
	owl/iao.owl \
	owl/sepio.owl \
	owl/ero.owl \
	owl/pw.owl \
	owl/oban_core.ttl \
	owl/pco.owl \
	owl/xco.owl \
	owl/foaf.rdf \
	owl/dcelements.rdf \
	owl/metazoa.owl \
	owl/clo_core.owl \
    owl/monarch-merged.owl

owl/: ; mkdir $@

owl/geno.owl: FORCE
	cd owl; $(WGET) $(GITMON)/GENO-ontology/develop/src/ontology/geno.owl
owl/sepio.owl: FORCE
	cd owl; $(WGET) $(GITMON)/SEPIO-ontology/master/src/ontology/sepio.owl
owl/faldo.ttl: FORCE
	cd owl; $(WGET) $(GITRAW)/OBF/FALDO/master/faldo.ttl
owl/eco.owl: FORCE
	cd owl; $(WGET) $(OBO)/eco.owl
owl/iao.owl: FORCE
	cd owl; $(WGET) $(OBO)/iao.owl
owl/ero.owl: FORCE
	cd owl; $(WGET) $(OBO)/ero.owl
owl/pw.owl: FORCE
	cd owl; $(WGET) $(OBO)/pw.owl
owl/oban_core.ttl: FORCE
	cd owl; $(WGET) $(GITRAW)/jamesmalone/OBAN/master/ontology/oban_core.ttl
owl/pco.owl: FORCE
	cd owl; $(WGET) $(OBO)/pco.owl
owl/xco.owl: FORCE
	cd owl; $(WGET) $(OBO)/xco.owl


owl/foaf.rdf: FORCE
	cd owl; $(WGET) http://xmlns.com/foaf/spec/index.rdf -O foaf.rdf
	# Last-modified header missing workaround

owl/dcelements.rdf:  ## had local name dc.rdf
	cd owl; $(WGET) https://dublincore.org/2012/06/14/dcelements.rdf
########
# these next few are far too circular
owl/metazoa.owl: FORCE
	cd owl; $(WGET) https://data.monarchinitiative.org/owl/metazoa.owl
owl/clo_core.owl: FORCE
	cd owl; $(WGET) https://data.monarchinitiative.org/owl/clo_core.owl
# this one shouls be getting rebuilt based on if any if the previous are new
owl/monarch-merged.owl: FORCE
	cd owl; $(WGET) https://data.monarchinitiative.org/owl/monarch-merged.owl

owl_clean: ; $(RM) owl/*
########################################################
PNTH = ftp://ftp.pantherdb.org/ortholog
PNTHDL = ftp://ftp.pantherdb.org/ortholog/current_release

panther: panther/ \
	panther/current_release.ver \
	panther/RefGenomeOrthologs.tar.gz \
	panther/Orthologs_HCOP.tar.gz

panther/: ; mkdir $@

panther/current_release.ver: panther/
	/usr/bin/curl -s $(PNTH)/ |\
	/bin/sed -n 's/.*current_release -> \([0-9.]\+\)/\1/p' > $@

panther/RefGenomeOrthologs.tar.gz: FORCE
	cd panther; $(WGET) $(PNTHDL)/RefGenomeOrthologs.tar.gz
panther/Orthologs_HCOP.tar.gz: FORCE
	cd panther; $(WGET) $(PNTHDL)/Orthologs_HCOP.tar.gz

panther_clean: ; $(RM) panther/*
##########################################
RCTDL = https://www.reactome.org/download/current
reactome: reactome/ \
		reactome/Ensembl2Reactome.txt \
		reactome/ChEBI2Reactome.txt \
		reactome/gaf-eco-mapping.txt \
		reactome/gaf-eco-mapping.yaml

reactome/: ; mkdir $@
reactome/Ensembl2Reactome.txt: FORCE
	cd reactome; $(WGET) $(RCTDL)/Ensembl2Reactome.txt
reactome/ChEBI2Reactome.txt: FORCE
	cd reactome; $(WGET) $(RCTDL)/ChEBI2Reactome.txt

reactome/gaf-eco-mapping.txt:  eco/gaf-eco-mapping.txt
 	#cd reactome; unlink $(notdir $@); ln -s ../$< $(notdir $@)
	$(SYMLINK)
reactome/gaf-eco-mapping.yaml: eco/gaf-eco-mapping.yaml
	# unlink $@;cd reactome/; ln -s ../$< $(notdir $@)
	$(SYMLINK)
# curl -X GET "https://reactome.org/ContentService/data/diseases/doid" -H  "accept: text/plain"

reactome_clean: ; $(RM) rectome/*
##########################################

RGDFTP = ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology
rgd: rgd/ \
	rgd/rattus_genes_mp

rgd/: ; mkdir $@
rgd/rattus_genes_mp: FORCE
	cd rgd/; $(WGET) $(RGDFTP)/rattus_genes_mp

rgd_clean: ; $(RM) rdg/*
##########################################
SGDDL = https://downloads.yeastgenome.org/curation/literature
sgd: sgd/ sgd/phenotype_data.tab

sgd/: ; mkdir $@
sgd/phenotype_data.tab: FORCE
	cd sgd; $(WGET) $(SGDDL)/phenotype_data.tab
sgd_clean: ; $(RM) sgd/*
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

STRTAX = 9606 10090 7955 7227 6239 4932 10116
STRSPC = celegans fly human mouse yeast zebrafish
# no rat
STRFP = api/tsv-no-header

string: string/ \
		string/$(STRFP)/version \
		string/version \
		$(foreach txid, $(STRTAX), string/$(txid).$(SRTPTH).txt.gz) \
		$(foreach species, $(STRSPC), string/$(species).entrez_2_string.$(STRYR).tsv.gz)

string/: ; mkdir $@

string/$(STRFP)/version: FORCE
    # Last-modified header missing workaround
	cd string ; $(WGET) $(FULLPTH) $(STRING)/$(STRFP)/version
string/version: string/$(STRFP)/version
	$(COPYCHANGED); \
	if [ "$(STRVER)" != "$$(cut -f 1 version)" ] ;then echo "NEW VERSION of STRING!" ; \
	else  echo "same version of STRING" ; fi

$(foreach txid, $(STRTAX), string/$(txid).$(SRTPTH).txt.gz): FORCE
	cd string; $(WGET) $(STRDL)/$(SRTPTH)/$(notdir $@)
	#cd string; $(WGET) $(STRDL)/$(SRTPTH)/$(subst string/,,$@)

$(foreach species, $(STRSPC), string/$(species).entrez_2_string.$(STRYR).tsv.gz): FORCE
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
		wormbase/pub_xrefs.txt \
		wormbase/gaf-eco-mapping.yaml

wormbase/: ; mkdir $@
wormbase/CHECKSUMS: FORCE
	$(CDWB) $(WGET) $(WBFTP)/$(WBPROD)/$(notdir $@)
wormbase/letter: wormbase/CHECKSUMS
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)); \
	$(WGET) $(WBFTP)/$(WBPROD)/letter.$$wsnum ; \
	ln -s letter.$$wsnum $(notdir $@)
wormbase/c_elegans.PRJNA13758.geneIDs.txt.gz:  wormbase/CHECKSUMS
	#species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS273.geneIDs.txt.gz
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/species/$(WBSPC)/annotation/$(subst /,.,$(WBSPC)).$$wsnum.geneIDs.txt.gz;\
	ln -s $(WBPROD)/species/$(WBSPC)/annotation/$(subst /,.,$(WBSPC)).$$wsnum.geneIDs.txt.gz $(notdir $@)
wormbase/c_elegans.PRJNA13758.annotations.gff3.gz: wormbase/CHECKSUMS
	# species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS273.annotations.gff3.gz
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/species/$(WBSPC)/$(subst /,.,$(WBSPC)).$$wsnum.annotations.gff3.gz;\
	ln -s $(WBPROD)/species/$(WBSPC)/$(subst /,.,$(WBSPC)).$$wsnum.annotations.gff3.gz $(notdir $@)
wormbase/c_elegans.PRJNA13758.xrefs.txt.gz: wormbase/CHECKSUMS
	# species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS273.xrefs.txt.gz
	unlink $@ ; $(CDWB) wsnum=$$($(WSNUM)) ; \
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/species/$(WBSPC)/annotation/$(subst /,.,$(WBSPC)).$$wsnum.xrefs.txt.gz;\
	ln -s $(WBPROD)/species/$(WBSPC)/annotation/$(subst /,.,$(WBSPC)).$$wsnum.xrefs.txt.gz $(notdir $@)
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
	$(WGET) $(FULLPTH) $(WBFTP)/$(WBPROD)/ONTOLOGY/disease_association.$$wsnum.wb;\
	ln -s $(WBPROD)/ONTOLOGY/disease_association.$$wsnum.wb $(notdir $@)
# api call so no date or file version
wormbase/pub_xrefs.txt:
	$(CDWB) $(WGET) -O $(notdir $@) \
	http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref
	# Last-modified header missing workaround
wormbase/gaf-eco-mapping.yaml: eco/gaf-eco-mapping.yaml
	# unlink $@;cd wormbase/; ln -s ../$< $(notdir $@)
	$(SYMLINK)

wormbase_clean: ;  $(RM) wormbase/*
##########################################
ZFDL = https://zfin.org/downloads
# not        https://purl.obolibrary.org/obo/zp/src/curation
ZPCURATION = $(OBO)/zp
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

$(addprefix zfin/, $(ZFFILES)): FORCE
	cd zfin/; $(WGET) $(ZFDL)/$(notdir $@)

zfin/id_map_zfin.tsv: FORCE
	cd zfin/; $(WGET) $(ZPCURATION)/$(notdir $@)

zfin_clean: ;  $(RM) zfin/*
##########################################
zfinslim:  zfin/ zfinslim/ \
		zfinslim/phenoGeneCleanData_fish.txt \
		zfinslim/id_map_zfin.tsv

zfinslim/: ; mkdir $@

zfinslim/phenoGeneCleanData_fish.txt: zfin/phenoGeneCleanData_fish.txt
	#cd  $(dir $@); id [ !L "$(notdir $@)" ] ; then ln -s "$(notdir $@)" "$<"; fi
	$(SYMLINK)
zfinslim/id_map_zfin.tsv: zfin/id_map_zfin.tsv
	# cd  $(dir $@); id [ ! L "$(notdir $@)" ] ; then ln -s "$(notdir $@)" "$<"; fi
	$(SYMLINK)

zfinslim_clean: ; $(RM) zfin/id_map_zfin.tsv zfin/phenoGeneCleanData_fish.txt
##########################################


clean: $(foreach target, $(SOURCES), $(target)_clean)
