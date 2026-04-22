# EcopanGEM data fetcher
# Downloads all data files from Zenodo record 17581962 and organizes them
# into the directory structure expected by the scripts in scripts/.
#
# Usage (from inside data/):
#   make fetch-data     # download + organize everything (~5.4 GB)
#   make clean          # remove all downloaded data
#
# Or from the repo root:
#   make -C data fetch-data

ZENODO := https://zenodo.org/records/17581962/files
CURL   := curl -L --fail --progress-bar

.PHONY: fetch-data clean

fetch-data: \
	ref_model_dir/marlbr2.mat \
	neighbors/locustags_genes_mapping.pkl \
	gapfilled_curated/.done \
	ecoli_gprs.csv \
	biolog.csv \
	biolog_panGEM.csv \
	Unique_ModelSEED_Reaction_Aliases.txt \
	header_to_allele.pickle \
	phylon_locustags_df.csv \
	curated_metadata_mash_filtered.pickle \
	pangenome.csv \
	iB21_1397.json \
	neighbors/all_reactions_gene_neighborhood.csv \
	pangenome_s/.done
	@echo ""
	@echo "========================================"
	@echo "  Data download and setup complete!"
	@echo "  See README.md for pipeline details."
	@echo "========================================"

# ---- Reference model ----
ref_model_dir/marlbr2.mat:
	mkdir -p ref_model_dir
	$(CURL) -o $@ "$(ZENODO)/marlbr2.mat?download=1"

# ---- Locus-tag mapping (2.4 GB) ----
neighbors/locustags_genes_mapping.pkl:
	mkdir -p neighbors
	$(CURL) -o $@ "$(ZENODO)/locustags_genes_mapping.pkl?download=1"

# ---- Curated GEMs (zipped) ----
gapfilled_curated/.done: Ecoli_GEMs_for_Complete_genomes.zip
	mkdir -p gapfilled_curated
	unzip -o -q $< -d gapfilled_curated/
	touch $@

Ecoli_GEMs_for_Complete_genomes.zip:
	$(CURL) -o $@ "$(ZENODO)/Ecoli_GEMs_for_Complete_genomes.zip?download=1"

# ---- GPR matrix (zipped) ----
ecoli_gprs.csv: ecoli_gprs.csv.zip
	unzip -o -q $<
	touch $@

ecoli_gprs.csv.zip:
	$(CURL) -o $@ "$(ZENODO)/ecoli_gprs.csv.zip?download=1"

# ---- Biolog phenotype data ----
biolog.csv:
	$(CURL) -o $@ "$(ZENODO)/biolog.csv?download=1"

# Scripts expect biolog_panGEM.csv; Zenodo has biolog.csv
biolog_panGEM.csv: biolog.csv
	cp $< $@

# ---- ModelSEED aliases ----
Unique_ModelSEED_Reaction_Aliases.txt:
	$(CURL) -o $@ "$(ZENODO)/Unique_ModelSEED_Reaction_Aliases.txt?download=1"

# ---- Allele mapping (zipped) ----
header_to_allele.pickle: header_to_allele.pickle.zip
	unzip -o -q $<
	touch $@

header_to_allele.pickle.zip:
	$(CURL) -o $@ "$(ZENODO)/header_to_allele.pickle.zip?download=1"

# ---- Phylon locus tags (zipped) ----
phylon_locustags_df.csv: phylon_locustags_df.csv.zip
	unzip -o -q $<
	touch $@

phylon_locustags_df.csv.zip:
	$(CURL) -o $@ "$(ZENODO)/phylon_locustags_df.csv.zip?download=1"

# ---- Metadata pickle ----
curated_metadata_mash_filtered.pickle:
	$(CURL) -o $@ "$(ZENODO)/curated_metadata_mash_filtered.pickle?download=1"

# ---- Pangenome matrix (zipped) ----
pangenome.csv: pangenome.csv.zip
	unzip -o -q $<
	touch $@

pangenome.csv.zip:
	$(CURL) -o $@ "$(ZENODO)/pangenome.csv.zip?download=1"

# ---- Universal model; scripts expect iB21_1397.json ----
iB21_1397.json: universal_model.json
	cp $< $@

universal_model.json:
	$(CURL) -o $@ "$(ZENODO)/universal_model.json?download=1"

# ---- Gene neighborhood (zipped) ----
neighbors/all_reactions_gene_neighborhood.csv: all_reactions_gene_neighborhood.csv.zip
	mkdir -p neighbors
	unzip -o -q $< -d neighbors/
	touch $@

all_reactions_gene_neighborhood.csv.zip:
	$(CURL) -o $@ "$(ZENODO)/all_reactions_gene_neighborhood.csv.zip?download=1"

# ---- Pangenome sensitivity data (zipped, ~1.6 GB) ----
pangenome_s/.done: pangenome_s.zip
	mkdir -p pangenome_s
	unzip -o -q $< -d pangenome_s/
	touch $@

pangenome_s.zip:
	$(CURL) -o $@ "$(ZENODO)/pangenome_s.zip?download=1"

# ---- Cleanup ----
clean:
	rm -rf ref_model_dir neighbors pangenome_s gapfilled_curated
	rm -f marlbr2.mat locustags_genes_mapping.pkl
	rm -f biolog.csv biolog_panGEM.csv universal_model.json iB21_1397.json
	rm -f ecoli_gprs.csv header_to_allele.pickle phylon_locustags_df.csv
	rm -f pangenome.csv curated_metadata_mash_filtered.pickle
	rm -f Unique_ModelSEED_Reaction_Aliases.txt
	rm -f *.zip
