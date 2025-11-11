# Fetch Data.zip from Zenodo into data/
DATA_DIR := data
ZIP := $(DATA_DIR)/Data.zip
URL := https://zenodo.org/records/17581962
.PHONY: data clean

data: $(ZIP)

$(DATA_DIR):
	mkdir -p $(DATA_DIR)

$(ZIP): | $(DATA_DIR)
	curl -L --fail -o $(ZIP) "$(URL)"

clean:
	rm -rf $(DATA_DIR)
