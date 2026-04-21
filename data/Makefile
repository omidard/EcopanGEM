DATA_DIR := data
ZIP := $(DATA_DIR)/Data.zip
URL := https://zenodo.org/records/17581962/files/Data.zip

.PHONY: fetch-data clean

fetch-data: $(ZIP)

$(DATA_DIR):
	mkdir -p $(DATA_DIR)

$(ZIP): | $(DATA_DIR)
	curl -L --fail -o $(ZIP) "$(URL)"

clean:
	rm -rf $(DATA_DIR)
