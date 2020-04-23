WORK=work

all: $(WORK)/output.csv

# A (NCBI 2015)

A=$(WORK)/ncbi/2015-01-01
$(A)/dump.zip:
	mkdir -p $(A)/dump
	wget -O $@ ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-01-01.zip
$(A)/dump/names.dmp: $(A)/dump.zip
	mkdir -p $(A)/dump
	unzip -d $(A)/dump $(A)/dump.zip
	touch $(A)/dump/*
# Convert to DwC form
$(A)/converted.csv: src/ncbi_to_dwc.py $(A)/dump/names.dmp
	python3 src/ncbi_to_dwc.py $(A)/dump --out $@

# B (NCBI 2020)

B=$(WORK)/ncbi/2020-01-01
$(B)/dump.zip:
	mkdir -p $(B)/dump
	wget -O $@ ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-01-01.zip
$(B)/dump/names.dmp: $(B)/dump.zip
	mkdir -p $(B)/dump
	unzip -d $(B)/dump $(B)/dump.zip
	touch $(B)/dump/*
# Convert to DwC form
$(B)/converted.csv: src/ncbi_to_dwc.py $(B)/dump/names.dmp
	python3 src/ncbi_to_dwc.py $(B)/dump --out $@

# C (GBIF)

C=$(WORK)/gbif/2019-09-16
$(C)/backbone.zip:
	mkdir -p $(C)/dump
	wget -O $@ http://rs.gbif.org/datasets/backbone/2019-09-06/backbone.zip
$(C)/dump/Taxon.tsv: $(C)/backbone.zip
	mkdir -p $(C)/dump
	unzip -d $(C)/dump $<
	touch $(C)/dump/*
# No conversion needed

# Extract Primates from each

$(A)/primates.csv: $(A)/converted.csv
	python3 src/subset_dwc.py $< 9443 --out $@

$(B)/primates.csv: $(B)/converted.csv
	python3 src/subset_dwc.py $< 9443 --out $@

$(C)/primates.csv: $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $< 798 --out $@
foo: $(C)/primates.csv

$(WORK)/diff.csv: src/cldiff.py $(A)/primates.csv $(B)/primates.csv
	python3 src/cldiff.py $(A)/primates.csv $(B)/primates.csv --out $@.new
	mv $@.new $@

$(WORK)/output.csv: src/diff.py $(A)/primates.csv $(B)/primates.csv
	python3 src/diff.py $(A)/primates.csv $(B)/primates.csv --out $@.new
	mv $@.new $@
