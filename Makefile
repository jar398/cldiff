WORK=work

all: $(WORK)/ncbi-2015-2020.ex $(WORK)/ncbi-2015-2020.csv 

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

$(A)/primates.csv: src/subset_dwc.py $(A)/converted.csv
	python3 src/subset_dwc.py $(A)/converted.csv 9443 --out $@

$(B)/primates.csv: src/subset_dwc.py $(B)/converted.csv
	python3 src/subset_dwc.py $(B)/converted.csv 9443 --out $@

$(C)/primates.csv: src/subset_dwc.py $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $(C)/dump/Taxon.tsv 798 --out $@
p: $(C)/primates.csv

SOURCES=src/report.py src/alignment.py src/articulation.py src/relation.py src/checklist.py src/diff.py src/merge.py src/eulerx.py

$(WORK)/ncbi-2015-2020.csv: $(SOURCES) $(A)/primates.csv $(B)/primates.csv
	python3 src/report.py $(A)/primates.csv $(B)/primates.csv \
	  --out $@.new 
	mv $@.new $@

ex: $(WORK)/ncbi-2015-2020.ex
$(WORK)/ncbi-2015-2020.ex: $(SOURCES) $(A)/primates.csv $(B)/primates.csv
	python3 src/report.py $(A)/primates.csv $(B)/primates.csv \
	  --out $@.new --format eulerx
	mv $@.new $@

$(WORK)/ncbi-gbif.csv: $(SOURCES) $(B)/primates.csv $(C)/primates.csv
	python3 src/report.py $(B)/primates.csv --low-tag=B \
			      $(C)/primates.csv --high-tag=C --out $@.new
	mv $@.new $@

publish: doc/ncbi-2015-2020.csv doc/ncbi-gbif.csv

doc/ncbi-2015-2020.csv: $(WORK)/ncbi-2015-2020.csv
	cp -p $< $@

doc/ncbi-gbif.csv: $(WORK)/ncbi-gbif.csv
	cp -p $< $@

gbif: $(WORK)/ncbi-gbif.csv

test:
	python3 src/report.py "(l(pab)c)" "(l(qab)c)" --out -


# ----------------------------------------------------------------------

# Trillium

$(A)/trillium.csv: src/subset_dwc.py $(A)/converted.csv
	python3 src/subset_dwc.py $(A)/converted.csv 49674 --out $@

$(B)/trillium.csv: src/subset_dwc.py $(B)/converted.csv
	python3 src/subset_dwc.py $(B)/converted.csv 49674 --out $@

$(C)/trillium.csv: src/subset_dwc.py $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $(C)/dump/Taxon.tsv 2742182 --out $@

$(WORK)/trillium-ncbi-2015-2020.csv: $(SOURCES) $(A)/trillium.csv $(B)/trillium.csv
	python3 src/report.py $(A)/trillium.csv $(B)/trillium.csv \
	  --out $@.new 
	mv $@.new $@
t: $(WORK)/trillium-ncbi-2015-2020.csv


$(A)/mag.csv: src/subset_dwc.py $(A)/converted.csv
	python3 src/subset_dwc.py $(A)/converted.csv 3401 --out $@

$(B)/mag.csv: src/subset_dwc.py $(B)/converted.csv
	python3 src/subset_dwc.py $(B)/converted.csv 3401 --out $@

$(C)/mag.csv: src/subset_dwc.py $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $(C)/dump/Taxon.tsv 4690 --out $@

$(WORK)/mag-ncbi-2015-2020.csv: $(SOURCES) $(A)/mag.csv $(B)/mag.csv
	python3 src/report.py $(A)/mag.csv $(B)/mag.csv \
	  --out $@.new 
	mv $@.new $@

m: $(WORK)/mag-ncbi-2015-2020.csv

