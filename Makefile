WORK=work

SOURCES=src/report.py src/alignment.py src/articulation.py src/relation.py \
        src/checklist.py src/diff.py src/merge.py src/report.py src/eulerx.py \
	src/intension.py src/table.py src/dribble.py src/property.py

all: $(WORK)/primates-ncbi-2015-2020.csv 

# N15 (NCBI 2015)

N15=$(WORK)/ncbi/2015-05-01
$(N15)/dump.zip:
	mkdir -p $(N15)/dump
	wget -O $@ ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-05-01.zip
$(N15)/dump/names.dmp: $(N15)/dump.zip
	mkdir -p $(N15)/dump
	unzip -d $(N15)/dump $(N15)/dump.zip
	touch $(N15)/dump/*
# Convert to DwC form
$(N15)/converted.csv: src/ncbi_to_dwc.py $(N15)/dump/names.dmp
	python3 src/ncbi_to_dwc.py $(N15)/dump --out $@

# N20 (NCBI 2020)

N20=$(WORK)/ncbi/2020-08-01
$(N20)/dump.zip:
	mkdir -p $(N20)/dump
	wget -O $@ ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-08-01.zip
$(N20)/dump/names.dmp: $(N20)/dump.zip
	mkdir -p $(N20)/dump
	unzip -d $(N20)/dump $(N20)/dump.zip
	touch $(N20)/dump/*
# Convert to DwC form
$(N20)/converted.csv: src/ncbi_to_dwc.py $(N20)/dump/names.dmp
	python3 src/ncbi_to_dwc.py $(N20)/dump --out $@

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

# --------------------
# Extract Primates from each, and compare

$(N15)/primates.csv: src/subset_dwc.py $(N15)/converted.csv
	python3 src/subset_dwc.py $(N15)/converted.csv 9443 --out $@

$(N20)/primates.csv: src/subset_dwc.py $(N20)/converted.csv
	python3 src/subset_dwc.py $(N20)/converted.csv 9443 --out $@

$(C)/primates.csv: src/subset_dwc.py $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $(C)/dump/Taxon.tsv 798 --out $@.new
	mv $@.new $@
p: $(C)/primates.csv

$(WORK)/primates-ncbi-2015-2020.csv: $(SOURCES) $(N15)/primates.csv $(N20)/primates.csv
	python3 src/report.py $(N15)/primates.csv \
	                      $(N20)/primates.csv \
	  --out $@.new 
	mv $@.new $@
	mv $@.new.log $@.log
#  --low-tag=N15 --high-tag=N20

$(WORK)/primates-ncbi-2015-2020.ex: $(SOURCES) $(N15)/primates.csv $(N20)/primates.csv
	python3 src/report.py $(N15)/primates.csv \
	                      $(N20)/primates.csv \
	  --out $@.new --format eulerx
	mv $@.new $@
	mv $@.new.log $@.log

$(WORK)/primates-ncbi-gbif.csv: $(SOURCES) $(N20)/primates.csv $(C)/primates.csv
	python3 src/report.py $(N20)/primates.csv \
			      $(C)/primates.csv --out $@.new
	mv $@.new $@
	mv $@.new.log $@.log

pri: $(WORK)/primates-ncbi-2015-2020.csv

# ----------------------------------------------------------------------
# Mammalia = NCBI 40674

$(N15)/mammalia.csv: src/subset_dwc.py $(N15)/converted.csv
	python3 src/subset_dwc.py $(N15)/converted.csv 40674 --out $@

$(N20)/mammalia.csv: src/subset_dwc.py $(N20)/converted.csv
	python3 src/subset_dwc.py $(N20)/converted.csv 40674 --out $@

$(WORK)/mammalia-ncbi-2015-2020.csv: $(SOURCES) $(N15)/mammalia.csv $(N20)/mammalia.csv
	python3 src/report.py $(N15)/mammalia.csv $(N20)/mammalia.csv \
	  --out $@.new 
	mv $@.new $@
	mv $@.new.log $@.log

mam: $(WORK)/mammalia-ncbi-2015-2020.csv

# --------------------

publish: doc/primates-ncbi-2015-2020.csv \
	 doc/primates-ncbi-2015-2020.ex \
	 doc/mammalia-ncbi-2015-2020.csv \
	 doc/primates-ncbi-gbif.csv

doc/primates-ncbi-2015-2020.csv: $(WORK)/primates-ncbi-2015-2020.csv
	cp -p $< $@
	cp -p $<.log $@.log
doc/primates-ncbi-2015-2020.ex: $(WORK)/primates-ncbi-2015-2020.ex
	cp -p $< $@

doc/mammalia-ncbi-2015-2020.csv: $(WORK)/mammalia-ncbi-2015-2020.csv
	cp -p $< $@
	cp -p $<.log $@.log

doc/primates-ncbi-gbif.csv: $(WORK)/primates-ncbi-gbif.csv
	cp -p $< $@
	cp -p $<.log $@.log

gbif: $(WORK)/primates-ncbi-gbif.csv

test:
	python3 src/report.py "(l(pab)c)" "(l(qab)c)" --out -

# ----------------------------------------------------------------------
# Other groups to play with

# Trillium

$(N15)/trillium.csv: src/subset_dwc.py $(N15)/converted.csv
	python3 src/subset_dwc.py $(N15)/converted.csv 49674 --out $@

$(N20)/trillium.csv: src/subset_dwc.py $(N20)/converted.csv
	python3 src/subset_dwc.py $(N20)/converted.csv 49674 --out $@

$(C)/trillium.csv: src/subset_dwc.py $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $(C)/dump/Taxon.tsv 2742182 --out $@

$(WORK)/trillium-ncbi-2015-2020.csv: $(SOURCES) $(N15)/trillium.csv $(N20)/trillium.csv
	python3 src/report.py $(N15)/trillium.csv $(N20)/trillium.csv \
	  --out $@.new 
	mv $@.new $@
t: $(WORK)/trillium-ncbi-2015-2020.csv


$(N15)/mag.csv: src/subset_dwc.py $(N15)/converted.csv
	python3 src/subset_dwc.py $(N15)/converted.csv 3401 --out $@

$(N20)/mag.csv: src/subset_dwc.py $(N20)/converted.csv
	python3 src/subset_dwc.py $(N20)/converted.csv 3401 --out $@

$(C)/mag.csv: src/subset_dwc.py $(C)/dump/Taxon.tsv
	python3 src/subset_dwc.py $(C)/dump/Taxon.tsv 4690 --out $@

$(WORK)/mag-ncbi-2015-2020.csv: $(SOURCES) $(N15)/mag.csv $(N20)/mag.csv
	python3 src/report.py $(N15)/mag.csv $(N20)/mag.csv \
	  --out $@.new 
	mv $@.new $@

mag: $(WORK)/mag-ncbi-2015-2020.csv

clean:
	rm -rf $(WORK)
