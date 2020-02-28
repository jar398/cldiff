# cldiff
Checklist difference


## Comparing two checklists: example: NCBI vs. GBIF primates

### Get NCBI taxonomy from FTP site

    mkdir -p ncbi/2020-02-01/dump
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-02-01.zip
    unzip -d ncbi/2020-02-01/dump taxdmp_2020-02-01.zip

### Get GBIF taxonomy from GBIF site

The GBIF site doesn't expose the download URL.  You'll have to [go
there
yourself](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c)
and figure out what to do.

    mkdir -p gbif/2019-09-16/dwca
    wget something-something/backbone-current.zip
    unzip -d gbif/2019-09-16/dwca -q backbone-current.zip

### Convert NCBI taxonomy to Darwin Core form

    time python3 src/ncbi_to_dwca.py ncbi/2020-01-01/dump ncbi/2020-01-01/dwca

### From each, extract Primates only 

    python3 src/subset_dwca.py gbif/2019-09-16/dwca gbif/2019-09-16/dwca 798 gbif/2019-09-16/primates
    python3 src/subset_dwca.py ncbi/2020-01-01/dwca ncbi/2020-01-01/dwca 9443 ncbi/2020-01-01/primates

### Compare them

    python3 src/cldiff.py gbif/2019-09-16/primates ncbi/2020-01-01/primates
