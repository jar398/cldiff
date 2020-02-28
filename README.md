# Checklist difference

## Comparing two checklists: example: NCBI vs. GBIF primates

### Get NCBI taxonomy from FTP site

    mkdir -p ncbi/2020-02-01/dump
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-02-01.zip
    unzip -d ncbi/2020-02-01/dump taxdmp_2020-02-01.zip

### Get GBIF taxonomy from GBIF site

    mkdir -p gbif/2019-09-16/dwca
    wget https://files.opentreeoflife.org/gbif/gbif-20190916/gbif-20190916.zip
    unzip -d gbif/2019-09-16/dwca -q gbif-20190916.zip

The GBIF site doesn't expose the download URL.  You'll have to [go
there
yourself](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c)
and figure out what to do.  The downloaded file will be called
`backbone-current.zip`.

### Convert NCBI taxonomy to Darwin Core form

    time python3 src/ncbi_to_dwca.py ncbi/2020-01-01/dump ncbi/2020-01-01/dwca

### From each, extract Primates only 

    python3 src/subset_dwca.py gbif/2019-09-16/dwca gbif/2019-09-16/dwca 798 gbif/2019-09-16/primates
    python3 src/subset_dwca.py ncbi/2020-01-01/dwca ncbi/2020-01-01/dwca 9443 ncbi/2020-01-01/primates

You can find the identifiers for the taxon of interest via each site's
web interface.

### Compare them

    python3 src/cldiff.py gbif/2019-09-16/primates ncbi/2020-01-01/primates

The first checklist/taxonomy is the "A checklist" and the second is
the "B checklist".

The output file name is hardwired as 'diff.csv'.  The file is in CSV
format with a header row.

## Diff file format

This is going to be highly in flux... as of today we have:

 * `nesting` - increases with hierarchical nesting depth, so in
   general (for example) the value in this row will be greater for
   species than for families.
 * `A_id` - TNU id of the TNU in the A taxonomy
 * `A_name` - canonical name of the TNU in the A taxonomy
 * `B_id` - TNU id of the TNU in the B taxonomy
 * `B_name` - canonical name of the TNU in the B taxonomy
 * `how` - describes matching
 * `mode` - describes matching

The row order follows the hierarchy of A, except where there is a
graft of a subtree of B, in which case the B subtree follows the B
hierarchy.

## DwC fields of interest:

The inputs are unzipped Darwin Core archives (DwCA), but the
`meta.xml` file is not examined.  Instead the TNU file (usually called
something like taxon.tsv or Taxon.txt) only is found and is assumed to
have column headers of the sort used in GBIF DwCAs.  Only the
following columns/fields are used here:

 * `taxonID`  - this is a misnomer, it's a taxon usage identifier (TNU) id
 * `scientificName`  - if available; includes authority + year 
 * `canonicalName`   - without authority + year
 * `acceptedNameUsageID` - TNU id for corresponding accepted TNU
 * `parentNameUsageID`  - TNU id for parent TNU
 * `taxonomicStatus`  - value is `accepted` or something else
 * `nomenclaturalStatus` - for NCBI, gives synonym class

