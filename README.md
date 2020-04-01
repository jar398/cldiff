# Checklist difference

The 'checklist diff' idea developed out of conversations between Nico, Beckett, and Jonathan in Feb 2020.  I thought it might be fun to write a quick and dirty prototype as way to help develop a specification for what we might want eventually.

I think there are probably several use cases, e.g.:
* General user of checklists (especially from GBIF) who wants to know how they compare
* Someone interested in taxonomy (or 'taxon concepts') who wants a starting point for developing a an RCC-5 alignment of two checklists or taxonomies

The inputs are two checklists (or taxonomies) in Darwin Core Archive
format, that have been unzipped to directories.  The output is an ad
hoc report file.  See below for restrictions.

## Example: NCBI Primates vs. GBIF Primates

### Get NCBI taxonomy from FTP site

We'll put everything related to the February 2020 version of NCBI
taxonomy under `ncbi/2020-01-01`.  Start with the release (`dump`).

    mkdir -p ncbi/2020-01-01/dump
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-01-01.zip
    unzip -d ncbi/2020-01-01/dump taxdmp_2020-01-01.zip

Of course you can do this with any version you like, by substituting the date.

### Get GBIF taxonomy from GBIF site

Similarly `gbif/2019-09-16`.  The release is a DwCA file (`dwca`).

    mkdir -p gbif/2019-09-16/dwca
    wget http://rs.gbif.org/datasets/backbone/2019-09-06/backbone.zip
    unzip -d gbif/2019-09-16/dwca -q backbone.zip

For futher information see [the backbone taxonomy landing
page](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c).

Other archived versions of GBIF are available on their site.

### Convert NCBI taxonomy to Darwin Core form

NCBI has its own taxonomy dump format, which needs to be converted to
DwCA.  There's a tool for this purpose:

    python3 src/ncbi_to_dwca.py ncbi/2020-01-01/dump ncbi/2020-01-01/dwca

(The GBIF files are DwCA format already, so they can be used directly.)

### From each, extract Primates only 

Dealing with the whole taxonomies would be time consuming and painful
to analyze, so it's best to just look at reasonably sized taxa.
The `subset_dwca.py` tool extracts just the part we want from the whole
taxonomy's DwCA.

    python3 src/subset_dwca.py gbif/2019-09-16/dwca gbif/2019-09-16/dwca 798 gbif/2019-09-16/primates
    python3 src/subset_dwca.py ncbi/2020-01-01/dwca ncbi/2020-01-01/dwca 9443 ncbi/2020-01-01/primates

You can find the identifiers for the taxon of interest via each site's
web interface.

### Compare them

    python3 src/cldiff.py gbif/2019-09-16/primates ncbi/2020-01-01/primates

The first checklist (or taxonomy) is the "A checklist" and the second is
the "B checklist".

The output file (the comparison report) defaults to 'diff.csv'.  To
specify the output file name, use the `--out` parameter: `python3
src/cldiff.py A B --out out.csv`

The output file is in CSV format with a header row.

If you put the two in the opposite report you'll get the comparison
ordered by the other taxonomy:

    python3 src/cldiff.py gbif/2019-09-16/primates ncbi/2020-01-01/primates

## Report file format

This is highly in flux... as of today we have:

 * `nesting` - increases with hierarchical nesting depth, so in
   general (for example) the value in this row will be greater for
   species than for families.
 * `A_id` - TNU id of the TNU in the A taxonomy
 * `A_name` - canonical name of the TNU in the A taxonomy
 * `B_id` - TNU id of the TNU in the B taxonomy
 * `B_name` - canonical name of the TNU in the B taxonomy
 * `how` - describes matching status, usually via topology
 * `mode` - describes matching via synonymy

The row order follows the hierarchy of A.  Where there is a graft of
an unmatched subtree of B into the A hierarchy, the B grafted subtree
is shown according to the B hierarchy.

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

eventually: epithet and a few others

tbd: handle incertae sedis / unclassified containers sensibly
