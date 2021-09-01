# Checklist difference

Two rewrites in progress.  See the 'alignment' branch, which I have
abandoned for now, and the [list
tools](https://github.com/jar398/listtools) repository, which is
intended to be a code base shared between ASU BioKIC Taxon Concepts
projects and Encyclopedia of Life.





The inputs are two checklists in TSV or CSV format (extension .tsv or
.csv).  The output is an ad hoc report file.

Examples
 * [NCBI Taxonomy 2015 to 2020](doc/ncbi-2015-2020.csv) (Primates only)
 * [NCBI Taxonomy 2020 to GBIF](doc/ncbi-gbif.csv) (Primates only)

Documentation:

* [Method](doc/method.md)
* [User guide](doc/user-guide.md)


## Example: NCBI Primates vs. GBIF Primates

You can run a similar example (comparing two versions of NCBI
Taxonomy) by simply saying `make`.  But in detail, the procedure has
multiple steps.

### Get NCBI taxonomy from FTP site

We'll put everything related to the February 2020 version of NCBI
taxonomy under `work/ncbi/2020-01-01`.  Start with the release (`dump`).

    mkdir -p work/ncbi/2020-01-01/dump
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-01-01.zip
    unzip -d work/ncbi/2020-01-01/dump taxdmp_2020-01-01.zip

Of course you can do this with any version of NCBI you like, by substituting the date.

### Get GBIF taxonomy from GBIF site

Similarly `work/gbif/2019-09-16`.  The release is a DwCA file (thus `dwca`).

    mkdir -p work/gbif/2019-09-16/dwca
    wget http://rs.gbif.org/datasets/backbone/2019-09-06/backbone.zip
    unzip -d work/gbif/2019-09-16/dwca -q backbone.zip

For futher information see [the backbone taxonomy landing
page](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c).

Other archived versions of GBIF are available on their site.

### Compare them

    python3 src/cldiff.py work/gbif/2019-09-16/primates.csv \
      work/ncbi/2020-01-01/primates.csv --out diff.out

The first checklist (or taxonomy) is the "A checklist" and the second is
the "B checklist".

If you put the two in the opposite report you'll get the comparison
ordered by the other taxonomy:

    python3 src/cldiff.py work/gbif/2019-09-16/primates.csv work/ncbi/2020-01-01/primates.csv
