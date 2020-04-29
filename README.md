# Checklist difference

The inputs are two checklists in TSV or CSV format (extension .tsv or
.csv).  The output is an ad hoc report file.

Examples
 * [NCBI Taxonomy 2015 to 2020](doc/ncbi-2015-2020.csv) (Primates only)
 * [NCBI Taxonomy 2020 to GBIF](doc/ncbi-gbif.csv) (Primates only)

See below for
restrictions.

A "checklist" is a table with one row per taxon name.  A "taxonomy" is
a checklist that provides a hierarchy (parent pointers).

## Report file format

This is highly in flux... as of today we have:

 * `operation` - what needs to be done to A to make it more like B
     * `NO CHANGE`
     * Intensional properies:
         * `CHANGE ID`   - change the taxon id of this node (similarly RANK)
         * `CHANGE NAME` - change the canonical name of this node
     * Comings and goings:
         * `GRAFT`   - add a B subtree to A
         * `REMOVE`  - remove an A subtree so it's not in B
         * `INSERT`  - insert a new B node into the A tree, with
                       descendants mapped to A
         * `ELIDE`   - splice this A node, one with B descendants, out of the tree
     * Reorganization:
         * `MOVE`    - this B subtree now has a parent mapped to a
                       different place in the A tree
         * `BREAK`   - remove an A node that has descendants mapped to B
         * `REFORM`  - add new B node as part of reorganization
         * `CHANGE PARENT` - change the parent; it is now in another place in the tree
     * `MULTIPLE` - this A node is the favorite of multiple B nodes
         * `OPTION`  - several B nodes have this A node as their best
                       match.  Could be either a genuine split, some kind
                       of error, or something else.
     * `...` - this node's children are all unchanged so they have
          not been listed.
     * Indentation level reflects parent/child relationships in A
 * `A node` - name of node in A
 * `RCC5` - one of `< > = >< !` for proper containment, properly
    contains, same extension, conflict (overlap without equality or
    containment), disjoint
 * `B node`
 * `steps` - the steps, as a sequence separated by '→', taken in
    finding a path from the A node to the 
    B node; alternatively, the relationship between the two nodes.
    Sometimes this column just holds a comment.

The row order is a preorder traversal of the hierarchy of A, with the
children of each node sorted according to the B hierarchy preorder.

Where there is a graft of an unmatched subtree of B into the A
hierarchy, only the root of the B subtree is shown, not its
descendants.

## Tools

### NCBI to Darwin Core CSV

This converts the NCBI taxonomy format into something more similar to
what GBIF and EOL use, which I'll call "Darwin Core CSV".

    python3 src/ncbi_to_dwc.py {directory} --out {path}

e.g.

    python3 src/ncbi_to_dwc.py ncbi-dump --out ncbi.csv

You can also generate a TSV file by changing the specified output file name

    python3 src/ncbi_to_dwc.py ncbi-dump --out ncbi.tsv

### Extract a subset checklist

Many of these taxonomy files (NCBI, GBIF, etc) are large and difficult
to work with during development where you want rapid turnaround on
changes to code.  For this reason it's useful to be able to work with
some particular subset, for example mammals, primates, birds, mosses,
etc.

    python3 src/subset_dwc.py {source} {root} --out {path}

Source is a CSV or TSV file containing a larger checklist, and root is
a `taxonID` occurring in the source.  The checklist stored at the
given path is a subset of the source where a row is kept if the root
node is reachable from it via parent links.  If source is not a
taxonomy (a source of parent pointers), a taxonomy over the same
taxonIDs can be provided with `--taxonomy {taxonomy}`.


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

### Convert NCBI taxonomy to Darwin Core form

NCBI has its own taxonomy dump format, which needs to be converted to
DwCA.  There's a tool for this purpose:

    python3 src/ncbi_to_dwca.py work/ncbi/2020-01-01/dump \
      --out work/ncbi/2020-01-01/converted.csv

(The GBIF files are DwCA format already, so they can be used directly.)

### From each, extract Primates only 

Dealing with the whole taxonomies would be time consuming and painful
to analyze, so it's best to just look at reasonably sized taxa.
The `subset_dwca.py` tool extracts just the part we want from the whole
taxonomy's DwCA.

    python3 src/subset_dwca.py work/gbif/2019-09-16/dwca 798 \
      --out work/gbif/2019-09-16/primates.csv
    python3 src/subset_dwca.py work/ncbi/2020-01-01/dwca 9443 \
      --out work/ncbi/2020-01-01/primates.csv

You can find the identifiers for the taxon of interest via each site's
web interface, or using `grep`.

### Compare them

    python3 src/cldiff.py work/gbif/2019-09-16/primates.csv \
      work/ncbi/2020-01-01/primates.csv --out diff.out

The first checklist (or taxonomy) is the "A checklist" and the second is
the "B checklist".

If you put the two in the opposite report you'll get the comparison
ordered by the other taxonomy:

    python3 src/cldiff.py work/gbif/2019-09-16/primates.csv work/ncbi/2020-01-01/primates.csv

## Darwin Core fields of interest:

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
