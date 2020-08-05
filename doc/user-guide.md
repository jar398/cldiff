# User guide to the 'diff' tool

## Input

Command line:

    python3 src/report.py --help
    usage: report.py [-h] [--low-tag LOW_TAG] [--high-tag HIGH_TAG] [--out OUT] [--format FORMAT] low high

    positional arguments:
      low                  lower priority checklist
      high                 higher priority checklist

    optional arguments:
      -h, --help           show this help message and exit
      --low-tag LOW_TAG
      --high-tag HIGH_TAG
      --out OUT            file name for report
      --format FORMAT      report format

The two checklists are given in files with either TSV (tab separated)
or CSV (comma separated) format.  The file names should end in .tsv or
.csv to signal the format.

The first row of each checklist should give column headings.  Certain
headings (mostly Darwin Core) are known to the program:

    taxonID              identifies record in accepted or parent link
    canonicalName        e.g. Plagiaulax minor
    scientificName       e.g. Plagiaulax minor Falconer, 1857
    parentNameUsageID    parent in hierachy
    acceptedNameUsageID  the record for the taxon of which this record gives a synonym
    taxonRank
    nomenclaturalStatus
    taxonomicStatus


### Extracting a subset of a checklist

A checklist such as the GBIF backbone can be quite large and it is
useful to be able to work with smaller parts of it.  There is a
utility for extracting a single clade from a checklist.

[Documentation to be written!]

    python3 src/subset_dwc.py --help
    usage: subset_dwc.py [-h] [--taxonomy TAXONOMY] [--out OUT] source id

    positional arguments:
      source               taxonomy or checklist from which to extract subset
      id                   taxon id of subset's root

    optional arguments:
      -h, --help           show this help message and exit
      --taxonomy TAXONOMY  taxonomy from which to extract hierarchy; defaults to source
      --out OUT            where to store the subset

### Converting an NCBI Taxonomy dump to CSV

[Documentation to be written!]

    python3 src/ncbi_to_dwc.py --help
    usage: ncbi_to_dwc.py [-h] [--out OUT] dump

    positional arguments:
      dump        directory containing taxdump files

    optional arguments:
      -h, --help  show this help message and exit
      --out OUT   where to store the DwC version

## Output

The 'diff' operation generates a single CSV file that is intended to
be human readable.  Each row gives information about a node (taxon) in the merged taxonomy.

I continue to experiment with the columns quite a bit but as of 5
August the diff report has these columns:

 1. `indent` - shows the hierarchical structure of the merged taxonomy;
    this can be interpreted as specifying the parent pointer of each node
    (the parent is the closest row going upwards whose indentation is less than that of the current row)
 1. `operation` - `KEEP` (taxon in both low and high priority checklists),
    `DELETE` (low only), `ADD` (high only)
 1. The next five columns show either an articulation, or a single low or high priority taxon.
     * `dom` - domain - the left hand node of the articulation (low priority)
     * `dom id` - the taxonID of the domain
     * `relation` - an RCC-5 relation, defaulting to `=`
     * `cod id` - the taxonID of the codomain
     * `cod` - codomain - the right hand node of the articulat (high priority)
    The taxon in the merged taxonomy is indicated by the `dom` for `DELETE` rows, by the `cod` 
    for `ADD` rows, and by either the domain or codomain for the `KEEP` rows (they are 
    judged equivalent).  In either case, the other taxon is a sample with the 
    articulation provided just for documentation.
 1. `unchanged` - if the value in this column is `subtree=` then all
    descendants have been suppressed because the subtrees of the two
    checklists at the point are identical
 1. `changed_props` - a list of properties (record fields) that differ between the two checklists
 1. `reason` - the justification for this match.
     * `extension` means they have the same subtended particles
     * `name` means they have the same name
     * `synonym+name` means the domain has a synonym with the same name as the codomain [I might have this backwards]
     * `name+synonym` means the codomain has a synonym with the same name as the domain
     * `synonym+name+synonym` means the domain and codomain have synonyms with the same name


### Commentary

As the tool runs it generates a running commentary.  The most useful
of these comments are the explanations of how taxa were determined to
be inconsistent.  Example:

    # B.Hominoidea conflicts with C.Primates because:
    #   child C.Galagonidae is disjoint while C.Hylobatidae isn't
    # B.Homininae conflicts with C.Hominidae because:
    #   child C.Pongo#5219531 is disjoint while C.Pan isn't

Here B is low priority and C is high, and these comments explain why
B.Hominoidea and B.Homininae are determined to be excluded from the
merged hierarchy (`DELETE`).  For example, low priority Homininae is
excluded because the high priority checklist has a taxon Hominidae
with children Pongo and Pan, and Pongo is contained in Hominoidea
while Pan is not (by analysis of extensions based on 'particle'
matches).

