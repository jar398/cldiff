# User guide to the 'diff' tool

The following modes of operation are supported:
[2020-08-04 Only the first of these is implemented; the second is
partially implemented; the rest should be straightforward given 
the existing infrastucture]

 1. Report on how the two checklists compare ("diff")
 1. Provide an alignment in Euler/X format for use with reasoning tools
 1. "Underlay mode":
    Generate commands to modify the high-priority checklist
    to add information provided by the low-priority checklist
    (high priority takes precedence in case of conflict)
 1. "Overlay mode":
    Generate commands to modify the low-priority checklist
    to add information provided by the high-priority checklist
    (high priority takes precedence in case of conflict)
 1. "Update mode":
    Generate commands to transform an instance of the low-priority checklist
    into the high-priority checklist, removing information that's 
    only in the low-priority checklist

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

Other columns may be present, but they are ignored by the program.


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

Source is a CSV or TSV file containing a larger checklist, and root is
a `taxonID` for a record in the source.  The checklist stored at the
given path is a subset of the source where a row is kept if the root
node is reachable from it via parent links.  If source is not a
taxonomy (a source of parent pointers), a taxonomy over the same
taxonIDs can be provided with `--taxonomy {taxonomy}`.

### Converting an NCBI Taxonomy dump to CSV

NCBI has its own taxonomy dump format, which needs to be converted to
DwCA.  There's a tool for this purpose:

    python3 src/ncbi_to_dwc.py --help
    usage: ncbi_to_dwc.py [-h] [--out OUT] dump

    positional arguments:
      dump        directory containing taxdump files

    optional arguments:
      -h, --help  show this help message and exit
      --out OUT   where to store the DwC version

E.g.

    python3 src/ncbi_to_dwca.py work/ncbi/2020-01-01/dump \
      --out work/ncbi/2020-01-01/converted.csv

(The GBIF files are DwCA format already, so they can be used directly.)

## Output

The 'diff' operation generates a single CSV file that is intended to
be human readable.  Each row gives information about a node (taxon) in the merged taxonomy.

I continue to experiment with the columns quite a bit but as of 2020-08-05
the diff report has these columns:

 1. `indent` - shows the hierarchical structure of the merged taxonomy;
    this can be interpreted as specifying the parent pointer of each node
    (the parent is the closest row going upwards whose indentation is less than that of the current row)
 1. `operation` - `KEEP` (taxon in both low and high priority checklists),
    `DELETE` (low only), `ADD` (high only)
 1. The next five columns show either an articulation, or a single low or high priority node.
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
    checklists at this point are identical
 1. `changed_props` - a list of properties (record fields) that differ 
    between the two records (`KEEP` only)
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
merged hierarchy (`DELETE`).  For example, low priority _Homininae_ is
excluded because the high priority checklist has a taxon _Hominidae_
with children _Pongo_ and _Pan_, and _Pongo_ is contained in _Hominoidea_
while _Pan_ is not (by analysis of extensions based on 'particle'
matches).



## Output modes

All outputs are driven off the merged hierarchy.

[2020-08-05 TBD: only the diff report is ready right now.]

### Diff

The output is an abridged version of the merged hierarchy in a human readable CSV form.  See above.

Subtrees that are identical between the low and high priority
checklists are elided in order to focus attention where there is
taxonomic action.  This can result in a file that is substantially smaller than one covering the entire merged hierarchy would be.


### Alignment

The output is simply the articulations of the alignment, with reversals left out (there is no reason to have both x = y and y = x).

[2020-08-04 TBD: Also need to render the two checklists in Euler/X taxonomy format.]

### 'Underlay' - adding low-priority nodes to high-priority checklist

The command set consists of

  1. Node additions for the 'low only' nodes in the merged hierarchy,
     excluding those that are inconsistent (`><`) with the high priority checklist
  1. Node field additions ?  for situations where a low only node has
     property values that are missing from the high only node that is
     matched?

### 'Overlay'

The command set consists of

  1. Node additions for the 'high only' nodes
  1. Node field changes to modify matched low nodes so that they have information from high nodes
  1. Removal of inconsistent 'low only' nodes, taking care to update parent pointers

### Update

  1. Node additions for the 'high only' nodes
  1. Node field changes to modify matched low nodes so that they have information from high nodes
  1. Removal of all 'low only' nodes
