# Overview

The "checklist diff" tool operates on a pair of checklists
('checklist' is construed so as to include taxonomies), and combines
and compares them in a variety of ways.

The tool has three inputs: two checklists and a set (perhaps empty) of
"curated" relationships between records in the two checklists [not
currently implemented 2020-08-04].  (The curated relationships allow
interventions to augment or override what the automated process would
come up with.)  One checklist is "low priority" and the other is "high
priority".  The following operations are supported:

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

[2020-08-04 Only the first of these is implemented; the second is
partially implemented]

The method has some features in common with the "smasher" tool from
the Open Tree of Life project [ref], whose merge operation is an
"underlay".

All operations proceed as follows:

 1. Ingest the two checklists (see section 'Checklists')
 1. Align the two checklists (see 'Alignments')
 1. Merge the two checklists into a single annotated hierarchy (see 'Merge')
 1. Generate output: report, script, etc.


## Checklists

A checklist is a set of taxon records or "nodes" (I use the two words
more or less interchangeably).  Each node is intended to correspond to some
taxon (in the extensional sense - set of organisms), but multiple
nodes can correspond to the same taxon (synonyms).

A taxonomy connotes a checklist where hierarchy information
(parent/child) is important.  I'll use the words "checklist" and
"taxonomy" more or less interchangeably.

A node consists of a set of Darwin Core fields appropriate to
`dwc:Taxon` records, and perhaps others.  It contains

 * optional node identifier (`dwc:taxonID`)
 * one or more names (canonical name, scientific
   name, optionally scientific name split up into genus, species epithet,
   etc.)
 * other name- or taxon-related properties (rank, nomenclatural status, etc)
 * an optional 'parent' link which if present gives a node for a
   taxon that strictly contains taxon for this node (position in hierarchy)
   (`parentNameUsageID` when serialized)
 * an optional 'accepted' link which if present signals that this
   node is for a synonym (`acceptedNameUsageID` when serialized)

Links can be interpreted in either direction, so that each node also
has, by implication, a set of children (nodes in the same checklist
that have the node as their parent) and a set of synonyms.

I'll call a node that has no 'accepted' link (i.e. one that is not a
synonym) an 'accepted' node.


## Alignments

All operations begin with the determination of an alignment between
two checklists/taxonomies, say A (low priority) and B (high priority).  An
alignment is a set of 'articulations' each relating a node of A to a node of
B.  An articulation is a triple (a, r, b) where a is a node in A, b is
a node in B, and r is one of the five RCC-5 relations {=, <, >, !, <>}
(same extension, strictly contained in, strictly contains, disjoint,
overlap without containment).

Articulations are advanced based on names, synonymies, and/or
membership (for higher taxa).  Alignment generation is 'eager' in that
it can easily create articulations that do not hold when information
not captured in the checklists is brought to bear.  For this reason it
is expected that some automatic articulations will be overridden by
providing curated alignments that override them.

To be able to produce alignments as textual output (as well as curated
articulations), there has to be a way to indicate particular nodes
(records) in checklists.  For now, assume this is via the identifier in the node's
`taxonID` column, accompanied, when not clear from context, by a tag
or name saying which of the checklists the node is taken from.

An alignment is determined in two stages: 'particle' analysis and
'extension' analysis.


### Particle analysis

In the first stage, a set 'particles' is chosen.  ('Particle' is a
silly word, but I couldn't come up with anything better.)
Each particle is an `=`
articulation (a, =, b) between two accepted nodes, one in each checklist.
Particles correspond to hypothetical taxa for which nodes plausibly exist in both
checklists.  Particles cover only the taxa that
are 'terminal' or 'most tipward' in both checklists' hierarchies, and provide a starting point for
a membership-based ('extensional') comparison of the two checklists.

A particle has the following properties.
  
     1. No particle in A or B is the descendant of another particle.
     1. A particle's nodes are a 'best mutual intensional match' i.e.
         1. nodes a and b 'match' intensionally 
         1. a and b are one another's best intensional match, mutually, under a quality
            ordering.

The nodes of a particle match 'intensionally', according to their
properties (most importantly their names) and the properties of their
synonyms.

The quality ordering is a bit in flux, but roughly involves

     1. relying on the fewest possible synonyms
     1. having node fields in common; the fields have a
        priority order (e.g. matching on 'scientific name' is
        better than matching on 'canonical name', and more generic 
        properties such as rank can be used as tie-breakers)


### Extension (membership) analysis

The checklists induce hierarchies whose tips are the members of the
particle set.  These two hierarchies may be different, but we can
compare them by comparing the sets of particles subtended by interior
nodes in the checklists.  Two nodes that subtend the same particle
set might plausibly correspond to the same taxon; or more generally,
an RCC-5 relation between particle sets is a good hypothesis for an
RCC-5 relation between the taxa whose nodes subtend those particle
sets.

In particular, extension analysis can find those situations where a
node in A must be discarded because it conflicts with B's hierarchy.
This usually happens because of an update in classification, where B
is more correct taxonomically than A, but it can also come from errors
in either A or B, or errors in the alignment.

In the case of a 'tie' among extensional matches, 


## Alignment

Given the particle and extension analyses, an alignment can be
proposed.

There are certain desirable properties of alignments (when considered
in logical conjunction with the two checklists):

 1. An alignment is _consistent_ if at most one RCC-5 articulation can be
    inferred between any pair of nodes (one from each checklist).
 1. An alignment is _complete_ if between any pair of nodes (one in
    each checklist) exactly one RCC-5 articulation can be inferred.
 1. A complete alignment is a "basis" if removing any alignment would
    make it incomplete.

Because some information is missing, and the biological meaning of
information in the checklists is not always clear, this ideal is not
always possible, but we can try to approach it as best we can through
an automated process, based on the intensional and hierarchical
information in the two checklists.

The generated alignment will never contain a disjointness
articulation, although the curation list might.

### Aligning particles

The alignment we seek contains all of the particles (a, =, b), except
where one of these would be inconsistent with the curated
articulations.

Descendants of particles (with exceptions induced by the curated
articulations) will not participate in the alignment.

### Co-extensional node sets

Nodes drawn from one or both checklists are said to be
'co-extensional' if they subtend the same (nonempty) particle set.
As a special case, the nodes of a particle are co-extensional.

Co-extension is not computed exactly but rather by using a
'cross-mrca' calculation.  The cross-mrca of a node is the mrca in the
opposite checklist of the node's subtended particles.  Co-extensional
nodes are one another's cross-mrcas; cases where the converse fails to
hold can be calculated via a quick search, making extension analysis
linear when the two checklists have similar hierarchies.

Co-extension is strong evidence that nodes should be aligned, but
co-extension in terms of particle sets doesn't automatically imply
identity of the intended taxa.  The taxa may be distinguishable
according to members belonging to the unmatched (non-particle)
descendants, or descendants that are simply not recorded in a
checklist at all.  But since the checklist representation is missing
this information, there is little to be done to account for it.

CHECK THIS

Nodes form equivalence classes according to co-extensionality (one
equivalence class per distinct particle set).  If the checklists are
similar enough, each equivalence class will contain one node from each
checklist, and we can posit an articulation between the two nodes.
But in case any equivalence class contains three or more nodes, care
must be taken.

Within a single checklist, co-extensional non-synonym nodes should be
presumed distinct.  They can be distinguished either by their other
children (those not in the co-extensional set), which do not share
particles with the other checklist, or simply by intension
(properties).

For an equivalence class C, ideally we would like those nodes of A
that are in C to match those nodes of B that are in C.  This can be
often done via intension.

Some of the nodes in chains may not have any match by intension, or
may not have an unambiguous match.  If a low priority node in a chain
has only one child, and doesn't match a high priority node
unambiguously, it can be left unaligned (it appears to be a synonym,
and perhaps is too troublesome to warrant a merge).  But in general
these situations probably require manual intervention.

### Splitting and merging

DELETE ME

A node x may have an intensional best match, say y, that is not y's
intensional best match.  There is then an opportunity for lumping or
splitting, similar to that for co-extensional node sets.

This situation is quite unclear.  A low-priority node, say x, could
have high priority y as its best match, while high priority z
(different from y) could have x as its best match.  It is hard to come
up with such an example, but it should probably be detected and
flagged as a curation demand.

### Synonyms

DELETE ME

A merged checklist should contain all synonyms from both sources.
If b is the best match of a, then a and a's synonyms should all be
represented by synonyms of b in the merged checklist.  Therefore,
synonyms in b need to be matched to nodes in a, perhaps other
synonyms, with articulations added to the alignment as needed.

An identity articulation ("same node in both checklists") needs to be
distinguished from a synonym articulation ("a in A is a synonym of b
in B"); the basis for such a judgment is not completely clear, but it
can probably be done by simply comparing canonical names, assuming
that any changes in authority are corrections, and that there is no
harm in considering corrected canonical name misspellings to synonyms
new to checklist A.

### Conflict between intensional and extensional matches

If x has y has its best intensional match, and z different from y as
an extensional match, we have an exceptional situation that should
probably get manual review.  While in many cases this means that one
higher priority checklist is a better reflection of evolutionary and
taxonomic truth than the other checklist, it might also reflect an
error introduced into in the high priority checklist by a typo,
clerical error, or overzealous match elsewhere in the trees.

[TBD: the tool should detect and report on these situations]

### Non-equivalence articulations

Extension analysis creates opportunities to record articulations other
than equivalences.

As mentioned above, disjointness is so pervasive that it does not help
much to express it - it is almost alway already implied by sibling
disjointness in the hierarchies.  But the other three RCC-5 relations
all express useful information that cannot be derived from equivalence
articulations.

We obtain an articulation (a, >, b) when the high priority taxon b "splits" the 
low priority taxon a.  In the merged checklist, b is the parent of a.  In a
split, b is maximal, so (a, >, b') shouldn't go in the
alignment if b' is a descendant of b.

Also useful is the articulation (a, ><, b) where a is maximal (in the
same manner as above).  This says that the low priority node a cannot
be included in a merged checklist because it is inconsistent with a
hierarchy containing the high priority node.

## Merge

With an alignment in hand, it is possible to internally construct a
merged hierarchy.  Each node in the merged hierarchy is either 'kept'
(mutual best match including all particles), 'low only' or 'deleted' (present in low
priority checklist but not in high), or 'high only' or 'added' (in high but not in
low).  Parent/child relationships can hold between any of these
kinds of merged node (other than between a low only and a high only).

The central task of merging is determining what the parent ought to be.
Each 'kept' node may have two candidate parents, namely the parents of
the component nodes in each checklist.  The 'smaller' of the two
parents is chosen according to the extension analysis.

Once parents are determined, the child-parent relationship can be
inverted to obtain the set of children of each node in the merged hierarchy.


## Outputs

All outputs are driven off the merged hierarchy.

[2020-08-04 TBD: only the diff report is ready right now.]

### Diff

The output is an abridged version of the merged hierarchy in a human readable CSV form.

Subtrees that are identical between the low and high priority
checklists are elided in order to focus attention where there is
taxonomic action.  This can result in a substantially smaller file.


### Alignment

The output is simply the articulations of the alignment, with reversals removed (there is no reason to have both x = y and y = x).

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
