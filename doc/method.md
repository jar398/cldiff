# Overview

The "checklist diff" tool operates on a pair of checklists
('checklist' is construed so as to include taxonomies), and combines
and compares them in a variety of ways.

The tool has three inputs: two checklists and a set (perhaps empty) of
"curated" relationships between records in the two checklists [not
currently implemented 2020-08-04].  (The curated relationships allow
interventions to augment or override what the automated process would
come up with.)  One checklist is "low priority" and the other is "high
priority".

The method has some features in common with the "smasher" tool from
the Open Tree of Life project [ref], whose merge operation is an
"underlay".

All operations on checklists proceed as follows:

 1. Ingest the two checklists (see section 'Checklists')
 1. Align the two checklists (see 'Alignments')
 1. Merge the two checklists into a single annotated hierarchy (see 'Merge')
 1. Generate output: report, alignment, script, etc.


## Checklists

A checklist is a set of taxon records or "nodes" (I use the two words
more or less interchangeably).  Each node is intended to correspond to
some taxon (in the extensional sense - set of organisms), but multiple
nodes can correspond to the same taxon.  In this case one of the
taxon's nodes is 'accepted' and all the others are 'synonyms'.

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


### Alignment ideal

Given the particle and extension analyses, a set of matches can be
proposed.

There are certain desirable properties of alignments (when considered
in logical conjunction with the two checklists):

 1. An alignment is _consistent_ if at most one RCC-5 articulation can be
    inferred between any pair of nodes (one from each checklist).
 1. An alignment is _complete_ if between any pair of nodes (one in
    each checklist) exactly one RCC-5 articulation can be inferred.
 1. A complete alignment is a "basis" if removing any alignment would
    make it incomplete.

Because some biological information is missing from the checklists,
and the biological meaning of information in the checklists is not
always clear, this ideal is not always possible, but we can try to
approach it as best we can through an automated process, based on the
intensional and hierarchical information in the two checklists.

The generated alignment will never contain a disjointness
articulation, although the curation list might.

### Matching particles

The alignment we seek contains all of the particles (a, =, b), except
where one of these would be inconsistent with the curated
articulations.

Descendants of particles (with exceptions induced by the curated
articulations) will not participate in the alignment.

### Co-extensional node sets

Nodes drawn from one or both checklists are said to be
'co-extensional' if they subtend the same (nonempty) particle set.
As a special case, the nodes of a particle are co-extensional.

Co-extension is computed via a 'cross-MRCA' calculation.  The
cross-MRCA of a node X is the MRCA in the opposite checklist of the
nodes corresponding to X's subtended particles.  Co-extensional nodes
are one another's cross-MRCAs, but the converse doesn't always hold.
However, cases where the converse fails to hold
can be calculated via a quick search, making extension analysis nearly linear
when the two checklists have similar hierarchies.

Co-extension is strong evidence that nodes should be aligned, but
co-extension in terms of particle sets doesn't automatically imply
identity _of the intended taxa_.  The intended taxa may be
distinguishable or identifiable according to members belonging to the
unmatched (non-particle) descendants, or in some other way.  But since
the checklist representation is missing this information, there is
little to be done to resolve such situations automatically.

A node's 'co-extensionality class' is the set of nodes in either
checklist that have the same particle extension.  If the checklists
are similar enough, each such class will contain one node from each
checklist (two total), and we can posit an articulation between the
two nodes.  But in case any equivalence class contains three or more
nodes, care must be taken.

Within a single checklist, co-extensional accepted nodes should be
presumed distinct.  They can be distinguished either by their
non-particle descendants (e.g. tips not matched to the other
checklist) or by their properties (intension) such as canonical name.

Co-extensional nodes in a checklist are related by ancestry; they form
a partial lineage.  A similarly co-extensional node in the other
checklist can potentially be matched to any member of the chain.  An
attempt is made to match the node to a chain node intensionally.  If
this fails, an ambiguity remains that may require manual intervention.


### Splitting and merging

If multiple low priority nodes match a high priority node, we have
'lumping' or 'merging' and the 'extra' low priority nodes can be
safely dropped.

If multiple high priority nodes match a low priority node, we have
an instance of 'splitting'.


### Conflict between intensional and extensional matches

If x has y has its best intensional match, and z (different from y) as
an extensional match, we have an exceptional situation that should
probably get manual review.  While in many cases this means that one
higher priority checklist is a better reflection of evolutionary and
taxonomic truth than the other checklist, it might also reflect an
error introduced into the high priority checklist by a typo,
clerical error, or overzealous match elsewhere in the trees.

[TBD: the tool should detect and report on these situations]

### Non-equivalence articulations

Extension analysis creates opportunities to record articulations other
than equivalences.  The other three RCC-5 relations all express useful
information that cannot be derived from equivalence articulations.

We obtain an articulation (a, >, b) when the high priority taxon b "splits" the 
low priority taxon a.  In the merged checklist, b is the parent of a.  In a
split, b is maximal, so (a, >, b') shouldn't go in the
alignment if b' is a descendant of b.

Also useful is the articulation (a, ><, b) where a is maximal (in the
same manner as above).  This says that the low priority node a cannot
be included in a merged checklist because it is inconsistent with a
hierarchy containing the high priority node.

As mentioned above, disjointness is so pervasive that it does not help
much to express it - it is almost always already implied by automatic
sibling disjointness captured in the hierarchies.

## Merged hierarchy

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
