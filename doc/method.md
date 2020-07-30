# Overview

The "checklist diff" (or more appropriately "checklist patch") tool
operates on checklists and taxonomies, and generates sequences of
commands to merge some or all of the information from one of the
checklists into the other checklist.

The tool has three inputs: two checklists and a set (perhaps empty) of
"curated" relationships between records in the two checklists.  The
curated relationships allow interventions to augment or override what
the automated process would come up with.

The method has some features in common with the 'smasher' tool from
the Open Tree project [ref].

Following is a description of checklists in general, followed by a
description of the merge operations.  After that is a closer
examination of alignments and articulations, which are essential
players in the operation of the tool.


## Checklists

A checklist is a set of taxon records or "nodes" (I use the two words
interchangeably).  Each node is intended to correspond to some
taxon (in the extensional sense - set of organisms), but multiple
nodes can correspond to the same taxon (synonyms).

A taxonomy connotes a checklist where hierarchy (parent/child)
information is important.  I'll use the words "checklist" and
"taxonomy" more or less interchangeably.

A node consists of a set of Darwin Core fields appropriate to
dwc:Taxon records, and perhaps others.  It contains

 * one or more names (canonical name, scientific
   name, optionally scientific name split up into genus, species epithet,
   etc.)
 * other name- or taxon-related properties (rank, nomenclatural status, etc)
 * an optional 'accepted' link which if present signals that this
   node is for a synonym
 * an optional 'parent' link which if present gives a node for a
   taxon that strictly contains taxon for this node (position in hierarchy)

Links can be interpreted in either direction, so that each node also
has, by implication, a set of children (those that have the node as
their parent) and a set of synonyms.

I'll call a node that has no 'accepted' link (i.e. one that is not a
synonym) an 'accepted' node.

I'll call this information, other than the parent and child
information, 'intentional' information.


## Alignments

All operations begin with the determination of an alignment between
the two checklists (or taxonomies) A and B.  An alignment is a set of
'articulations' relating nodes of A to nodes of B.  An
articulation is a triple (a, r, b) where a is a node in A, b is a
node in B, and r is one of the five RCC5 relations {=, <, >, !, <>}
(same extension, strictly contained in, strictly contains, disjoint,
overlap without containment).

Articulations are proposed based on names, synonymies, and/or
membership (for higher taxa).  Alignment generation is 'eager' in that
it can easily create articulations that do not hold when information
not captured in the checklists is brought to bear.  For this reason it
is expected that some automatic articulations will be overridden by
providing curated alignments that override them.

To be able to write down alignments (as well as curated
articulations), there has to be a way to designate nodes (records) in
the two checklists.  For now assume this is via the node's 'taxonID'
column, accompanied, when not clear from context, by a tag or name
saying which of the two checklists the node is taken from.

An alignment is determined in two stages: 'particle' (or 'fringe')
analysis and extension analysis.


### Particle analysis

In the first, a set 'particles' is chosen.  Each particle is an `=`
articulation (a, =, b) between two accepted nodes, one in each checklist.
Particles correspond to taxa for which nodes plausibly exist in both
checklists.  Particles cover not all such taxa, but only the ones that
are 'terminal' or 'most tipward', and provide a starting point for
a membership based comparison of the two checklists.

('Particle' is a completely made-up usage here and must be taken only
as defined for present purposes.  It is a kind of silly word but I
haven't come up with anything better.  'Tip' and 'shared tip' would
work much of the time but there are cases where one node of a particle
is not a tip in its checklist, so it would be confusing to say 'tip'.)

A particle has the following properties.
  
     1. No particle in A or B is the descendant of another particle.
     1. Every node in A and B has a particle as either a descendant or an ancestor.
     1. A particle's nodes are a 'best mutual intensional match' i.e.
         1. a and b 'match' intensionally i.e. they have a name and/or
            identifier in common, or there is such a match involving a
            synonym of one or synonyms of both
         1. a and b are the best intensional match, mutually, under a quality
            ordering.  Matches requiring synonyms are not as
            good as those that don't, and matches a greater number of node
            properties are preferred.  By 'mutual' is meant that a
            must be b's best match and b must be a's best match.

The quality ordering is a bit in flux, but roughly involves

     1. relying on the fewest possible synonyms
     2. having as many node fields in common as possible (with 
        a field priority order to break ties)

A particle may be declared as one of the curated articulations.


### Extension (membership) analysis

Each checklist induces a hierarchy whose tips are the members of the
particle set.  These two hierarchies may be different, but we can
compare them by comparing the sets of particles subtended by interior
nodes in the checklists.  Two nodes that subtend the same particle
set might plausibly correspond to the same taxon; or more generally,
an RCC5 relation between particle sets is a good hypothesis for an
RCC5 relation between the taxa whose nodes subtend those particle
sets.

In particular, extension analysis can find those situations where a
node in A must be discarded because it conflicts with B's hierarchy.
This usually happens because of an update in classification, where B
is more correct taxonomically than A, but it can also come from errors
in either A or B, or errors in the alignment.


## Alignment formulation

Given the particle and extension analyses, an alignment can be
proposed.

There are certain desirable properties of alignments (when considered
in logical conjunction with the two checklists):

 1. An alignment is _consistent_ if at most one RCC5 articulation can be
    inferred between any pair of nodes (one from each checklist).
 1. An alignment is _complete_ if between any pair of nodes (one in
    each checklist) exactly one RCC5 articulation can be inferred.
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

Co-extension is strong evidence that nodes should be aligned, but
co-extension in terms of particle sets doesn't automatically imply
identity of the intended taxa.  The taxa may be distinguishable
according to members belonging to the unmatched (non-particle)
descendants, or descendants that are simply not recorded in a
checklist at all.

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
often done via intension.  When it can't, and an ambiguity remains, we
have choices:

  1. lumping (when there are multiple low priority nodes matching a
     single high priority node)
  1. splitting (when there are multiple high priority nodes matching
     a single low priority node)
  1. refuse to align, and demand that matches be done as curated inputs

Some of the nodes in chains may not have any match by intension, or
may not have an unambiguous match.  If a low priority node in a chain
has only one child, and doesn't match a high priority node
unambiguously, it can be left unaligned (it appears to be a synonym,
and perhaps is too troublesome to warrant a merge).  But in general
these situations probably require manual intervention.

### Splitting and merging

A node x may have an intensional best match, say y, that is not y's
intentional best match.  There is then an opportunity for lumping or
splitting, similar to that for co-extensional node sets.

This situation is quite unclear.  A low-priority node, say x, could
have high priority y as its best match, while high priority z
(different from y) could have x as its best match.  It is hard to come
up with such an example, but it should probably be detected and
flagged as a curation demand.

### Synonyms

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


### Conflict between intentional and extensional matches

If x has y has its best intensional match, and z different from y as
an extensional match, we have an exceptional situation that should
probably get manual review.  While in many cases this means that one
higher priority checklist is a better reflection of evolutionary and
taxonomic truth than the other checklist, it might also reflect an
error introduced into in the high priority checklist by a typo,
clerical error, or overzealous match elsewhere in the trees.

### Non-equivalence articulations

Extension analysis creates opportunities to record articulations other
than equivalences.

As mentioned above, disjointness is so pervasive that it does not help
much to express it - it is almost alway already implied by sibling
disjointness in the hierarchies.  But the other three RCC-5 relations
all express useful information that cannot be derived from equivalence
articulations.

We obtain an articulation (a, <, b) when the a taxon "splits" the b
taxon: in a merged checklist, b ought to be the parent of a.  In a
split, the a taxon is maximal, so (a', <, b) shouldn't go in the
alignment if a' is an ancestor or a descendant of a.

Also useful is the articulation (a, <>, b) where a is maximal (in the
same manner as above).  This says that the low priority node a cannot
be included in a merged checklist because it is incompatible with a
high priority node.

### Adding low-priority nodes to high-priority checklist

Low priority nodes that do not subtend particles represent parts of
the lower priority hierarchy (A) that, on a merge, would need to be
added to the high priority checklist (B).  A node a whose parent node
subtends a particle should be 'attached' or 'grafted' to B.  But in
some cases it is not clear where a should be attached in B.

At the very least, it seems reasonable that the following property
should hold:

 * a's attachment point should be a B-node whose children subtend all
   particles subtended by the siblings (in A) of the node to be attached.

This B-node would be some node that is co-extensional with the MRCA in
B of the sibling's particles:

 - a node co-extensional with the MRCA that matches a's parent, or
 - the MRCA itself, if that MRCA doesn't match any of a's siblings, or 
 - the parent of the MRCA, otherwise.





## Varieties of merge

All merge operations take a "low priority" checklist A and a "high
priority" checklist B, and produce a merged checklist C.  Whenever
there are conflicts between A and B, e.g. in hierarchy (a node can
only have one parent) or properties, the relationships in C are made
to be compatible with those in B.

The generated commands can be either

  - commands to transform A into C 
    (by adding overriding information from B), or
  - commands to transform B into C 
    (by adding consistent information from A)

An additional execution mode is 

  - the command set simply
    transforms A into B, discarding information in A that is not also in B.

Use cases:

  - Simple comparison of A and B: doesn't much matter which direction 
    the commands go; that would be a matter of user preference.  (overlay)
  - 'Smasher-like' use case, high to low priority source taxonomy processing:
    Add new taxonomic information in A to B
    without losing existing information:
    B = version in database (all life), A = new information (perhaps 
    just a small group), C = mostly B but with A added when consistent
  - Update:
    Add new taxonomic information in B to A, overriding 
    conflicting information in B if A conflicts
  - Updating a dynamic hierarchy in place: A = the version in the
    database, B = newly synthesized offline version, commands to alter 
    the database so the database ends up containing B, not A

## Generating commands from an alignment

(a, =, b)  becomes a command to 'keep' node a.  Should be accompanied
           by any needed field changes (name, etc).
(a, <, b)  becomes a command to 'graft'

