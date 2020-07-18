# Overview

The "checklist diff" (or more appropriately "checklist patch") tool
operates on checklists and taxonomies, and generates sequences of
commands to merge some or all of the information from one of the
checklists into the other checklist, producing a third checklist
(which in some circumstances is the same as one of the inputs).

Following is a description of checklists followed by a description of
the merge operations.  After that is a closer examination of
alignments, which are an essential artifact in this process.


## Checklists

A checklist is a list of taxon records.  Each taxon record is intended
to correspond to some taxon (in the extensional sense - set of
organisms), but multiple records can correspond to the same taxon
(synonyms).

A taxonomy is a checklist that happens to include hierarchy
(parent/child) information.  I will use the words "checklist" and
"taxonomy" more or less interchangeably.

A taxon record follows Darwin Core.  It contains 
 * one or more names (canonical name, scientific
   name, perhaps scientific name split up into genus, species epithet,
   etc.)
 * an optional 'accepted' link which if present signals that this
   record is for a synonym
 * an optional 'parent' link which if present gives a record for a
   taxon that strictly contains the one for this record
 * perhaps other properties (nomenclatural status, etc)


## Varieties of merge

All merge operations take a "low priority" checklist A and a "high
priority" checklist B, and produce a merged checklist C.  Whenever
there are conflicts between A and B, e.g. in hierarchy (a record can
only have one parent) or properties, the relationships in C are made
to be compatible with B.

The generated commands can be one of
  - commands to transform A into C, or
  - commands to transform B into C.

We can [or rather will be able to] do one of the following:
  - overlay A on B, with hierarchical relationships in A overriding those in B.
  - replace subtrees of A with corresponding corresponding subtrees from B
    (in particular: when A and B have the same root, C = B)

In each case an optional input is a 'seed' set of articulations (see
below) that override mistakes or omissions made by the automatic
articulation generator (see below).

Use cases:
  - Simple comparison of A and B: doesn't much matter which direction 
    the commands go; a matter of user preference.  (overlay)
  - Add new taxonomic information to EOL dynamic hierarchy in database, 
    without losing existing information:
    A = version in database (all life), B = new information (perhaps 
    just a small group), C = mostly A but with B added, overriding 
    conflicting A-information if there are conflicts
  - Updating an EOL dynamic hierarchy in place: A = the version in the
    database, B = newly synthesized version, get commands to alter 
    database replacing A with B, C = B  (replace)

Generating commands from an alignment

  - Hmm


## Alignments

All operations begin with the determination of an alignment between
the two checklists (or taxonomies) A and B.  An alignment is a set of
'articulations' relating records of A to records of B.  An
articulation is a triple (a, r, b) where a is a record in A, b is a
record in B, and r is one of the five RCC5 relations {=, <, >, !, <>}
(same extension, strictly contained in, strictly contains, disjoint,
overlap without containment).

To be able to write down alignments, there has to be a way to identify
records.  For now assume this is via the 'taxonID' column and, when
not clear from context, some designator saying which checklist the
record is taken from.

An alignment is determined in two stages: 'particle' (or 'fringe')
analysis, and extension analysis.

### Particle analysis

In the first, a set 'particles' is chosen.  Each particle is an `=`
articulation (a, =, b) between two records, one in each checklist.
Particles correspond to taxa for which records plausibly exist in both
checklists.  Particles cover not all such taxa, but only the ones that
are 'terminal' or 'most tipward', and provide a starting point for
a membership based comparison of the two checklists.

('Particle' is a kind of silly word but I haven't come up with
anything better.  'Tip' and 'shared tip' work much of the time
but there are cases where one record or the other is not a tip in the
checklist's hierarchy, so it would be confusing to use that word.)

A particle has the following properties.
  
     1. No particle in A or B is the descendant of another particle.
     1. Every record in A and B has a particle as either a descendant or an ancestor.
     1. The records are a 'best mutual name match' i.e.
         1. a and b 'match' i.e. they have a name in common, perhaps after 
            taking synonyms into consideration
         1. a and b are the best such match, mutually, under a quality ordering

The quality ordering is a bit in flux, but roughly involves

     1. relying on the fewest possible synonyms
     2. having as many record fields in common as possible (with 
        a field priority order to break ties)

### Extension (membership) analysis

Each checklist induces a hierarchy whose tips are the members of the
particle set.  These two hierarchies may be different, but we can
compare them by comparing the sets of particles subtended by interior
records in the checklists.  Two records that subtend the same particle
set might plausibly correspond to the same taxon; or more generally,
an RCC5 relation between particle sets is a good hypothesis for an
RCC5 relation between the taxa whose records subtend those particle
sets.

In particular, extension analysis can find those situations where a
record in A must be discarded because it conflicts with B's hierarchy.
This usually happens because of an update in classification, where B
is more correct taxonomically than A, but it can also come from errors
in either A or B, or errors in the alignment.


## Alignment formulation

Given the particle and extension analysis, an alignment can be
proposed.

Clearly the alignment contains all of the particles (a, =, b).

Descendants of particles will not participate in the alignment.

For every ancestor of a particle a in A, we can find a 'most
appropriate' extensional match in B (perhaps selecting between
extensional matches using name matching).

TBD

The generated alignment will never contain a disjointness
articulation, although the exceptions list might.

Just because we have co-extension in terms of particle sets doesn't
mean we have co-extension of the intended taxa.


### Co-extensional sets

A difficulty arises with 'chains', i.e. maximimal parent/child
sequences for which the subtended particle set is the same for all
records (equivalently: each parent only has one child that subtends a
particle).  For a given chain there may be a chain in the other
checklist with the same subtended particle set.  Together these chains
form I'll call a 'co-extensional' set of records.

If a match is unique (an A-record has only one particle-set match in
B, even though there is no name match), there is no ambiguity and an
articulation based on particle set can be proposed.

In this case, we first revert to name-based matching (see above) to
sort out the potentially many-to-many correspondence between all the
records subtending the same particle set.

Some of the records in chains may not have a match by name, or not
have an unambiguous match.  If an A-record in a chain has only one child, and doesn't
match a B-record unambiguously, it can be dropped.  However, if one or
more unmatched A-records in a chain.

### Splitting and merging

This requires some explicit thought.  TBD
