# How it works

The tool is essentially creating an alignment, or a set of RCC5
articulations, that are sufficient to specify a transform from tree A
into tree B.

Definition: a 'match' is an articulation bridging different checklists.
 ... ... definitions ...

The method works in two stages.  

* In the first, a set 'particles' is chosen.  Each particle is an `=`
  articulation between a node in each checklist, as similar as
  possible to one another (mainly: same name), with the property
  that no ancestor of a particle is a particle - there are particles
  only near the tips, not in the interior.  Then the set of
  particles is treated as individuals in RCC5 'regions' that can be
  related in the 5 different usual ways.

* Nodes are then related according to their 'regions'.  If the sets are
  the same, the nodes denote equivalent groups and are candidates for
  an `=` articulation assignment.

then, details
