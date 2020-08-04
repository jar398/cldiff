debug = False

import table
import property
import checklist as cl

# Difference report.
# Each bit position is 1 for a mismatch, 0 for a match.
# Bit positions correspond to the URIs listed in property.py.
# The bit mask is a 'badness' measure and can be used to sort
# options for matching.

no_diffs = (0, 0)
all_diffs = (-1, -1)

# ---------- Record comparison.

# Degrees of freedom:
#  - columns to consider: only shared, or all of them?
#  - how many columns to report: only most specific, or all of them?
#  - provide commands to
#       add info from A to B, overriding B?
#       add info from A to B, only when consistent with B?
#  - should diff report be threaded?

# TBD: filter out taxonID if idspaces are different

def differences(uid1, uid2):  # mask
  drop = 0
  add = 0
  for prop in table.get_table(uid2).properties:
    if prop:
      v1 = cl.get_value(uid1, prop)
      v2 = cl.get_value(uid2, prop)

      if v1 != v2:
        if v1 != None:
          drop |= 1 << prop.specificity
        if v2 != None:
          add |= 1 << prop.specificity
  if len(cl.get_children(uid1)) != len(cl.get_children(uid2)):
    drop |= 1 << number_of_children.specificity
    add |= 1 << number_of_children.specificity
  return (drop, add)

number_of_children = property.by_name("number_of_children")
number_of_synonyms = property.by_name("number_of_synonyms")

def same(comparison):
  return comparison == (0,0)

# Returns a list of properties ... ?

def unpack(comparison):
  if same(comparison):
    return []
  else:
    props = []
    for prop in property.properties_by_specificity:
      spec = prop.specificity
      (drop, add) = comparison
      if ((drop | add) & (1 << spec)) != 0:
        props.append(prop)
    return props
