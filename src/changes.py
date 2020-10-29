debug = False

import table
import property
import checklist as cl

# Difference report, comparing two nodes.
# The result is a "comparison" which is a triple (drop, change, add)
# where each element of the triple is an integer mask.

# Each bit position is 1 for a mismatch, 0 for a match.
# Bit positions correspond to the URIs listed in property.py.
# The bit mask is a 'badness' measure and can be used to sort
# options for matching.

no_diffs = (0, 0, 0)
all_diffs = (0, -1, 0)

taxonID = property.by_name("taxonID")
parentNameUsageID = property.by_name("parentNameUsageID")

# ---------- Record comparison.

# Degrees of freedom:
#  - columns to consider: only shared, or all of them?
#  - how many columns to report: only most specific, or all of them?
#  - provide commands to
#       add info from A to B, overriding B?
#       add info from A to B, only when consistent with B?
#  - should diff report be threaded?

# TBD: filter out taxonID if idspaces are different

def differences(uid1, uid2, props = None):  # mask
  (add, change, drop) = differences_in_record(uid1, uid2, props = None)
  if len(cl.get_children(uid1)) != len(cl.get_children(uid2)):
    change |= 1 << number_of_children.specificity
  return (add, change, drop)

def differences_in_record(uid1, uid2, props = None):  # mask
  if props == None:
    props = table.get_table(uid2).properties
  drop = 0
  change = 0
  add = 0
  for prop in props:
    if prop and prop != taxonID and prop != parentNameUsageID:
      v1 = cl.get_value(uid1, prop)
      v2 = cl.get_value(uid2, prop)
      if v1 != v2:
        if v1 == None:
          add |= 1 << prop.specificity
        elif v2 == None:
          drop |= 1 << prop.specificity
        else:
          change |= 1 << prop.specificity
  # TBD: compare parents ??
  return (drop, change, add)

number_of_children = property.by_name("number_of_children")
number_of_synonyms = property.by_name("number_of_synonyms")

def same(comparison):
  return comparison == no_diffs

# Returns a list of properties ... ?

def unpack(comparison):
  if same(comparison):
    return []
  else:
    (_, change, _) = comparison
    return unpack1(change)

def unpack1(mask):
  props = []
  for prop in property.properties_by_specificity:
    spec = prop.specificity
    if (mask & (1 << spec)) != 0:
      props.append(prop)
  return props
