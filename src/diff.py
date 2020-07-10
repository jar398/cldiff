debug = False

import table
import property
import checklist as cl

# Difference report.
# Each bit position is 1 for a mismatch, 0 for a match.
# Bit positions correspond to the URIs listed in property.py.
# The bit mask is a 'badness' measure and can be used to sort
# options for matching.

no_diffs = 0

# Among d1 and d2, the worst of both
# (this still isn't right)

def compose(d1, d2):
  return d1 | d2

# Among d1 and d2, whichever is better

def conjoin(d1, d2):
  return min(d1, d2)

tip_prop = property.by_name("tips")

# Fix later

def combine(d1, d2):
  return d1 | d2

# ---------- Record comparison.

# Degrees of freedom:
#  - columns to consider: only shared, or all of them?
#  - how many columns to report: only most specific, or all of them?
#  - provide commands to
#       add info from A to B, overriding B?
#       add info from A to B, only when consistent with B?
#  - should diff report be threaded?

def differences(uid1, uid2):  # mask
  (r1, t1) = table.record_and_table(uid1)
  (r2, t2) = table.record_and_table(uid2)

  comparison = 0
  for pos1 in range(len(t1.properties)):
    prop = t1.properties[pos1]
    pos2 = t2.position_index[prop.uid]
    if pos2 != None:
      # Property is provided in both tables
      v1 = r1[pos1]
      v2 = r2[pos2]
      if v1 != None and v2 != None:
        if v1 != v2:
          if debug:
            print("difference: %s %s %s" %
                  (cl.get_unique(uid1), prop.pet_name,
                   cl.get_unique(uid2)))
          comparison |= 1 << prop.specificity
  if len(cl.get_children(uid1)) != len(cl.get_children(uid2)):
    comparison |= 1 << number_of_children.specificity
  return comparison

number_of_children = property.by_name("number_of_children")
number_of_synonyms = property.by_name("number_of_synonyms")

# Returns a list of properties ... ?

def unpack(comparison):
  if comparison == 0:
    return []
  else:
    props = []
    for prop in property.properties_by_specificity:
      spec = prop.specificity
      if comparison & (1 << spec) != 0:
        props.append(prop)
    return props
