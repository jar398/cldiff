"""
Command line arguments:
  - checklist 1
  - checklist 2

Read the two checklists, call them A and B.
Make indexes by the values that are to be used in matching
  (all names, or just scientific names).
Go through A in preorder.
Find all matches for all children.
Collate into three groups:
  - unique match between A and B
  - ambiguous match, many As to many Bs
  - ambiguous match, one A to many Bs
  - ambiguous match, one B to many As
  - in A with no match in B
Followed by
  - children(?) of A-parent's match in B that have match in A (wait... ?)
and in the 'A and B' group distinguish unique matches from ambiguous matches.
Recur to the children of all three types.
"""

import sys, csv
import argparse
from checklist import *

def main(c1, c2, out):
  A = read_checklist(c1, "A.")
  B = read_checklist(c2, "B.")
  print ("counts:", len(get_all_tnus(A)), len(get_all_tnus(B)))

  analyze_fringes(A, B)
  print ("fringe matches:", len(fringe_matches))

  analyze_topologies(A, B)
  print ("cross_mrcas:", len(cross_mrcas))

  (best_A_for_B, reasons) = choose_best_matches(A, B)

  report(A, B, best_A_for_B, reasons, out)

# For every B TNU, whether accepted or synonym, determine a unique A
# TNU (either accepted or synonym) that's a best match for it.

def choose_best_matches(A, B):
  # For each higher taxon in B, find the higher taxon t in A that is the MRCA of 
  # all taxa in A that are matched to descendants of B.
  # That t is then a place to put taxa in B that are unassigned.

  # So... we need to be able take mrcas of nodes in A...
  # so, need to keep track of depths...

  best_A_for_B = {}
  reasons = {}

  def process(B_tnu):           # B_tnu is accepted

    # Look for match
    links = all_matches(B_tnu, A)
    if len(links) > 0:
      (A_tnu, synonym, comment) = links[0]
      best_A_for_B[B_tnu] = A_tnu
      reasons[B_tnu] = (25, comment)
    else:
      A_tnu = None
      reasons[B_tnu] = (99, "no match in A")

    # Recur, and find orphans/grafts
    for B_child in get_children(B_tnu):
      A_child = process(B_child)
      if not A_child:
        if A_tnu:
          best_A_for_B[B_child] = A_tnu   # parent/child, though
          reasons[B_child] = (80, "graft")
        else:
          reasons[B_child] = (85, "part of graft")

    return A_tnu

  for root in get_roots(B):
    process(root)
  print ("As that are best for some B:", len(best_A_for_B))
  return (best_A_for_B, reasons)


# Find unique match by text (name) in A for B_tnu (ignore ambiguous matches)

def text_match(A, B_tnu):
  A_candidates = get_tnus_with_value(A, canonical_name_field, get_name(B_tnu))
  if len(A_candidates) == 1:
    return A_candidates[0]
  else:
    return None

def match_order(B_tnu, reasons):
  (order, _) = reasons[B_tnu]
  return order

# Write checklist comparison report

def report(A, B, best_A_for_B, reasons, outpath):
  Bs_for_A = invert_dict(best_A_for_B)

  seen = {}
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)

    writer.writerow(["nesting", "A_taxon", "basis",
                     "B_taxon", "remark"])

    # Print in order of A hierarchy
    # Rel = B_tnu's relationship to A_tnu
    # A_tnu is accepted

    def descend(A_tnu, depth):
      B_tnus = Bs_for_A.get(A_tnu, ())

      # Show matches
      # Show A children
      # Show B grafts

      show_matches(A_tnu, B_tnus, depth)
      subdepth = depth + 1
        
      # Show children
      for child in get_children(A_tnu):
        descend(child, subdepth)

      # Show grafts (inferred children)
      # [B_tnu for B_tnu in B_tnus if is_graft(B_tnu)]
      for B_tnu in B_tnus:
        if is_graft(B_tnu):
          descend_graft(B_tnu, subdepth, True)

    def descend_graft(B_tnu, depth, graftp):
      # part of graft
      write_row(None, B_tnu, "placed with sibling(s)" if graftp else "", 
                depth)
      subdepth = depth + 1
      for child in get_children(B_tnu):
        descend_graft(child, subdepth, False)

    def show_matches(A_tnu, B_tnus, depth):

      B_tnus = \
        [B_tnu for B_tnu in B_tnus if not is_graft(B_tnu)]
      B_tnus = \
        sorted(B_tnus, key=lambda B_tnu: match_order(B_tnu, reasons))

      if len(B_tnus) == 1:
        B_tnu = B_tnus[0]
        write_row(A_tnu, B_tnu, relationship(A_tnu, B_tnu, "%1"),
                  depth)
      elif len(B_tnus) > 1:
        write_row(A_tnu, None, "%s matches follow:" % len(B_tnus),
                  depth)
        for B_tnu in B_tnus:
          write_row(None, B_tnu,
                    relationship(A_tnu, B_tnu, "%2"),
                    depth)
      else:
        write_row(A_tnu, None, "no match in B", depth)

    def is_graft(B_tnu):
      return reasons[B_tnu][0] >= 80

    def write_row(A_tnu, B_tnu, rel, depth):
      remark = None
      B_name = ''
      if B_tnu:
        seen[B_tnu] = True
        B_name = get_name(B_tnu)
        if A_tnu and get_name(A_tnu) == B_name:
          B_name = '='
        if not is_accepted(B_tnu):
          B_accepted = get_accepted(B_tnu)
          seen[B_accepted] = True
          remark = ("synonym of %s" %
                    (display_name(B_accepted)))
      writer.writerow([str(depth),
                       display_name(A_tnu),
                       rel,
                       display_name(B_tnu),
                       remark])

    def relationship(A_tnu, B_tnu, rel):
      return rel + (reasons[B_tnu][1] or "= by name")

    for root in get_roots(A):
      descend(root, 1)

    for B_tnu in get_all_tnus(B):
      if is_accepted(B_tnu):
        if not B_tnu in seen:
          p = get_parent(B_tnu)
          probe = reasons.get(B_tnu)
          write_row(best_A_for_B.get(B_tnu, None), B_tnu,
                    "fell through cracks; parent = %s; %s" %
                    (get_value(p, tnu_id_field) if p else "?",
                     probe[1] if probe else "??"),
                    1)

def display_name(tnu):
  if tnu:
    return get_unique(tnu)
  return ''


# ----------------------------------------------------------------------
# NEW

def compose(link1, link2):
  (node1, synonym1, comment1) = link1
  (node2, synonym2, comment2) = link2
  if badness(synonym1) < badness(synonym2):
    return link2
  else:
    return (node2, synonym1, comment1)

def synonym_matches(node):
  return [(node, None, None)] + \
         [(synonym, synonym, "<= by synonym")
            for synonym in get_synonyms(node)]
# tbd: use status = get_value(synonym, nomenclatural_status_field)

def direct_name_matches(node, other):
  name = get_name(node)
  candidates = get_tnus_with_value(other, canonical_name_field, name)
  return [(candidate, None, None)
            for candidate in candidates]

def accepted_matches(node):
  if is_accepted(node):
    return [(node, None, None)]
  else:
    return [(get_accepted(node), node, ">= by synonym")]

def name_based_matches(node1, other):
  links = []
  assert node1 + 0 == node1
  for link1 in synonym_matches(node1):
    (node2, _, comment2) = link1
    for link2 in direct_name_matches(node2, other):
      link12 = compose(link1, link2)
      (node3, _, comment3) = link12
      for link3 in accepted_matches(node3):
        links.append(compose(link12, link3))
  links = sorted(links, key=lambda link: badness(link[1]))
  pruned = []
  seen = []
  for link in links:
    if not link[0] in seen:
      pruned.append(link)
      seen.append(link)
  return pruned

def all_matches(node, other):
  m = topological_match(node)
  m2 = name_based_matches(node, other)
  if m:
    return [m] + m2
  else:
    return m2

#

def analyze_topologies(A, B):
  analyze_topology(A, B)
  analyze_topology(B, A)

def analyze_topology(checklist, other):
  for root in get_roots(checklist):
    subanalyze_topology(root, other)

cross_mrcas = {}

def subanalyze_topology(tnu, other):
  children = get_children(tnu)
  if len(children) > 0:
    m = None
    for child in children:
      partner = subanalyze_topology(child, other)
      if partner:
        m = mrca(m, get_parent(partner))
    if m:
      cross_mrcas[tnu] = m
      return m
  return mutual_fringe_match(tnu)

def topological_match(tnu):
  partner = cross_mrcas.get(tnu)
  if partner:
    back = cross_mrcas.get(partner)
    if tnu == back:
      return (partner, None, "= by topology")
    else:
      return (partner, None, "> by topology")
  else:
    return None

#

def analyze_fringes(A, B):
  analyze_fringe(A, B)
  analyze_fringe(B, A)

def analyze_fringe(checklist, other):
  for root in get_roots(checklist):
    subanalyze_fringe(root, other)

fringe_matches = {}

def subanalyze_fringe(tnu, other):
  children = get_children(tnu)
  if len(children) > 0:
    found_match = False
    for child in children:
      if subanalyze_fringe(child, other):
        found_match = True
    if found_match:
      return found_match
  matches = name_based_matches(tnu, other)
  if len(matches) > 0:
    fringe_matches[tnu] = matches[0][0]
    return True
  return False

def mutual_fringe_match(tnu):
  partner = fringe_matches.get(tnu)
  if partner and fringe_matches.get(partner) == tnu:
    return partner
  return None

# For assigning priorities to synonyms

def badness(tnu):
  if tnu == None: return 0
  if is_accepted(tnu): return 0
  status = get_value(tnu, nomenclatural_status_field)
  if status is None:
    return 99
  badness = badnesses.get(status, None)
  if badness is None: badness = 99
  return badness

# Name / match classes, best to worst

badnesses = {
  "identical names": 1,        # namestrings are the same
  "authority": 1.5,
  "scientific name": 2,        # (actually canonical) exactly one per node
  "equivalent name": 3,        # synonym but not nomenclaturally
  "misspelling": 3.8,
  "genbank synonym": 4,        # at most one per node; first among equals
  "anamorph": 4.1,
  "teleomorph": 4.2,
  "unpublished name": 4.5,    # non-code synonym
  "id": 4.7,
  "merged id": 4.8,

  # above here: equivalence implied. below here: acc>=syn implied.
  # except in the case if 'in-part' which is acc<syn.

  "synonym": 5,
  "misnomer": 5.5,
  "includes": 6,
  "in-part": 6.5,              # this node is part of a polyphyly
  "type material": 7,
  "blast name": 8,             # large well-known taxa
  "genbank common name": 9,    # at most one per node
  "genbank acronym": 9.2,      # at most one per node
  "genbank anamorph": 9.4,     # at most one per node
  "common name": 10,
  "acronym": 10.5,
}

# When invoked from command line:

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
