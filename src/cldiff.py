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

    # 1. Try to put it with its siblings (n.b. only accepted TNUs have children)
    A_tnu = None
    pending = []
    for B_child in get_children(B_tnu):
      A_child = process(B_child)
      if A_child:
        A_tnu = mrca(get_parent(A_child), A_tnu)
      else:
        pending.append(B_child) # graft
    if A_tnu:
      best_A_for_B[B_tnu] = A_tnu   # Accepted
      reasons[B_tnu] = (10, "by topology")

    # 2. If no siblings got assigned, look for textual match to
    # accepted name.
    if A_tnu == None:
      A_tnu = text_match(A, B_tnu)
      if A_tnu:
        best_A_for_B[B_tnu] = A_tnu
        if is_accepted(A_tnu):
          reasons[B_tnu] = (20, "accepted/accepted")
        else:
          # This is the failing case??
          reasons[B_tnu] = (25, "synonym/accepted")

    # Now process the synonyms, looking for textual matches in A.
    # (? maybe only do this if no match so far ?)
    seen = [A_tnu]
    for B_syn in get_synonyms(B_tnu): # Sorted by quality
      A_match = text_match(A, B_syn)
      if A_match:
        A_accepted = get_accepted(A_match)
        if A_accepted in seen:   # Eliminate redundant low quality synonyms
          reasons[B_syn] = (98, "redundant B-synonym")
        else:
          seen.append(A_accepted)
          best_A_for_B[B_syn] = A_match
          if is_accepted(A_match):
            reasons[B_syn] = (30, "accepted/synonym")
          else:
            reasons[B_syn] = (35, "synonym/synonym")

    # 3. No match of accepted name; use mrca of synonym matches (if any)
    #  (alternatively, 'best' synonym ?)
    if A_tnu == None:
      for B_syn in get_synonyms(B_tnu):
        A_match = text_match(A, B_syn)
        A_tnu = mrca(A_match, A_tnu)          # ???
      if A_tnu:
        best_A_for_B[B_tnu] = A_tnu
        reasons[B_tnu] = (50, "mrca of B-synonyms")

    # Graft unmatched siblings to mrca of matched siblings
    if A_tnu:
      for p in pending:
        best_A_for_B[p] = A_tnu
        reasons[p] = (80, "graft")
    else:
      for p in pending:
        reasons[p] = (85, "part of graft")

    if not B_tnu in best_A_for_B:
      reasons[B_tnu] = (99, "no match/accepted")

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

    writer.writerow(["nesting", "A_id", "A_name", "relationship",
                     "B_id", "B_name", "remark"])

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
      for A_syn in get_synonyms(A_tnu):
        B_tnus += Bs_for_A.get(A_syn, ())
      B_tnus = \
        sorted(B_tnus, key=lambda B_tnu: match_order(B_tnu, reasons))

      if len(B_tnus) == 1:
        B_tnu = B_tnus[0]
        write_row(A_tnu, B_tnu, relationship(A_tnu, B_tnu, "unique "),
                  depth)
      elif len(B_tnus) > 1:
        write_row(A_tnu, None, "%s matches follow:" % len(B_tnus),
                  depth)
        for B_tnu in B_tnus:
          write_row(None, B_tnu,
                    relationship(A_tnu, B_tnu, ""),
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
      return rel + reasons[B_tnu][1]

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


# Common ancestor - utility

def mrca(tnu1, tnu2):
  if tnu1 == None: return tnu2
  if tnu2 == None: return tnu1
  if tnu1 == tnu2: return tnu1

  tnu1 = get_accepted(tnu1)
  tnu2 = get_accepted(tnu2)

  d1 = get_depth(tnu1)
  d2 = get_depth(tnu2)
  while d1 > d2:
    tnu1 = get_parent(tnu1)
    d1 -= 1
  while d2 > d1:
    tnu2 = get_parent(tnu2)
    d2 -= 1
  while tnu1 != tnu2:
    tnu1 = get_parent(tnu1)
    tnu2 = get_parent(tnu2)
  return tnu1

depth_cache = {}

def get_depth(tnu):
  depth = depth_cache.get(tnu, None)
  if depth: return depth
  parent = get_parent(tnu)
  if parent == None:
    d = 0
  else:
    d = get_depth(parent) + 1
  depth_cache[tnu] = d
  return d


# When invoked from command line:

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
