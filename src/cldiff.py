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

import os, sys, csv
import argparse

def main(c1, c2, out):
  A = read_checklist(c1)
  B = read_checklist(c2)
  print ("counts:", len(A), len(B))

  (As_for_Bs, routes) = match_by_name(A, B)
  (grafts, by_topology) = match_by_topology(A, B, As_for_Bs)

  report(A, B, As_for_Bs, routes, grafts, by_topology, out)

# TBD: Also match by scientificName

def match_by_name(A, B):
  A_text_index = index_by_column(A, "canonicalName")
  B_synonym_index = index_by_column(B, "acceptedNameUsageID")
  A_id_index = index_unique_by_column(A, "taxonID")
  print("A indexed by text:", len(A_text_index))
  As_for_Bs = {}
  routes = {}
  attempts = [0]
  # Match as many by name as possible
  for B_tnu in B:
    if is_accepted(B_tnu):
      B_candidates = [B_tnu] + sorted(get_synonyms(B_tnu, B_synonym_index),
                                      key=badness)
      for B_candidate in B_candidates:
        # Find all the A-tnus with same name as this B-tnu...
        attempts[0] += 1
        text = get_value(B_candidate, "canonicalName")
        A_candidates = A_text_index.get(text, ())
        if len(A_candidates) == 1:
          # Many different ways to design this.  Prioritize etc.
          A_candidate = A_candidates[0]
          # TBD: Keep track of nomenclaturalStatus of A_candidate and B_candidate
          A_accepted = A_candidate
          if not is_accepted(A_candidate):
            A_accepted_id = get_value(A_candidate, "acceptedNameUsageID")
            if A_accepted_id in A_id_index:
              A_accepted = A_id_index[A_accepted_id]
            else:
              print ("accepted but no accepted id", A_candidate)
          As_for_Bs[B_tnu] = A_accepted
          routes[B_tnu] = (A_accepted, A_candidate, B_candidate)
          # This candidate matched; no need to look for any others.
          break
  print ("Attempts to match a B:", attempts[0])
  print ("A's that got B matches:", len(As_for_Bs))
  return (As_for_Bs, routes)

# Place orphaned B nodes according to where their siblings got placed.

def match_by_topology(A, B, As_for_Bs):
  # For each higher taxon in B, find the higher taxon t in A that is the MRCA of 
  # all taxa in A that are matched to descendants of B.
  # That t is then a place to put taxa in B that are unassigned.

  # So... we need to be able take mrcas of nodes in A...
  # so, need to keep track of depths...

  B_hierarchy = index_hierarchy(B)
  A_id_index = index_unique_by_column(A, "taxonID")
  depth_cache = {}

  grafts = {}                   # B -> A (new parent where it gets grafted)
  by_topology = {}

  def process(B_tnu):
    m = None                    # m will be a node in A
    pending = []
    clumped = False
    for B_child in get_children(B_tnu, B_hierarchy):
      x = process(B_child)
      if x:
        # B_child has an A-mrca (contains descendants matched to A)
        if m:
          m2 = mrca(m, x, A_id_index, depth_cache)
          if m2 != m and m2 != x:
            clumped = True
          m = m2
        else:
          m = x
      else:
        # These might become grafts!!!
        pending.append(B_child)
    # end of loop to process children.
    A_tnu = None
    if m:
      if clumped:
        A_tnu = m
      else:
        A_tnu = get_parent(m, A_id_index)
      by_topology[B_tnu] = A_tnu
    else:
      # No topological match.  Look for name match.
      A_tnu = As_for_Bs.get(B_tnu, None)   # might be None
    if A_tnu:
      for p in pending:
        grafts[p] = A_tnu
    return A_tnu
  B_id_index = index_unique_by_column(B, "taxonID")
  for root in get_roots(B, B_id_index):
    process(root)
  print ("Grafts:", len(grafts))
  print ("By topology:", len(by_topology))
  return (grafts, by_topology)

def mrca(tnu1, tnu2, id_index, depth_cache):
  d1 = get_depth(tnu1, id_index, depth_cache)
  d2 = get_depth(tnu2, id_index, depth_cache)
  while d1 > d2:
    tnu1 = get_parent(tnu1, id_index)
    d1 -= 1
  while d2 > d1:
    tnu2 = get_parent(tnu2, id_index)
    d2 -= 1
  while tnu1 != tnu2:
    tnu1 = get_parent(tnu1, id_index)
    tnu2 = get_parent(tnu2, id_index)
  return tnu1

def get_depth(tnu, id_index, depth_cache):
  depth = depth_cache.get(tnu, None)
  if depth: return depth
  parent = get_parent(tnu, id_index)
  if parent == None:
    d = 0
  else:
    d = get_depth(parent, id_index, depth_cache) + 1
  depth_cache[tnu] = d
  return d

def report(A, B, As_for_Bs, routes, grafts, by_topology, outpath):
  Bs_for_As = invert_dict(As_for_Bs)
  print ("Bs for As:", len(Bs_for_As))
  B_grafts_for_As = invert_dict(grafts)
  B_topological_for_As = invert_dict(by_topology)
  A_id_index = index_unique_by_column(A, "taxonID")
  A_hierarchy = index_hierarchy(A)
  B_hierarchy = index_hierarchy(B)
  seen = {}
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)

    writer.writerow(["nesting", "A_id", "A_name", "B_id", "B_name", "how", "mode"])

    def write_row(A_tnu, B_tnu, how, depth):
      if B_tnu:
        seen[B_tnu] = True
        route = routes.get(B_tnu, None)
        if route:
          (A_accepted, A_candidate, B_candidate) = routes[B_tnu]
          if A_accepted == A_candidate and B_tnu == B_candidate:
            mode = "direct"
          elif A_accepted == A_candidate and B_tnu != B_candidate:
            mode = "via synonym in B"
          elif A_accepted != A_candidate and B_tnu == B_candidate:
            mode = "via synonym in A"
          else:
            mode = ("via synonym in both: %s" %
                    get_name(A_candidate))
        else:
          mode = None           # Not in A
      else:
        mode = None             # Not in B
      writer.writerow([str(depth),
                       ('A:' + get_value(A_tnu, "taxonID") if A_tnu  else ''),
                       (get_name(A_tnu) if A_tnu else ""),
                       ('B:' + get_value(B_tnu, "taxonID") if B_tnu else ''),
                       (get_name(B_tnu) if B_tnu else ''),
                       how,
                       mode])

    def tweak_how(B_tnu, how, by_topo, by_text):
      if B_tnu in by_topo:
        if B_tnu in by_text:
          return how
        else:
          return how + " (topology only)"
      else:
          return how + " (text only)"

    # Print in order of A hierarchy

    def descend(A_tnu, depth):
      subdepth = depth + 1
      A_id = get_value(A_tnu, "taxonID")

      B_tnus = []
      B_topos = B_topological_for_As.get(A_tnu, ())
      B_textuals = Bs_for_As.get(A_tnu, ())
      for B_tnu in B_topos:
        B_tnus.append(B_tnu)
      for B_tnu in B_textuals:
        if not B_tnu in B_tnus:
          B_tnus.append(B_tnu)
      if len(B_tnus) == 1:
        B_tnu = B_tnus[0]
        write_row(A_tnu, B_tnu, tweak_how(B_tnu, "one-one", B_topos, B_textuals),
                  depth)
      elif len(B_tnus) > 1:
        write_row(A_tnu, None, "one-many", depth)
        for B_tnu in B_tnus:
          write_row(A_tnu, B_tnu,
                    tweak_how(B_tnu, "one of many", B_topos, B_textuals),
                    subdepth)
      else:
        write_row(A_tnu, None, "no match in B", depth)
        
      for child in get_children(A_tnu, A_hierarchy):
        descend(child, subdepth)
      for B_tnu in B_grafts_for_As.get(A_tnu, ()):
        descend_graft(B_tnu, subdepth)

    def descend_graft(B_tnu, depth):
      subdepth = depth + 1
      if not B_tnu in seen:
        write_row(None, B_tnu, "no match in A", depth)
        for child in get_children(B_tnu, B_hierarchy):
          descend_graft(child, subdepth)

    for root in get_roots(A, A_id_index):
      descend(root, 1)

    for B_tnu in B:
      if is_accepted(B_tnu):
        if not B_tnu in seen:
          write_row(None, B_tnu, "fell through cracks", 1)

# The following is for display purposes

def get_name(tnu):
  name = get_value(tnu, "canonicalName")
  if name != None:
    return name
  return get_value(tnu, "scientificName")

# Before calling this, cache as follows:
#   synonym_index = index_by_column(tnus, "acceptedNameUsageID").

def get_synonyms(tnu, synonym_index):
  id = get_value(tnu, "taxonID")
  if id is None: return ()
  return synonym_index.get(id, ())

def is_accepted(tnu):
  return get_value(tnu, "taxonomicStatus") == "accepted"

# Roots - accepted tnus without parents

def get_roots(tnus, id_index):
  roots = []
  for tnu in tnus:
    if is_accepted(tnu):
      parent = get_parent(tnu, id_index)
      if parent is None:
        roots.append(tnu)
  print (len(roots), "roots")
  return roots

# Parent/children

def get_parent(tnu, id_index):
  id = get_value(tnu, "parentNameUsageID")
  if id != None:
    return id_index.get(id, None)
  else:
    return None

# List of child tnus, or () if none

def get_children(tnu, hierarchy):
  parent_id = get_value(tnu, "taxonID")
  return hierarchy.get(parent_id, ())

# For each parent tnu id, all tnus having that id as taxon id

def index_hierarchy(tnus):
  return index_by_column(tnus, "parentNameUsageID")

def badness(tnu):
  status = get_value(tnu, "nomenclaturalStatus")
  if status is None:
    status = get_value(tnu, "taxonomicStatus")
    if status is None:
      return 99
  badness = badnesses.get(status, None)
  if badness is None: badness = 99
  return badness

# Name classes, best to worst

badnesses = {
  "authority": 0,
  "scientific name": 1,        # exactly one per node
  "accepted": 1,
  "equivalent name": 2,        # synonym but not nomenclaturally
  "misspelling": 3,
  "genbank synonym": 4,        # at most one per node; first among equals
  "synonym": 5,
  "anamorph": 5.1,
  "teleomorph": 5.2,
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
  "unpublished name": 10.7,
  "id": 11,
  "merged id": 12,
}

# Totally general utilities from here down... I guess...

def invert_dict(d):
  inv = {}
  for (key, val) in d.items():
    if val in inv:
      inv[val].append(key)
    else:
      inv[val] = [key]
  return inv

registry = {0: "unused"}

def read_checklist(indir):
  uids = []
  with open(get_tnu_path(indir), "r") as infile:
    reader = csv.reader(infile, delimiter="\t", quotechar='\a', quoting=csv.QUOTE_NONE)
    head = next(reader)
    guide = {key: index for (key, index) in zip(head, range(len(head)))}
    for row in reader:
      uid = len(registry)
      registry[uid] = (guide, row)
      uids.append(uid)
  return uids

def get_value(uid, label):
  (guide, row) = registry[uid]
  text = row[guide.get(label, None)]
  if text == '':
    return None
  return text

# Create a dict mapping field values to lists of keys (uids, tnus)

def index_by_column(checklist, label):
  index = {}
  for uid in checklist:
    value = get_value(uid, label)
    if value != None:
      have = index.get(value, None)
      if have:
        have.append(uid)
      else:
        index[value] = [uid]
  print (label, "entries:", len(index))
  return index

# You could use this for, say, taxonID, or maybe authority (scientificName)

def index_unique_by_column(checklist, label):
  index = {}
  for uid in checklist:
    value = get_value(uid, label)
    if value != '':
      if label in index:
        raise ValueError("conflict in unique index: %s -> %s, %s" % (value, index[value], uid))
      index[value] = uid
  return index

# Utility - copied from another file - really ought to be shared

def get_tnu_path(dwca_dir):
  for name in ["taxon.tsv",
               "Taxon.tsv",
               "taxon.txt",
               "Taxon.txt"]:
    path = os.path.join(dwca_dir, name)
    if os.path.exists(path):
      return path
  raise ValueError("cannot find TNU file", dwca_dir)


# From command line

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
