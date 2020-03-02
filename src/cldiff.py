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

  (A_by_name_for_B, routes) = match_by_name(A, B)
  (grafts, A_by_topology_for_B) = match_by_topology(A, B, A_by_name_for_B)

  report(A, B, A_by_name_for_B, routes, A_by_topology_for_B, grafts, out)

# TBD: Also match by scientificName

def match_by_name(A, B):
  A_text_index = index_by_column(A, canonical_name_field)
  B_synonym_index = index_by_column(B, accepted_tnu_id_field)
  A_id_index = index_unique_by_column(A, tnu_id_field)
  print("A indexed by text:", len(A_text_index))
  A_by_name_for_B = {}
  routes = {}
  attempts = [0]
  # Match as many by name as possible

  # Find all B->A exact text matches
  A_by_text_for_B = {}
  for B_tnu in B:
    A_candidates = A_text_index.get(get_name(B_tnu), ())
    if len(A_candidates) == 1:
      A_by_text_for_B[B_tnu] = A_candidates[0]
  print ("B's with at least one A match:", len(A_by_text_for_B))

  A_by_name_for_B = {}

  # Look for synonym - text - synonym routes between accepted TNUs
  for B_tnu in B:
    if is_accepted(B_tnu):
      B_accepted = B_tnu
      B_candidates = [B_accepted] + get_synonyms(B_tnu, B_synonym_index)
      for B_candidate in B_candidates:
        # Find all the A-tnus with same name as this B-tnu...
        attempts[0] += 1
        A_candidate = A_by_text_for_B.get(B_candidate, None)
        if A_candidate:
          A_accepted = get_accepted(A_candidate, A_id_index)
          if False:
            A_by_name_for_B[B_accepted] = A_candidate
          else:
            A_by_name_for_B[B_accepted] = A_accepted
          routes[B_accepted] = (A_accepted, A_candidate, B_candidate, B_accepted)
          # This candidate matched; no need to look for any others.
          break

  for B_tnu, A_tnu in A_by_name_for_B.items():
    if B_tnu % 50 == 0: print (B_tnu, A_tnu)

  print ("Attempts to match an accepted B:", attempts[0])
  print ("A's that got B matches:", len(A_by_name_for_B))
  return (A_by_name_for_B, routes)

# Place orphaned B nodes according to where their siblings got placed.

def match_by_topology(A, B, A_by_name_for_B):
  # For each higher taxon in B, find the higher taxon t in A that is the MRCA of 
  # all taxa in A that are matched to descendants of B.
  # That t is then a place to put taxa in B that are unassigned.

  # So... we need to be able take mrcas of nodes in A...
  # so, need to keep track of depths...

  B_hierarchy = index_hierarchy(B)
  A_id_index = index_unique_by_column(A, tnu_id_field)
  depth_cache = {}

  grafts = {}                   # B -> A (new parent where it gets grafted)
  A_by_topology_for_B = {}

  def process(B_tnu):
    m = None                    # m will be a node in A
    pending = []
    clumped = False
    for B_child in get_children(B_tnu, B_hierarchy):
      x = process(B_child)
      if x:
        # B_child has an A-mrca (contains descendants matched to A)
        if m:
          # TBD: Toss out children that are in incertae sedis/unclassified containers?
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
      A_by_topology_for_B[B_tnu] = A_tnu
      for p in pending:
        grafts[p] = A_tnu
        # graft_kinds[p] = kind
      return A_tnu
    else:
      # No topological match.  Look for name match.
      return A_by_name_for_B.get(B_tnu, None)   # might be None
  B_id_index = index_unique_by_column(B, tnu_id_field)
  for root in get_roots(B, B_id_index):
    process(root)
  print ("Grafts:", len(grafts))
  print ("By topology:", len(A_by_topology_for_B))
  return (grafts, A_by_topology_for_B)

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

# Write checklist comparison report

def report(A, B, A_by_name_for_B, routes, A_by_topology_for_B, grafts, outpath):
  A_for_B = {}
  for B_tnu, A_tnu in A_by_name_for_B.items():
    A_for_B[B_tnu] = A_tnu
  for B_tnu, A_tnu in A_by_topology_for_B.items():
    A_for_B[B_tnu] = A_tnu  # override
  print ("As for Bs:", len(A_for_B))

  Bs_for_A = invert_dict(A_for_B)

  print ("Bs for As:", len(Bs_for_A))

  A_id_index = index_unique_by_column(A, tnu_id_field)
  A_hierarchy = index_hierarchy(A)
  B_hierarchy = index_hierarchy(B)
  A_synonym_index = index_by_column(A, accepted_tnu_id_field)
  B_grafts_for_A = invert_dict(grafts)
  seen = {}
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)

    writer.writerow(["nesting", "A_id", "A_name", "relationship", "B_id", "B_name", "match mode",
                     "synonymy route"])

    # Print in order of A hierarchy

    def descend(A_tnu, depth, mode):
      show_matches(A_tnu, depth, mode)
      subdepth = depth + 1
        
      # Also show synonym matches here
      for A_synonym in get_synonyms(A_tnu, A_synonym_index):
        # Does this synonym match any B?
        if A_synonym in Bs_for_A:
          # descend(A_synonym, subdepth, "synonym")
          show_matches(A_synonym, subdepth, "synonym")

      for child in get_children(A_tnu, A_hierarchy):
        descend(child, subdepth, "child")

      for B_tnu in B_grafts_for_A.get(A_tnu, ()):
        descend_graft(B_tnu, subdepth, True)

    def show_matches(A_tnu, depth, mode):
      B_tnus = Bs_for_A.get(A_tnu, ())
      if len(B_tnus) == 1:
        B_tnu = B_tnus[0]
        write_row(A_tnu, B_tnu, tweak_how(B_tnu, "unique match",
                                          A_by_topology_for_B, A_by_name_for_B),
                  mode, depth)
      elif len(B_tnus) > 1:
        write_row(A_tnu, None, "%s matches:" % len(B_tnus), mode, depth)
        subdepth = depth + 1
        for B_tnu in B_tnus:
          write_row(None, B_tnu,
                    tweak_how(B_tnu, "match",
                              A_by_topology_for_B, A_by_name_for_B),
                    "match", subdepth)
      else:
        write_row(A_tnu, None, "no name match in B", mode, depth)

    def descend_graft(B_tnu, depth, graftp):
      subdepth = depth + 1
      # part of graft
      write_row(None, B_tnu, "placed with sibling(s)" if graftp else "", "graft", depth)
      for child in get_children(B_tnu, B_hierarchy):
        descend_graft(child, subdepth, False)

    def write_row(A_tnu, B_tnu, how, mode, depth):
      if B_tnu:
        seen[B_tnu] = True
        route = routes.get(B_tnu, None)
        if route:
          (A_accepted, A_candidate, B_candidate, B_accepted) = route
          A_status = get_value(A_candidate, nomenclatural_status_field)
          B_status = get_value(B_candidate, nomenclatural_status_field)
          if A_status: A_status = " (%s)" % A_status
          else: A_status = ""
          if B_status: B_status = " (%s)" % B_status
          else: B_status = ""
          if A_accepted != A_tnu:
            path = None
          elif A_accepted == A_candidate and B_accepted == B_candidate:
            path = "same name"
          elif A_accepted == A_candidate and B_accepted != B_candidate:
            path = "A name is synonym in B" + B_status
          elif A_accepted != A_candidate and B_accepted == B_candidate:
            path = "B name is synonym in A" + A_status
          elif get_name(A_candidate) == get_name(B_candidate):
            if A_status == B_status: B_status = ""
            path = ("A and B both have synonym %s%s%s" %
                    (get_name(A_candidate), A_status, B_status))
          else:
            path = "not by name"
        else:
          path = None           # Not in A, or no A given
      else:
        path = None             # Not in B
      B_name = ''
      if B_tnu:
        B_name = get_name(B_tnu)
        if A_tnu and get_name(A_tnu) == B_name:
          B_name = '='
      writer.writerow([str(depth),
                       display_id(A_tnu, "A"),
                       (get_name(A_tnu) if A_tnu else ""),
                       mode,
                       display_id(A_tnu, "B"),
                       B_name,
                       how,
                       path])

    def display_id(tnu, prefix):
      if tnu:
        id = get_value(tnu, tnu_id_field)
        if id:
          return "%s:%s" % (prefix, id)
      return ''

    def tweak_how(B_tnu, how, by_topo, A_by_name_for_B):
      if B_tnu in by_topo:
        if how != "": how = how + " "
        if B_tnu in A_by_name_for_B:
          return how + "(by topology)"
        else:
          return how + "(by topology only)"
      else:
          return how

    for root in get_roots(A, A_id_index):
      descend(root, 1, "root")

    for B_tnu in B:
      if is_accepted(B_tnu):
        if not B_tnu in seen:
          write_row(None, B_tnu, "fell through cracks", "oops", 1)

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
  id = get_value(tnu, parent_tnu_id_field)
  if id != None:
    return id_index.get(id, None)
  else:
    return None

# List of child tnus, or () if none

def get_children(tnu, hierarchy):
  parent_id = get_value(tnu, tnu_id_field)
  return hierarchy.get(parent_id, ())

# For each parent tnu id, all tnus having that id as taxon id

def index_hierarchy(tnus):
  return index_by_column(tnus, parent_tnu_id_field)

def badness(tnu):
  status = get_value(tnu, nomenclatural_status_field)
  if status is None:
    status = get_value(tnu, taxonomic_status_field)
    if status is None:
      return 99
  badness = badnesses.get(status, None)
  if badness is None: badness = 99
  return badness

# Name classes, best to worst

badnesses = {
  "authority": 0,
  "scientific name": 1,        # (actually canonical) exactly one per node
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
  "unpublished name": 10.7,    # non-code synonym
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

# The following is for display purposes

def get_name(tnu):
  name = get_value(tnu, canonical_name_field)
  if name != None:
    return name
  return get_value(tnu, scientific_name_field)

# Before calling this, cache as follows:
#   synonym_index = index_by_column(tnus, accepted_tnu_id_field).

def get_synonyms(tnu, synonym_index):
  id = get_value(tnu, tnu_id_field)
  if id is None: return ()
  return sorted(synonym_index.get(id, ()),
                key=badness)

def is_accepted(tnu):
  return get_value(tnu, taxonomic_status_field) == "accepted"

def get_accepted(A_candidate, A_id_index):
  if is_accepted(A_candidate):
    A_accepted = A_candidate
  else:
    A_accepted_id = get_value(A_candidate, accepted_tnu_id_field)
    if A_accepted_id in A_id_index:
      A_accepted = A_id_index[A_accepted_id]
    else:
      print ("accepted but no accepted id", A_candidate)
      A_accepted = None
  return A_accepted

# Get the value of a field of a TNU record

def get_value(uid, label):
  (guide, row) = registry[uid]
  text = row[guide.get(label, None)]
  if text == '':
    return None
  return text

tnu_id_field = "taxonID"
accepted_tnu_id_field = "acceptedNameUsageID"
taxonomic_status_field = "taxonomicStatus"
nomenclatural_status_field = "nomenclaturalStatus"
canonical_name_field = "canonicalName"    # without authority
scientific_name_field = "scientificName"  # with authority
parent_tnu_id_field = "parentNameUsageID"

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
  print (label, "values:", len(index))
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
               "taxa.txt",
               "taxon.txt",
               "Taxon.txt"]:
    path = os.path.join(dwca_dir, name)
    if os.path.exists(path):
      return path
  raise ValueError("cannot find TNU file", dwca_dir)


# When invoked from command line:

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
