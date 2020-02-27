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

def main(c1, c2):
  A = read_checklist(c1)
  B = read_checklist(c2)
  print ("counts:", len(A), len(B))

  A_text_index = index_by_column(A, "canonicalName")
  B_synonym_index = index_by_column(B, "acceptedNameUsageID")
  A_id_index = index_unique_by_column(A, "taxonID")

  print("A indexed by text:", len(A_text_index))

  As_for_Bs = {}

  attempts = [0]
  successes = [0]

  # Match as many by name as possible
  for B_tnu in B:
    if is_accepted(B_tnu):

      def attempt_match(B_tnu):
        # Find all the A-tnus with same name as this B-tnu...
        attempts[0] += 1
        text = get_value(B_tnu, "canonicalName")
        A_synonyms = A_text_index.get(text, ())
        if len(A_synonyms) == 1:
          A_tnu = A_synonyms[0]
          A_accepted_id = get_value(A_tnu, "acceptedNameUsageID")
          if A_accepted_id != None:
            A_tnu = A_id_index[A_accepted_id]
          As_for_Bs[B_tnu] = A_tnu
          successes[0] += 1
          return True
        else:
          return False

      if attempt_match(B_tnu):
        pass
      else:
        for B_synonym in sorted(get_synonyms(B_tnu, B_synonym_index), key=badness):
          if attempt_match(B_synonym):
            break

  print ("Attempts:", attempts[0])
  print ("Successes:", successes[0])

  # TBD: match by topology
  #   (a) those that didn't get assignments previously
  #   (b) higher nodes including those that already have name matches

  Bs_for_As = invert_dict(As_for_Bs)
  print ("Bs for As:", len(Bs_for_As))

  # TBD: generate a report
  report(A, Bs_for_As, As_for_Bs)

def report(A, Bs_for_As, As_for_Bs):
  outpath = "diff.csv"
  A_id_index = index_unique_by_column(A, "taxonID")
  A_hierarchy = index_hierarchy(A)
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)

    writer.writerow(["A_id", "A_name", "B_id", "B_name", "how"])

    def write_row(A_tnu, B_tnu, how):
      writer.writerow([(get_value(A_tnu, "taxonID") if A_tnu  else ""),
                       (get_value(A_tnu, "canonicalName") if A_tnu else ""),
                       (get_value(B_tnu, "taxonID") if B_tnu else ""),
                       (get_value(B_tnu, "canonicalName") if B_tnu else ""),
                       how])

    def descend(A_tnu):
      A_id = get_value(A_tnu, "taxonID")
      B_tnus = Bs_for_As.get(A_tnu, None)
      if B_tnus:
        if len(B_tnus) == 1:
          write_row(A_tnu, B_tnus[0], "one-one")
        else:
          write_row(A_tnu, None, "one-many")
          # Iterate over A-synonyms ???
          # sorted(get_synonyms(A_tnu, A_synonym_index), key=badness)
          for B_tnu in B_tnus:
            write_row(None, B_tnu, "")
      else:
        write_row(A_tnu, None, "not in B")
      for child in get_children(A_tnu, A_hierarchy):
        descend(child)

    for root in get_roots(A, A_id_index):
      descend(root)

# Before calling this, cache synonym_index = index_by_column(tnus, "acceptedNameUsageID").

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
    (guide, row) = registry[uid]
    value = row[guide[label]]
    if value != '':
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
    (guide, row) = registry[uid]
    value = row[guide[label]]
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
  main(sys.argv[1], sys.argv[2])
