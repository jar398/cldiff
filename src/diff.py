
import sys, csv
import argparse

import checklist as cl
import relation as rel
import articulation as art

def main(c1, c2, out):
  A = cl.read_checklist(c1, "A.")
  B = cl.read_checklist(c2, "B.")
  print ("counts:", len(cl.get_all_tnus(A)), len(cl.get_all_tnus(B)))

  analyze_fringes(A, B)
  print ("number of fringe matches:", len(fringe_status))

  analyze_topologies(A, B)
  print ("number of cross-mrcas:", len(cross_mrcas))

  report(A, B, "diff.csv")

def report(A, B, outpath):
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)
    for root in cl.get_roots(A):
      subreport(root, A, B, writer)

def subreport(node, A, B, writer):
  art = mutual_fringe_match(node, B)
  if art:
    writer.writerow([cl.get_unique(art.dom),
                     art.rel.name,
                     cl.get_unique(art.cod),
                     art.comment])
  else:
    art = topological_match(node)
    if art:
      writer.writerow([cl.get_unique(art.dom),
                       art.rel.name,
                       cl.get_unique(art.cod),
                       art.comment])
    else:
        writer.writerow([cl.get_unique(node),
                         "",
                         "",
                         "no match in B"])
  for child in cl.get_children(node):
    subreport(child, A, B, writer)


# Topological matches (by tip containment)

def analyze_topologies(A, B):
  analyze_topology(A, B)
  analyze_topology(B, A)

def analyze_topology(checklist, other):
  for root in cl.get_roots(checklist):
    subanalyze_topology(root, other)

cross_mrcas = {}

def subanalyze_topology(tnu, other):
  children = cl.get_children(tnu)
  art = mutual_fringe_match(tnu, other)
  if art:
    return art.cod
  if len(children) > 0:
    m = None
    for child in children:
      m2 = subanalyze_topology(child, other)
      if m2:
        m = cl.mrca(m, cl.get_parent(m2))
    if m:
      cross_mrcas[tnu] = m
      return m
  return None

# Returns an articulation, as long as tnu is a fringe ancestor

def topological_match(tnu):
  partner = cross_mrcas.get(tnu)
  if partner:
    back = cross_mrcas.get(partner)
    if tnu == back:
      return art.art(tnu, partner, rel.eq, 0, "by topology")
    else:
      return art.art(tnu, partner, rel.lt, 100, "by topology")
  else:
    return None

# Set up 1-1 correspondence between matching fringe nodes

def analyze_fringes(A, B):
  analyze_fringe(A, B)
  analyze_fringe(B, A)

fringe_status = {}

def analyze_fringe(checklist, other):
  for root in cl.get_roots(checklist):
    subanalyze_fringe(root, other)

def subanalyze_fringe(tnu, other):
  children = cl.get_children(tnu)
  if len(children) > 0:
    found_match = False
    for child in children:
      if subanalyze_fringe(child, other):
        found_match = True
    if found_match:
      fringe_status[tnu] = 2  # Ancestor of fringe
      return True
  matches = name_based_matches(tnu, other)
  if len(matches) > 0:
    fringe_status[tnu] = 1    # On fringe
    return True
  else:
    fringe_status[tnu] = 0  # Descendant of fringe
    return False

# Returns a single best articulation

def mutual_fringe_match(tnu, other):
  if fringe_status[tnu] == 1:
    arts = name_based_matches(tnu, other)
    if arts:
      partner = arts[0].cod
      if fringe_status[partner] == 1:
        arts2 = name_based_matches(partner, cl.get_checklist(tnu))
        if arts2 and arts2[0].cod == tnu:
          return arts[0]
  return None

# Name and synonym based matching

name_based_matches_cache = {}

def name_based_matches(node1, other):
  kept = name_based_matches_cache.get(node1, None)
  if kept == None:
    arts = []
    assert node1 + 0 == node1
    for art1 in synonym_matches(node1):
      for art2 in direct_name_matches(art1.cod, other):
        art12 = art.compose(art1, art2)
        for art3 in accepted_matches(art12.cod):
          arts.append(art.compose(art12, art3))
    arts = sorted(arts, key=lambda art: art.badness)
    kept = []
    seen = []
    for a in arts:
      if not a.cod in seen:
        kept.append(a)
        seen.append(a.cod)
    name_based_matches_cache[node1] = kept
  return kept

def synonym_matches(acc):
  return [art.art(acc, acc, rel.eq, 0, None)] + \
         [syn_articulation(syn, acc) for syn in cl.get_synonyms(acc)]

def syn_articulation(syn, node):
  (rel, badness, status) = analyze_badness(syn)
  return art.art(node, syn, rel, badness, status)

def direct_name_matches(node, other):
  name = cl.get_name(node)
  candidates = cl.get_tnus_with_value(other, cl.canonical_name_field, name)
  return [art.art(node, candidate, rel.eq, 0, "bridge")
            for candidate in candidates]

def accepted_matches(node):
  if cl.is_accepted(node):
    return [art.art(node, node, rel.eq, 0, None)]
  else:
    return [art.reverse(syn_articulation(node, cl.get_accepted(node)))]

# For assigning priorities to synonyms

def analyze_badness(syn):
  status = syn_status(syn)
  rel_and_badness = badnesses.get(status)
  if rel_and_badness is None:
    # unrecognized nomenclatural status, assume it's bad
    rel_and_badness = badnesses.get("synonym")
  (rel, badness) = rel_and_badness
  return (rel, badness, status)

def syn_status(synonym):
  return cl.get_value(synonym, cl.nomenclatural_status_field) or \
         cl.get_value(synonym, cl.taxonomic_status_field) or \
         "synonym"

# Name / match classes, best to worst

badnesses = {
  "identical names": (rel.eq, 1),        # namestrings are the same
  "authority": (rel.eq, 1.5),
  "scientific name": (rel.eq, 2),        # (actually canonical) exactly one per node
  "equivalent name": (rel.eq, 3),        # synonym but not nomenclaturally
  "misspelling": (rel.eq, 3.8),
  "genbank synonym": (rel.eq, 4),        # at most one per node; first among equals
  "anamorph": (rel.eq, 4.1),
  "genbank anamorph": (rel.eq, 4.15),    # at most one per node
  "teleomorph": (rel.eq, 4.2),
  "unpublished name": (rel.eq, 4.5),    # non-code synonym
  "id": (rel.eq, 4.7),
  "merged id": (rel.eq, 4.8),
  "acronym": (rel.eq, 4.9),

  # above here: equivalence implied. below here: acc>=syn implied.
  # except in the case if 'in-part' which is acc<syn.

  "synonym": (rel.approx, 5),
  "misnomer": (rel.approx, 5.5),
  "includes": (rel.gt, 6),
  "in-part": (rel.lt, 6.5),      # this node is part of a polyphyly
  "type material": (rel.ge, 7),
  "blast name": (rel.approx, 8),             # large well-known taxa
  "genbank common name": (rel.approx, 9),    # at most one per node
  "genbank acronym": (rel.approx, 9.2),      # at most one per node
  "common name": (rel.approx, 10),
}

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
