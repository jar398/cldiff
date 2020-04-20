
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

  analyze_grafts(A, B)

  report(A, B, out)

def report(A, B, outpath):
  grafts = cl.invert_dict(graft_points)
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)
    for root in cl.get_roots(A):
      subreport(root, A, B, grafts, writer, "")

def subreport(node, A, B, grafts, writer, indent):
  ar = best_match(node, B)
  if ar:
    if ar.relation == rel.eq or ar.relation == rel.similar:
      if cl.get_name(ar.dom) == cl.get_name(ar.cod):
        tag = "NO CHANGE"
      else:
        tag = "RENAME"
    else:
      tag = "REORG"
    writer.writerow([indent + tag,
                     cl.get_unique(node),
                     ar.relation.name,
                     cl.get_unique(ar.cod),
                     ar.comment])
  else:
    ar2 = best_match(cl.get_parent(node), B)
    writer.writerow([indent + "REMOVE",
                     cl.get_unique(node),
                     "<",
                     cl.get_unique(ar2.cod),
                     "A has no match in B"])
  for child in cl.get_children(node):
    subreport(child, A, B, grafts, writer, indent + "  ")
  for child in grafts.get(node, ()):
    assert child > 0
    if not best_match(child, A):
      writer.writerow([indent + "ADD",
                       cl.get_unique(node),
                       ">",
                       cl.get_unique(child),
                       "B has no match in A"])

def best_match(node, B):
  assert node > 0
  return mutual_fringe_match(node, B) or \
         best_topological_match(node, B)

graft_points = {}

def analyze_grafts(A, B):
  def process(tnu):
    parent = cl.get_parent(tnu)
    if parent:
      ar = best_match(parent, A)
      parent_in_A = ar.cod
    for child in cl.get_children(tnu):
      other = best_match(child, A)
      if other == None and parent:
        graft_points[child] = parent_in_A
      else:
        process(child)
  for root in cl.get_roots(B):
    process(root)

# Returns an articulation, as long as tnu is a fringe ancestor

def best_topological_match(tnu, other):
  matches = topological_matches(tnu, other)
  if matches:
    return matches[0]
  else:
    return None

def topological_matches(tnu, other):
  ar = articulate(tnu, other)
  if ar:
    if ar.relation != rel.eq: print("gotcha", cl.get_name(tnu))
    options = topological_options(tnu, ar.cod, other) # in other
    if len(options) > 1:
      print("%s options for %s" % (len(options), cl.get_unique(tnu)))
    arts = [inform_option_by_name(tnu, option, other, ar.relation, 0)
            for option in options]
    return sorted(arts, key=lambda ar: ar.badness)
  else:
    return []

# Determine how tnu compares to its cross-mrca

def articulate(tnu, other):
  partner = cross_mrca_or_fringe(tnu, other)
  if partner:
    return art.art(tnu, partner, rel.eq, 0.5, "on fringe")
  partner = cross_mrcas.get(tnu)
  if not partner:
    return None
  back = cross_mrca_or_fringe(partner, cl.get_checklist(tnu))
  if cl.are_disjoint(tnu, back):
    relation = rel.disjoint
  else:
    if cl.mrca(tnu, back) == tnu:
      relation = rel.eq
    else:
      relation = rel.lt
    for child in cl.get_children(partner):
      back = cross_mrca_or_fringe(child, other)
      if back and mrca(back, tnu) == tnu and cross_disjoint(tnu, partner):
        relation = rel.conflict
        break
  return art.art(tnu, partner, relation, 0.5, "fringe comparison")

def cross_disjoint(tnu, partner):
  assert tnu > 0
  back = cross_mrca_or_fringe(partner, cl.get_checklist(partner))
  if not back: return True
  if cl.are_disjoint(tnu, back):
    return True
  for child in partner.children:
    if not cross_disjoint(tnu, child):
      return False
  return True

def cross_mrca_or_fringe(tnu, other):
  assert tnu > 0
  ar = mutual_fringe_match(tnu, other)
  if ar:
    return ar.cod
  else:
    return cross_mrcas.get(tnu)

def mutual_mrca_partner(tnu, other):
  partner = cross_mrcas.get(tnu)     # in other checklist
  if partner:
    back = cross_mrcas.get(partner)  # in tnu's checklist
    # tnu could be an ancestor of back
    return (partner, cl.mrca(back, tnu) == tnu)
  return None

# What are possible matches of tnu in the other checklist?

def topological_options(tnu, partner, other):
  assert cl.get_checklist(tnu) != other
  assert cl.get_checklist(partner) == other
 
  options = []
  option = partner
  assert cl.get_checklist(option) == other
  # Scan lineage of partner until a larger group is found
  checklist = cl.get_checklist(tnu)
  while option:
    if reflection(option, checklist) != partner:
      break
    options.append(option)
    option = cl.get_parent(option)
  return options

# Some near descend of tnu via monotypics.

def reflection(tnu, other):
  assert tnu > 0
  partner = cross_mrca_or_fringe(tnu, other)
  if partner:
    assert partner > 0
    return cross_mrca_or_fringe(partner, cl.get_checklist(tnu)) or tnu
  else:
    return tnu

def inform_option_by_name(tnu, option, other, relation, extra):
  assert cl.get_checklist(option) == other
  matches = name_based_matches(tnu, other) # these come sorted
  if len(matches) != 1:
    print("matches for %s = %s" % (cl.get_name(tnu), [m.cod for m in matches]))
  if len(matches) > 0:
    # tnu has offered many other-nodes with various names.
    # see if one of them is the one we have in hand
    for match in matches:    # match : tnu.checklist -> other
      assert match.dom == tnu
      assert cl.get_checklist(match.cod) == other
      if match.cod == option:
        return art.art(tnu, option, rel.eq,
                       match.badness + extra,
                       "topology + name match")
  if cl.get_rank(tnu) == cl.get_rank(option):
    return art.art(tnu, option, rel.eq,
                        100 + extra,
                        "topology + rank match")
  return art.art(tnu, option, relation,
                 200 + extra,
                 "topology")

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
        m = cl.mrca(m, m2)
    if m:
      cross_mrcas[tnu] = m
      return m
  return None

# Set up 1-1 correspondence between matching fringe nodes

# Returns a single best mutual fringe articulation

def mutual_fringe_match(tnu, other):
  assert tnu > 0
  ar = best_fringe_match(tnu, other)
  if ar:
    art2 = best_fringe_match(ar.cod, cl.get_checklist(tnu))
    if art2 and art2.cod == tnu:
      return ar
  return None

def best_fringe_match(tnu, other):
  assert tnu > 0
  if fringe_status[tnu] == 1:
    matches = name_based_matches(tnu, other)
    if matches:
      return matches[0]
  return None

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
    kept = processed_matches(arts)
    name_based_matches_cache[node1] = kept
  return kept

def processed_matches(arts):
  arts = sorted(arts, key=lambda art: art.badness)
  kept = []
  seen = []
  for a in arts:
    if not a.cod in seen:
      kept.append(a)
      seen.append(a.cod)
  return kept

def synonym_matches(acc):
  return [art.art(acc, acc, rel.eq, 0, None)] + \
         [syn_articulation(syn, acc) for syn in cl.get_synonyms(acc)]

def syn_articulation(syn, node):
  (rel, badness, status) = analyze_badness(syn)
  return art.art(node, syn, rel, badness, status)

def direct_name_matches(node, other):
  assert node > 0
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
  "revised under": (rel.lt, 1),          # hmm
  "identical names": (rel.eq, 1.2),      # namestrings are the same
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

  "synonym": (rel.similar, 5),
  "misnomer": (rel.similar, 5.5),
  "includes": (rel.gt, 6),
  "in-part": (rel.lt, 6.5),      # this node is part of a polyphyly
  "type material": (rel.ge, 7),
  "blast name": (rel.similar, 8),             # large well-known taxa
  "genbank common name": (rel.similar, 9),    # at most one per node
  "genbank acronym": (rel.similar, 9.2),      # at most one per node
  "common name": (rel.similar, 10),
}

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
