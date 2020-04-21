
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
  print ("number of fringe matches:", len(is_fringe))

  analyze_cross_mrcas(A, B)
  print ("number of cross-mrcas:", len(cross_mrcas))

  analyze_grafts(A, B)

  report(A, B, out)

def report(A, B, outpath):
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)
    for root in cl.get_roots(A):
      subreport(root, B, writer, "")

reported = {}

def subreport(node, B, writer, indent):
  assert cl.get_checklist(node) != B
  match = best_match(node, B)      # cod is accepted
  if match:
    report_on_match(match, writer, indent)
    for sub in get_subarticulations(node):
      subreport(sub.cod, B, writer, indent + "  ")
  elif cl.is_accepted(node):
    writer.writerow([indent + "REMOVE",
                     cl.get_unique(node),
                     "",
                     "",
                     "no match in B"])

def get_subarticulations(node):
  return (child_articulations(node) +
          graft_articulations(node))

def report_on_match(match, writer, indent):
  if rel.is_variant(match.relation, rel.eq):
    if cl.get_name(match.dom) == cl.get_name(match.cod):
      if cl.get_tnu_id(match.dom) == cl.get_tnu_id(match.cod):
        tag = "NO CHANGE"
      else:
        tag = "CHANGE ID"
    else:
      tag = "RENAME"
  elif rel.is_variant(match.relation, rel.lt):
    tag = "REMOVE"
  elif rel.is_variant(match.relation, rel.gt):
    tag = "ADD"
  elif rel.is_variant(match.relation, rel.disjoint):
    tag = "MOVE"    # shouldn't happen ...
  else:
    tag = "REORGANIZE"
  writer.writerow([indent + tag,
                   cl.get_unique(match.dom),
                   rel.variant(match.relation, 0).name,
                   cl.get_unique(match.cod),
                   art.get_comment(match)])

# For any node in A, return the most appropriate related node in B

def best_match(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  return best_topological_match(node, other) or \
    best_name_based_match(node, other)

def best_topological_match(node, other):
  matches = topological_matches(node, other)
  if len(matches) > 0:
    return matches[0]
  return None

# Grafts

def graft_articulations(node):
  return [art.art(node, graft, rel.lt)
          for graft in graft_points.get(node, ())]

grafts = {}
graft_points = {}

def analyze_grafts(A, B):
  global grafts
  def process(tnu, sup):
    if sup:
      sup_match = best_match(sup, A)
    for inf_art in child_articulations(tnu):
      inf = inf_art.cod
      if sup:
        other = best_match(inf, A)
        if other == None:
          graft_points[inf] = art.compose(sup_match,
                                          rel.reverse(inf_art))
      else:
        process(inf, tnu)
  for root in cl.get_roots(B):
    process(root, None)
  grafts = cl.invert_dict(graft_points)

# ---------- TOPOLOGY

# Sorted

def topological_matches(tnu, other):
  assert cl.get_checklist(tnu) != other

  match = compare_fringes(tnu, other)    # Single topo match
  if not match: return []
  matches = [match]

  if rel.is_variant(match.relation, topo_eq):
    # Scan upwards looking for nodes whose cross_mrca is us... not
    # quite right
    scan = match.cod    # tnu -> something
    while True:
      scan = cl.get_superior(scan)
      if scan == None: break
      if cross_mrcas.get(scan) != tnu: break
      matches.append(art.art(tnu, scan, topo_eq))
  return sorted([score_topo_match(match) for match in matches],
                key=lambda match: match.relation.badness)

# Determine how tnu compares to its cross-mrca

topology_badness = 1000
topo_eq       = rel.variant(rel.eq, topology_badness, "fringe=")
topo_lt       = rel.variant(rel.lt, topology_badness, "fringe<", "fringe>")
topo_gt       = rel.reverse(topo_lt)
topo_conflict = rel.variant(rel.conflict, topology_badness, "fringe-conflict")
topo_disjoint = rel.variant(rel.disjoint, topology_badness, "fringe-disjoint")

def compare_fringes(tnu, other):
  assert tnu > 0
  match = direct_fringe_match(tnu, other)
  if match: return match
  partner = cross_mrcas.get(tnu)    # another TNU I think
  if not partner:
    return None
  here = cl.get_checklist(tnu)
  back = cross_mrca_or_fringe(partner, here)
  if not back:
    relation = topo_disjoint
  elif cl.are_disjoint(tnu, back):
    relation = topo_disjoint
  else:
    if cl.mrca(tnu, back) == tnu:
      relation = topo_eq
    else:
      relation = topo_lt
    for sub in cl.get_inferiors(partner):
      back = cross_mrca_or_fringe(sub, here)
      if back:
        assert cl.get_checklist(tnu) == cl.get_checklist(back)
        if cl.mrca(tnu, back) == tnu and cross_disjoint(tnu, partner):
          relation = topo_conflict
          break
  return art.art(tnu, partner, relation)

def direct_fringe_match(tnu, other):
  matches = direct_fringe_matches(tnu, other)
  if len(matches) == 1:
    return matches[0]
  else:
    return None

def cross_disjoint(tnu, partner):
  assert tnu > 0
  assert partner > 0
  assert cl.get_checklist(tnu) != cl.get_checklist(partner)
  back = cross_mrca_or_fringe(partner, cl.get_checklist(tnu))
  if not back: return True
  assert back > 0
  assert cl.get_checklist(back) == cl.get_checklist(tnu)
  if cl.are_disjoint(tnu, back):
    return True
  for inf in cl.get_inferiors(partner):
    assert inf > 0
    if not cross_disjoint(tnu, inf):
      return False
  return True

# returns new articulation with same domain and codomain

def score_topo_match(match):
  tnu = match.dom
  option = match.cod
  other = cl.get_checklist(option)

  relation = match.relation
  
  if cl.get_rank(tnu) == cl.get_rank(option):
    relation = rel_fringe_and_rank
  
  matches = name_based_matches(tnu, other) # these come sorted
  # tnu has offered many other-nodes with various names.
  # see if one of them is the one we have in hand
  for match in matches:    # match : tnu.checklist -> other
    assert match.dom == tnu
    assert cl.get_checklist(match.cod) == other
    if match.cod == option:
      relation = rel_fringe_and_name
      break

  return art.art(tnu, option, relation)

rel_fringe_and_name = rel.variant(rel.eq, 2, "=fringe + =name")
rel_fringe_and_rank = rel.variant(rel.eq, 3, "=fringe + =rank")

#    assert match.dom == tnu
#    assert cl.get_checklist(match.cod) == other

# ---------- Cross-MRCAs

# Returns a tnu - we never care about the reason for the match here

def cross_mrca_or_fringe(tnu, other):
  assert tnu > 0
  match = mutual_fringe_match(tnu, other)
  if match:
    return match.cod
  else:
    return cross_mrcas.get(tnu)

def analyze_cross_mrcas(A, B):
  analyze_topology(A, B)
  analyze_topology(B, A)

def analyze_topology(checklist, other):
  for root in cl.get_roots(checklist):
    subanalyze_topology(root, other)

cross_mrcas = {}

def subanalyze_topology(tnu, other):
  inferiors = cl.get_inferiors(tnu)
  art = mutual_fringe_match(tnu, other)
  if art:
    return art.cod
  m = None
  for inf in inferiors:
    m2 = subanalyze_topology(inf, other)
    if m2:
      m = cl.mrca(m, m2)
  if m:
    assert cl.get_checklist(tnu) != cl.get_checklist(m)
    cross_mrcas[tnu] = m
    return m
  return None

# ---------- Fringe

# Returns a single best mutual fringe match

def mutual_fringe_match(tnu, other):
  assert tnu > 0
  match = best_fringe_match(tnu, other)
  if match:
    art2 = best_fringe_match(match.cod, cl.get_checklist(tnu))
    if art2 and art2.cod == tnu:
      return match
  return None

def best_fringe_match(tnu, other):
  assert tnu > 0
  if is_fringe.get(tnu):
    matches = direct_fringe_matches(tnu, other)
    if matches:
      if len(matches) == 1:
        return matches[0]
      else:
        print("Ambiguous", cl.get_name(tnu))
  return None

def direct_fringe_matches(tnu, other):
  return [match for match in direct_matches(tnu, other)
          if is_fringe.get(tnu)]

# ---------- Fringe determination

def analyze_fringes(A, B):
  print('A')
  print(cl.get_roots(A))
  analyze_fringe(A, B)
  print('B')
  analyze_fringe(B, A)

is_fringe = {}

def analyze_fringe(checklist, other):
  for root in cl.get_roots(checklist):
    subanalyze_fringe(root, other)

def subanalyze_fringe(tnu, other):
  found_match = False
  for inf in cl.get_inferiors(tnu):
    if subanalyze_fringe(inf, other):
      found_match = True
  if found_match:
    return True
  name = cl.get_name(tnu)
  partners = direct_matches(tnu, other)
  if len(partners) == 1:
    is_fringe[tnu] = True
    return True
  elif len(partners) > 1:
    return False

# ---------- Matches based on name and synonym

# Three components: synonym-or-self o direct o synonym-of-or-self
# from_accepted_articulations o direct_matches o [to_accepted_articulation]

def best_name_based_match(accepted, other):
  m1 = name_based_matches(accepted, other)
  if m1:
    winner = m1[0]
    m2 = name_based_matches(winner.cod, cl.get_checklist(accepted))
    if len(m2) > 0:
      if len(m2) > 1:
        if m2[0] == m1[0]:
          return art.compose(winner,
                             art.art(winner.cod, winner.cod, rel_best_eq))
        else:
          return art.compose(winner,
                             art.art(winner.cod, winner.cod, rel_alternate_eq))
      else:
        return winner
  return None

rel_best_eq = rel.variant(rel.eq, 0.1, "best=")
rel_alternate_eq = rel.variant(rel.eq, 0.2, "alternate=")


name_based_matches_cache = {}

def name_based_matches(tnu, other):
  probe = name_based_matches_cache.get(tnu)
  if probe: return probe

  matches = [art.compose(a, to_accepted_match(direct))
             for a in from_accepted_articulations(tnu)
             for direct in direct_matches(a.cod, other)]
  matches = prune_matches(matches)

  name_based_matches_cache[tnu] = matches
  return matches

# not used (yet?)
def accepted_direct_matches(tnu, other):
  return [to_accepted_match(direct)
          for direct in direct_matches(tnu, other)]

def to_accepted_match(m):
  if cl.is_accepted(m.cod):
    return m
  else:
    return art.compose(m, accepted_articulation(m.cod))

# Direct matches by name (no synonym following)

def direct_matches(node, other):
  assert node > 0
  name = cl.get_name(node)
  candidates = cl.get_tnus_with_value(other, cl.canonical_name_field, name)
  matches = [art.art(node, candidate, same_namestring) 
             for candidate in candidates]
  partner = cl.get_tnu_with_id(other, cl.get_tnu_id(node))
  if partner and not partner in candidates:
    matches = [art.art(node, partner, same_id)] + matches
  return matches

same_namestring = rel.variant(rel.eq, 2, "namestring=", "namestring=")
same_id = rel.variant(rel.eq, 2.1, "id=", "id=")

# ---------- Within-checklist articulations

# Synonym-or-self = to-accepted

def from_accepted_articulations(node):
  return [art.identity(node)] + synonym_articulations(node)

# Superior/inferior

def superior_articulation(tnu):
  return parent_articulation(tnu) or accepted_articulation(tnu)

def inferior_articulations(tnu):
  return child_articulations(tnu) + synonym_articulations(tnu)

# Parent/child

def parent_articulation(child):
  parent = cl.get_parent(child)
  if parent:
    return art.art(child, parent, rel.parent_child)
  else:
    return None

def child_articulations(node):
  return [art.reverse(parent_articulation(child))
          for child in cl.get_children(node)]

# Accepted/synonym

def accepted_articulation(syn):   # goes from synonym to accepted
  assert syn > 0
  accepted = cl.get_accepted(syn)
  if accepted:
    return art.art(syn, accepted, synonym_relation(syn_status(syn)))
  else:
    print("Shouldn't happen", cl.get_name(syn))
    return None

def syn_status(synonym):
  return cl.get_value(synonym, cl.nomenclatural_status_field) or \
         cl.get_value(synonym, cl.taxonomic_status_field) or \
         "synonym"

def synonym_articulations(tnu):
  if cl.is_accepted(tnu):
    return prune_matches([art.reverse(accepted_articulation(syn))
                          for syn in cl.get_synonyms(tnu)])
  else:
    return []

# ---------- Pruning

# Reduce a set of articulations so that all the codomains are
# different

def prune_matches(arts):
  arts = sorted(arts, key=lambda art: art.relation.badness)
  kept = []
  seen = []
  for a in arts:
    good = a.cod
    if not a.cod in seen:
      kept.append(a)
      seen.append(good)
  return kept

# ---------- NCBI nomenclatural statuses

def synonym_relation(status):
  return badnesses[status]

badnesses = {}
badness = 1

def declare_badnesses():
  def b(revname, re, name = None):
    assert re
    global badness
    name = name or (revname + " of")
    badnesses[revname] = rel.variant(re, badness, name, revname)
    badness += 1

  b("authority", rel.eq)
  b("scientific name", rel.eq)        # (actually canonical) exactly one per node
  b("equivalent name", rel.eq)        # synonym but not nomenclaturally
  b("misspelling", rel.eq)
  b("genbank synonym", rel.eq)        # at most one per node; first among equals
  b("anamorph", rel.eq)
  b("genbank anamorph", rel.eq)    # at most one per node
  b("teleomorph", rel.eq)
  b("unpublished name", rel.eq)    # non-code synonym
  b("merged id", rel.eq)
  b("acronym", rel.eq)

  # above here: equivalence implied. below here: acc>=syn implied.
  # except in the case if 'in-part' which is acc<syn.

  b("synonym", rel.eq)
  b("misnomer", rel.eq)
  b("includes", rel.gt, "included in")
  b("in-part", rel.lt, "part of")      # this node is part of a polyphyly
  b("type material", rel.eq)
  b("blast name", rel.eq)             # large well-known taxa
  b("genbank common name", rel.eq)    # at most one per node
  b("genbank acronym", rel.eq)      # at most one per node
  b("common name", rel.eq)

declare_badnesses()

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
