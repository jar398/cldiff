
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

  assign_matches(A, B)
  print ("number of besties:", len(anti_best))

  analyze_unmatched(A, B)
  print ("number of joins:", len(joins))
  print ("number of grafts:", len(grafts))

  report(A, B, out)

def report(A, B, outpath):
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)
    for root in cl.get_roots(A):
      subreport(root, B, writer, "")

reported = {}

def subreport(node, B, writer, indent):
  A = cl.get_checklist(node)
  assert A != B
  match = best_match(node, B)      # cod is accepted
  if match:
    report_on_match(match, False, writer, indent)
    indent = indent + "  "
    def for_seq(node):
      b = best_match(node)
      return b.cod if b else node
    def sort_key(triple):
      (B_node, which, arg) = triple
      if isinstance(B_node, art.Articulation):
        print("Loser", which, cl.get_name(B_node.cod))
      return cl.get_sequence_number(B_node)
    agenda = \
      [(for_seq(child), 0, child) for child in cl.get_children(node)] + \
      [(match.cod, 1, match) for match in join_matches(node)] + \
      [(B_node, 2, B_node) for B_node in get_graftees(node)]
    for (B_node, which, arg) in \
       sorted(agenda, key=sort_key):
      if which == 0:
        subreport(arg, B, writer, indent)
      elif which == 1:
        report_on_match(arg, True, writer, indent)    # split
      elif which == 2:
        writer.writerow([indent + "GRAFT",
                         "",
                         ">",
                         cl.get_unique(arg),
                         "not a match target"])

"""
    for child in cl.get_children(node):
      subreport(child, B, writer, indent)
    for match in join_matches(node):
      report_on_match(match, True, writer, indent)
    for B_node in get_graftees(node):
      writer.writerow([indent + "GRAFT",
                       "",
                       ">",
                       cl.get_unique(B_node),
                       "not a match target"])
"""

def report_on_match(match, 
                    splitp, writer, indent):
  A_unique = "" if splitp else cl.get_unique(match.dom)
  tag = tag_for_match(match, splitp)
  writer.writerow([indent + tag,
                   A_unique,
                   rel.variant(match.relation, rel.goodness).name,
                   cl.get_unique(match.cod),
                   art.get_comment(match)])

def tag_for_match(match, splitp):
  tag = "?"
  if rel.is_variant(match.relation, rel.eq):
    if splitp: tag = "SPLIT"
    elif parent_changed(match):
      tag = "MOVE"
    else:
      if cl.get_name(match.dom) == cl.get_name(match.cod):
        if cl.get_tnu_id(match.dom) == cl.get_tnu_id(match.cod):
          tag = "NO CHANGE"
        else:
          tag = "CHANGE ID"
      else:
        tag = "RENAME"
  elif rel.is_variant(match.relation, rel.lt):
    tag = "INSERT" if splitp else "REMOVE"
  elif rel.is_variant(match.relation, rel.gt):
    tag = "GRAFT" if splitp else "ADD"
  elif rel.is_variant(match.relation, rel.conflict):
    tag = "REFORM" if splitp else "BREAK"
  elif rel.is_variant(match.relation, rel.disjoint):
    tag = "PEPPERONI"    # shouldn't happen ...
  return tag

def parent_changed(match):
  parent = cl.get_parent(match.dom)
  coparent = cl.get_parent(match.cod)
  if parent == None and coparent == None:
    return False
  if parent == None or coparent == None:
    return True
  match = best_match(parent)
  if not match: return True
  return match.cod != coparent

# Fill the cache

def assign_matches(A, B):
  global anti_best
  anti_best = {}
  def process(tnu):
    bestie = best_match(tnu, B)  # A -> B
    if bestie:
      anti_best[bestie.cod] = art.reverse(bestie)
    for child in cl.get_children(tnu):
      process(child)
    return bestie
  for root in cl.get_roots(A):
    process(root)

def anti_best_match(B_tnu, parent = None):
  return anti_best.get(B_tnu)

# ---------- UNMATCHED

# Unmatched
# Find nodes in B that are not mutually matched to nodes in A

def join_matches(A_node):
  A = cl.get_checklist(A_node)
  B_nodes = (joins.get(A_node) or [])
  return [join_match(B_node, A) for B_node in B_nodes]

# A B_node that is not the best_match of any A_node

def analyze_unmatched(A, B):
  global joins
  global grafts
  join_points = {}
  graft_points = {}
  def process(B_tnu):
    if anti_best_match(B_tnu) == None:
      match = best_match(B_tnu, A)  # B -> A
      if match:
        assert is_match(match)
        assert match.dom == B_tnu
        join_points[B_tnu] = match.cod    # in A
      else:
        point = get_graft_point(B_tnu, A)    # in A
        if point:
          graft_points[B_tnu] = point    # in A
    for child in cl.get_children(B_tnu):
      process(child)
  for root in cl.get_roots(B):
    process(root)
  joins = cl.invert_dict(join_points)
  grafts = cl.invert_dict(graft_points)

def join_match(B_tnu, A):
  assert cl.get_checklist(B_tnu) != A
  match = best_match(B_tnu, A)
  return art.reverse(art.compose(art.art(match.dom, match.dom, rel_alternative),
                                 match))

def get_graftees(A_node):
  return (grafts.get(A_node) or [])

def graft_match(B_tnu, A):
  A_parent = get_graft_point(B_tnu, A)
  if A_parent:
    return to_graft_match(A_parent, B_tnu)
  else:
    return None

def to_graft_match(A_parent, B_tnu):
  to_B_parent = parent_articulation(B_tnu)
  match = art.compose(to_B_parent, B_parent_match) # B -> A
  return art.compose(art.art(match.dom, match.dom, rel_graft),
                     match)

def get_graft_point(B_tnu, A = None):
  A_match = anti_best_match(B_tnu)  # B -> A
  B_parent = cl.get_parent(B_tnu)
  if B_parent:
    B_parent_match = anti_best_match(B_parent, A)
    if B_parent_match:
      return B_parent_match.cod
  return None

rel_alternative = rel.variant(rel.eq, 0.5, "split", "merge")
rel_graft = rel.variant(rel.eq, 0, "graft", "prune")

# ---------- One-sided best match

best_match_cache = {}

def best_match(node, other = None):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in best_match_cache:
    return best_match_cache[node]

  assert other

  match = choose_least_bad(best_matches(node, other))
    
  best_match_cache[node] = match
  return match

def best_matches(node, other):
  assert node > 0
  matches = prune_matches([to_accepted_match(match)
                           for match in topological_matches(node, other)])
  if matches: return matches
  return name_based_matches(node, other)

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
      matches.append(bridge(tnu, scan, topo_eq))
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
  match = best_fringe_match(tnu, other)
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
  return bridge(tnu, partner, relation)

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

  return bridge(tnu, option, relation)

rel_fringe_and_name = rel.variant(rel.eq, 2, "fringe= + name=")
rel_fringe_and_rank = rel.variant(rel.eq, 3, "fringe= + rank=")

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
    return choose_least_bad(direct_fringe_matches(tnu, other))
  return None

def direct_fringe_matches(tnu, other):
  if not is_fringe.get(tnu): return []
  d = direct_matches(tnu, other)
  assert is_matches(d)
  return [match for match in d
          if is_fringe.get(match.cod)]

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

def analyze_hits(checklist, other):
  def subanalyze_hits(tnu, other):
    found_match = False
    for inf in cl.get_inferiors(tnu):
      if subanalyze_hits(inf, other):
        found_match = True
    if found_match:
      return True
    partners = direct_matches(tnu, other)
    if len(partners) == 1:
      is_fringe[tnu] = True
      return True
    elif len(partners) > 1:
      return False
  for root in cl.get_roots(checklist):
    subanalyze_hits(root, other)


# ---------- Matches based on name and synonym

# Three components: synonym-or-self o direct o synonym-of-or-self
# from_accepted_articulations o direct_matches o [to_accepted_articulation]

name_based_matches_cache = {}

def name_based_matches(tnu, other):
  probe = name_based_matches_cache.get(tnu)
  if probe: return probe

  matches = [to_accepted_match(art.compose(a, direct))
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
# This is what open tree was so concerned about...

def direct_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  hits = cl.get_tnus_with_value(other,
                                cl.canonical_name_field,
                                cl.get_name(node))
  id_hit = cl.get_tnu_with_id(other, cl.get_tnu_id(node))

  matches = [bridge(node,
                    hit,
                    same_namestring_and_id if hit == id_hit \
                    else same_namestring)
             for hit in hits]
  assert is_matches(matches) and True
  if id_hit and not id_hit in hits:
    matches = [bridge(node, id_hit, same_id)] + matches
  assert is_matches(matches)
  return matches

same_namestring_and_id = rel.variant(rel.eq, 2,
                                     "name= + id=")
same_namestring = rel.variant(rel.eq, 2.1, "name=", "name=")
same_id = rel.variant(rel.eq, 2.2, "id=", "id=")

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
    return art.art(child, parent, rel.child)
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
  if len(arts) <= 1: return arts
  assert arts[0].dom
  arts = sorted(arts, key=lambda art: art.relation.badness)
  kept = []
  seen = []
  for a in arts:
    good = a.cod
    if not a.cod in seen:
      kept.append(a)
      seen.append(good)
  return kept

def choose_least_bad(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  besties = prune_further(arts)
  assert besties[0].dom
  if len(besties) == 1: return besties[0]
  print("** This shouldn't happen. Need to look at other criteria",
        [cl.get_unique(a.cod) for a in besties])
  return None

def prune_further(arts):
  if len(arts) == 0: return []
  arts = sorted(arts, key=prune_ordering)
  # if len(arts) == 1: return [arts[0]]
  key = prune_ordering(arts[0])
  return [a for a in arts
          if prune_ordering(a) == key]

def prune_ordering(ar):
  return (ar.relation.badness, cl.get_depth(ar.cod))

def is_match(ar):
  c1 = cl.get_checklist(ar.dom)
  c2 = cl.get_checklist(ar.cod)
  if c1 == c2:
    print("%s and %s both in %s" %
          (cl.get_name(c1),
           c2.get_name(c2),
           c1.prefix))
    return False
  return True

def is_matches(matches):
  if len(matches) == 0: return True
  return is_match(matches[0])

def bridge(dom, cod, re):
  assert cl.get_checklist(dom) != cl.get_checklist(cod)
  return art.art(dom, cod, re)

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
