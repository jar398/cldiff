
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
  print ("id_match_count:", id_match_count)
  print ("tnu_count:", A.tnu_count())

  analyze_cross_mrcas(A, B)
  print ("number of cross-mrcas:", len(cross_mrcas))

  assign_matches(B, A)
  print ("number of besties:", len(anti_best))

  analyze_unmatched(A, B)
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
  multiple = report_on_matches(node, B, writer, indent)
  def for_seq(node):
    b = best_partner(node, B)
    return b.cod if b else node
  def sort_key(triple):
    (B_node, which, arg) = triple
    return cl.get_sequence_number(B_node)
  agenda = \
    [(for_seq(child), 0, child) for child in cl.get_children(node)] + \
    [(option.cod, 1, option) for option in multiple] + \
    [(B_node, 2, B_node) for B_node in get_graftees(node)]
  indent = indent + "__"
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

def report_on_matches(node, B, writer, indent):
  matches = best_partners(node, B)      # cod is accepted
  if len(matches) == 0:
    writer.writerow([indent + "REMOVE",
                     cl.get_unique(node),
                     "",
                     "",
                     "%s B nodes match this A node" % len(matches)])
    return []
  elif len(matches) == 1:
    report_on_match(matches[0], False, writer, indent)
    return []
  else:
    writer.writerow([indent + "MULTIPLE",
                     cl.get_unique(node),
                     "?",
                     "",
                     "%s B nodes match this A node" % len(matches)])
    return matches

def report_on_match(match, 
                    splitp, writer, indent):
  A_unique = "" if splitp else cl.get_unique(match.dom)
  tag = tag_for_match(match, splitp)
  writer.writerow([indent + tag,
                   A_unique,
                   rel.rcc5_name(match.relation),
                   cl.get_unique(match.cod),
                   art.get_comment(match)])

def tag_for_match(match, splitp):
  tag = "?"
  if rel.is_variant(match.relation, rel.eq):
    if splitp: tag = "OPTION"
    elif parent_changed(match):
      tag = "MOVE"
    else:
      if cl.get_name(match.dom) == cl.get_name(match.cod):
        if id_match_count * 2 < cl.get_checklist(match.dom).tnu_count() or \
           cl.get_tnu_id(match.dom) == cl.get_tnu_id(match.cod):
          tag = "NO CHANGE"
        else:
          tag = "CHANGE ID"
      else:
        tag = "RENAME"
  elif rel.is_variant(match.relation, rel.lt):
    tag = "ELIDE"
  elif rel.is_variant(match.relation, rel.gt):
    tag = "INSERT"
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
  other = cl.get_checklist(coparent)
  match = best_match(parent, other)
  if not match: return True
  return match.cod != coparent

# Partners list for reporting (A->B articulations)

def best_partner(tnu, other):  # for children sort order
  return choose_least_bad(best_partners(tnu, other))

def best_partners(tnu, other):
  matches = anti_best.get(tnu) or []
  if len(matches) == 0:
    investigate = best_match(tnu, other)
    if investigate: matches = [investigate]
  return matches

# Fill the cache

def assign_matches(here, other):
  global anti_best
  best_match_map = {}
  def process(tnu):
    m = best_match(tnu, other)  # Fill the cache
    if m: best_match_map[tnu] = m    # here -> other
    for child in cl.get_children(tnu):
      process(child)
  for root in cl.get_roots(here):
    process(root)
  anti_best = invert_dict_by_cod(best_match_map)

def invert_dict_by_cod(d):
  inv = {}
  for (node, ar) in d.items():
    if ar:
      ar = art.reverse(ar)
      if ar.dom in inv:
        inv[ar.dom].append(ar)
      else:
        inv[ar.dom] = [ar]
  return inv

# ---------- UNMATCHED

# Unmatched
# Find nodes in B that are not mutually matched to nodes in A

# A B_node that is not the best_match of any A_node

def analyze_unmatched(A, B):
  global grafts
  graft_points = {}
  def process(B_tnu):
    if not best_match(B_tnu, A):
      point = get_graft_point(B_tnu, A)    # in A
      if point:
        graft_points[B_tnu] = point    # in A
    for child in cl.get_children(B_tnu):
      process(child)
  for root in cl.get_roots(B):
    process(root)
  grafts = cl.invert_dict(graft_points)

def get_graftees(A_node):
  return (grafts.get(A_node) or [])

def get_graft_point(B_tnu, A = None):
  B_parent = cl.get_parent(B_tnu)
  if B_parent:
    B_parent_match = best_match(B_parent, A)
    if B_parent_match:
      return B_parent_match.cod
  return None

# ---------- One-sided best match

best_match_cache = {}

def best_match(node, other = None):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in best_match_cache:
    return best_match_cache[node]

  assert other
  match = choose_least_bad(good_matches(node, other))
    
  best_match_cache[node] = match
  return match

good_matches_cache = {}

def good_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in good_matches_cache: return good_matches_cache[node]

  topos = topological_matches(node, other)
  nameys = name_based_matches(node, other)
  matches = [art.conjoin(topo, namey) for topo in topos for namey in nameys]
  matches = [match for match in matches if match]
  matches = matches or topos or nameys
  # matches = [score_topo_match(topo, nameys) for topo in topos] + nameys
  matches = prune_matches([to_accepted_match(match) for match in matches])

  good_matches_cache[node] = matches
  return matches

# Returns new articulation with same domain and codomain

def score_topo_match(match, nameys):
  for namey in nameys:
    assert is_match(namey)
    if namey.cod == match.cod:
      return art.art(match.dom, match.cod, rel.fringe_and_name)
  return match

# ---------- TOPOLOGY

# Sorted

def topological_matches(tnu, other):
  assert cl.get_checklist(tnu) != other

  match = compare_fringes(tnu, other)    # Single topo match
  if not match: return []
  matches = [match]

  if rel.is_variant(match.relation, rel.topo_eq):
    # Scan upwards looking for nodes whose cross_mrca is us... not
    # quite right
    scan = match.cod    # tnu -> something
    while True:
      scan = cl.get_superior(scan)
      if scan == None: break
      if cross_mrcas.get(scan) != tnu: break
      matches.append(bridge(tnu, scan, rel.topo_eq))
  matches.reverse()    # hmm. for choose
  return matches

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
    relation = rel.topo_disjoint
  elif cl.are_disjoint(tnu, back):
    relation = rel.topo_disjoint
  else:
    if cl.mrca(tnu, back) == tnu:
      relation = rel.topo_eq
    else:
      relation = rel.topo_lt
    for sub in cl.get_inferiors(partner):
      back = cross_mrca_or_fringe(sub, here)
      if back:
        assert cl.get_checklist(tnu) == cl.get_checklist(back)
        if cl.mrca(tnu, back) == tnu and cross_disjoint(tnu, partner):
          relation = rel.topo_conflict
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

# One-sided fringe determination... 

id_match_count = 0

def analyze_fringe(checklist, other):
  def subanalyze_fringe(tnu, other):
    global id_match_count
    found_match = False
    for inf in cl.get_inferiors(tnu):
      if subanalyze_fringe(inf, other):
        found_match = True
    if found_match:
      return True
    partners = direct_matches(tnu, other)
    if len(partners) == 1:
      is_fringe[tnu] = True
      if cl.get_tnu_id(tnu) == cl.get_tnu_id(partners[0].cod):
        id_match_count += 1
      return True
    elif len(partners) > 1:
      return False
  for root in cl.get_roots(checklist):
    subanalyze_fringe(root, other)

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

def to_accepted_match(m):
  if cl.is_synonym(m.cod):
    return art.compose(m, accepted_articulation(m.cod))
  else:
    return m

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
                    rel.same_name_and_id if hit == id_hit \
                    else rel.same_name)
             for hit in hits]
  assert is_matches(matches) and True
  if id_hit and not id_hit in hits:
    matches = [bridge(node, id_hit, rel.same_id)] + matches
  assert is_matches(matches)
  return matches

# ---------- Within-checklist articulations

# Synonym-or-self = to-accepted

def from_accepted_articulations(node):
  return [art.identity(node)] + synonym_articulations(node)

# Superior/inferior

#def superior_articulation(tnu):
#  return parent_articulation(tnu) or accepted_articulation(tnu)

#def inferior_articulations(tnu):
#  return child_articulations(tnu) + synonym_articulations(tnu)

# Parent/child

# def parent_articulation(child):
#   parent = cl.get_parent(child)
#   if parent:
#     return art.art(child, parent, rel.child_of)
#   else:
#     return None

#def child_articulations(node):
#  return [art.reverse(parent_articulation(child))
#          for child in cl.get_children(node)]

# Accepted/synonym

def accepted_articulation(syn):   # goes from synonym to accepted
  assert syn > 0
  accepted = cl.get_accepted(syn)
  if accepted:
    return art.art(syn, accepted, rel.synonym_relation(syn_status(syn)))
  else:
    print("Shouldn't happen", cl.get_name(syn))
    return None

def syn_status(synonym):
  return cl.get_value(synonym, cl.nomenclatural_status_field) or \
         cl.get_value(synonym, cl.taxonomic_status_field) or \
         "synonym"

def synonym_articulations(tnu):
  if cl.is_synonym(tnu):
    return []
  else:
    return prune_matches([art.reverse(accepted_articulation(syn))
                          for syn in cl.get_synonyms(tnu)])

# ---------- Pruning

# Reduce a set of articulations so that all the codomains are
# different

def prune_matches(arts):
  if len(arts) <= 1: return arts
  assert arts[0].dom
  arts = sorted(arts, key=prune_ordering)
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
  print("** Multiple least-bad matches. Need to look at other criteria.")
  print("%s -> %s" %
        (cl.get_unique(arts[0].dom),
         [cl.get_unique(a.cod) for a in besties]))
  print(', '.join(["%o" % ar.relation.goodness for ar in besties]))
  return None

def prune_further(arts):
  if len(arts) == 0: return []
  arts = sorted(arts, key=prune_ordering)
  # if len(arts) == 1: return [arts[0]]
  key = prune_ordering(arts[0])
  return [a for a in arts
          if prune_ordering(a) == key]

def prune_ordering(ar):
  return (rel.sort_key(ar.relation), -cl.get_rank(ar.cod))

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

# ----

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
