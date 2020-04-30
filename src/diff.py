
import sys, csv
import argparse

import checklist as cl
import relation as rel
import articulation as art

def main(c1, c2, out):
  A = cl.read_checklist(c1, "A.")
  B = cl.read_checklist(c2, "B.")
  print ("TNU counts:", len(cl.get_all_tnus(A)), len(cl.get_all_tnus(B)))
  start(A, B)
  report(A, B, out)

def start(A, B):
  global particles
  global cross_mrcas

  particles = find_particles(A, B)
  print ("number of particles:", len(particles)>>1)

  cross_mrcas = analyze_cross_mrcas(A, B)
  print ("number of cross-mrcas:", len(cross_mrcas))

  assign_matches(B, A)
  print ("number of besties:", len(inverse_good_candidates))

  analyze_unmatched(A, B)
  print ("number of grafts:", len(grafts))

def alignment(A, B):
  start(A, B)

  # Unordered set of matches
  alignment = {}

  def subprepare(node):
    for match in good_candidates(node, B):      # cod is accepted
      alignment.append(match)
    for child in get_children(node):
      subprepare(child)
    for graftee in get_graftees(node):
      # Go from graftee's parent's match to graftee's parent, then graftee
      alignment.append([match_for_graftee(graftee) for graftee in graftee])

  for root in cl.get_roots(A):
    subprepare(root)

  return None

def report(A, B, outpath):
  with open(outpath, "w") as outfile:
    print ("Writing:", outpath)
    writer = csv.writer(outfile)
    sink = make_sink(writer, None)
    for root in cl.get_roots(A):
      subreport(root, B, sink, "")
    drain(sink)

def make_sink(writer, parent):
  return [False, [], writer, parent]

def subsink(sink):
  (changed, rows, writer, _) = sink
  return make_sink(writer, sink)

def drain(sink):
  (changed, rows, writer, parent) = sink
  if len(rows) == 0:
    pass
  elif changed:   # or len(rows) == 1:
    for row in rows:
      proclaim_row(row, parent)
    parent[0] = True
    sink[1] = []              # Ensure no double printing...
  else:
    indent = rows[0][0]
    proclaim_row((indent, "...", "", "", "",
                  "%s unchanged children" % len(rows)),
                 parent)

no_change_tags = ["NO CHANGE", "CHANGED ID", "..."]

def proclaim(sink, indent, tag, dom, re, cod, remark):
  if not tag in no_change_tags: sink[0] = True
  proclaim_row((indent, tag, dom, re, cod, remark), sink)

def proclaim_row(row, sink):
  (changed, rows, writer, parent) = sink
  if parent:
    rows.append(row)
  else:
    (indent, tag, dom, re, cod, remark) = row
    writer.writerow([indent + tag, dom, re, cod, remark])

# Report generation

def get_children(node):
  return [node
          for node in cl.get_children(node) + cl.get_synonyms(node)
          if not cl.get_accepted(node)]

def subreport(node, B, sink, indent):
  A = cl.get_checklist(node)
  assert A != B
  multiple = report_on_matches(node, B, sink, indent)
  sink = subsink(sink)
  def for_seq(node):
    b = good_candidate(node, B)
    return b.cod if b else node
  def sort_key(triple):
    (B_node, which, arg) = triple
    return cl.get_sequence_number(B_node)
  agenda = \
    [(for_seq(child), 0, child) for child in get_children(node)] + \
    [(option.cod, 1, option) for option in multiple] + \
    [(B_node, 2, B_node) for B_node in get_graftees(node)]
  indent = indent + "__"
  for (B_node, which, arg) in \
     sorted(agenda, key=sort_key):
    if which == 0:
      subreport(arg, B, sink, indent)
    elif which == 1:
      report_on_match(arg, True, sink, indent)    # split
    elif which == 2:
      proclaim(sink,
               indent, "GRAFT",
                       "",
                       "",
                       cl.get_unique(arg),
                       "")
  drain(sink)

def report_on_matches(node, B, sink, indent):
  matches = good_candidates(node, B)      # cod is accepted
  if len(matches) == 0:
    proclaim(sink, indent, "REMOVE",
                     cl.get_unique(node),
                     "",
                     "",
                     "%s B nodes match this A node" % len(matches))
    return []
  elif len(matches) == 1:
    report_on_match(matches[0], False, sink, indent)
    return []
  else:
    proclaim(sink, indent, "MULTIPLE",
                     cl.get_unique(node),
                     "?",
                     "",
                     "%s B nodes match this A node" % len(matches))
    return matches

def report_on_match(match, splitp, sink, indent):
  A_unique = "" if splitp else cl.get_unique(match.dom)
  tag = tag_for_match(match, splitp)
  proclaim(sink, indent, tag,
                   A_unique,
                   rel.rcc5_name(match.relation),
                   cl.get_unique(match.cod),
                   art.get_comment(match))

def tag_for_match(match, splitp):
  tag = "?"
  if rel.is_variant(match.relation, rel.eq):
    if splitp:
      if cl.get_accepted(match.cod):
        tag = "ADD SYNONYM"
      else:
        tag = "OPTION"
    elif parent_changed(match):
      tag = "MOVE"
    else:
      changes = []
      if cl.get_name(match.dom) != cl.get_name(match.cod):
        changes.append("NAME")
      if cl.get_nominal_rank(match.dom) != cl.get_nominal_rank(match.cod):
        changes.append("RANK")
      if cl.get_taxonomic_status(match.dom) != cl.get_taxonomic_status(match.cod):
        changes.append("STATUS")
      if len(changes) == 0 and cl.get_tnu_id(match.dom) != cl.get_tnu_id(match.cod):
        changes.append("ID")
      if changes:
        tag = "CHANGE " + " & ".join(changes)
      else:
        tag = "NO CHANGE"
  elif rel.is_variant(match.relation, rel.lt):
    tag = "ELIDE"
  elif rel.is_variant(match.relation, rel.gt):
    tag = "INSERT"
  elif rel.is_variant(match.relation, rel.conflict):
    tag = "REFORM" if splitp else "BREAK"
  elif rel.is_variant(match.relation, rel.disjoint):
    tag = "MUTEX"    # shouldn't happen ...
  return tag

def parent_changed(match):
  parent = cl.get_parent(match.dom)
  coparent = cl.get_parent(match.cod)
  if not parent and not coparent:
    return False
  if not parent or not coparent:
    return True
  other = cl.get_checklist(coparent)
  match = good_match(parent, other)
  if not match: return True
  return match.cod != coparent

# Partners list for reporting (A->B articulations)

def good_candidate(tnu, other):  # for children sort order
  return best_match(good_candidates(tnu, other))

def good_candidates(tnu, other):
  matches = inverse_good_candidates.get(tnu) or []
  if len(matches) == 0:
    investigate = good_match(tnu, other)
    if investigate: matches = [investigate]
  return [match for match in matches if not cl.get_accepted(match.cod)]


# Fill the cache

def assign_matches(here, other):
  global inverse_good_candidates
  good_match_map = {}
  def process(tnu):
    m = good_match(tnu, other)  # Fill the cache
    if m: good_match_map[tnu] = m    # here -> other
    for child in get_children(tnu):
      process(child)
  for root in cl.get_roots(here):
    process(root)
  inverse_good_candidates = invert_dict_by_cod(good_match_map)

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
    if not good_match(B_tnu, A) and not cl.get_accepted(B_tnu):
      point = get_graft_point(B_tnu, A)    # in A
      if point:
        graft_points[B_tnu] = point    # in A
    else:
      for child in cl.get_inferiors(B_tnu):
        process(child)
  for root in cl.get_roots(B):
    process(root)
  grafts = cl.invert_dict(graft_points)

def get_graftees(A_node):
  return (grafts.get(A_node) or [])

def get_graft_point(B_tnu, A = None):
  B_parent = cl.get_parent(B_tnu)
  if B_parent:
    B_parent_match = good_match(B_parent, A)
    if B_parent_match:
      return B_parent_match.cod
  return None

# ---------- One-sided best match

good_match_cache = {}

def good_match(node, other = None):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in good_match_cache:
    return good_match_cache[node]

  assert other
  matches = good_matches(node, other)
  match = best_match(matches)
  if match:
    if cl.get_accepted(match.cod) and not cl.get_accepted(match.dom):
      print("** Match is supposed to be terminal:\n  %s" %
            (art.express(match)))
    # Happens way often in GBIF
    if False and not cl.is_accepted(match.cod):
      print("** Match has taxonomic status %s\n  %s" %
            (cl.get_value(match.cod, cl.taxonomic_status_field),
             art.express(match)))

  good_match_cache[node] = match
  return match

good_matches_cache = {}

def good_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in good_matches_cache: return good_matches_cache[node]

  exties = extensional_matches(node, other)
  exties = [match_to_accepted(m) for m in exties]
  nameys = matches_to_accepted(node, other)
  matches = collapse_matches(exties + nameys)

  good_matches_cache[node] = matches
  return matches

# ---------- EXTENSIONALITY

# Starting with one match, extend to a set of matches by adding
# extensionally equivalent ancestors

def extensional_matches(tnu, other):
  assert cl.get_checklist(tnu) != other

  match = extensional_match(tnu, other)    # Single extensional match
  if not match: return []
  matches = [match]

  if rel.is_variant(match.relation, rel.extension_eq):
    # Scan upwards looking for nodes that return to us...
    scan = match.cod    # tnu -> partner
    here = cl.get_checklist(scan)
    while True:
      scan = cl.get_superior(scan)    # in other
      if not scan: break
      rev = extensional_match(scan, here)
      # rev: other -> here
      if not rev: break
      if not rel.is_variant(rev.relation, rel.extension_eq): break
      if rev.cod != tnu: break
      matches.append(rel.reverse(rev))
  matches.reverse()    # hmm. for choose ?
  return matches

# Guaranteed invertible, except for monotypic node chains

def extensional_match(tnu, other):
  assert tnu > 0
  here = cl.get_checklist(tnu)
  match = particle_match(tnu, other)
  if match:
    return match
  else:
    partner = cross_mrca(tnu, other)     # TNU in other checklist
    if not partner:
      return None
    match = bridge(tnu, partner, how_related_extensionally(tnu, partner))
    return match_to_accepted(match)

def how_related_extensionally(tnu, partner):
  here = cl.get_checklist(tnu)
  i = particle_match(tnu, partner)
  if i:
    return i.relation
  else:
    back = cross_mrca(partner, here)
    if not back:                # shouldn't happen
      assert False
      return rel.extension_disjoint
    elif cl.are_disjoint(tnu, back):
      return rel.extension_disjoint
    elif tnu == back:
      return rel.eq
    elif cl.mrca(tnu, back) == tnu:
      # Monotypic: back < tnu
      return rel.monotypic_in
    else:
      for sub in get_children(partner):
        back = cross_mrca(sub, here)
        if back:
          assert cl.get_checklist(tnu) == cl.get_checklist(back)
          if cl.mrca(tnu, back) == tnu and cross_disjoint(tnu, partner):
            return rel.extension_conflict
      return rel.extension_lt

# ---------- Determine disjointness across checklists

def cross_disjoint(tnu, partner):
  assert tnu > 0
  assert partner > 0
  assert cl.get_checklist(tnu) != cl.get_checklist(partner)

  back = cross_mrca(partner, cl.get_checklist(tnu))
  if not back:
    return True
  assert back > 0
  assert cl.get_checklist(back) == cl.get_checklist(tnu)
  if cl.are_disjoint(tnu, back):
    return True
  for inf in get_children(partner):
    assert inf > 0
    if not cross_disjoint(tnu, inf):
      return False
  return True

# ---------- Cross-MRCAs

# Returns a tnu - we never care about the reason for the match here

def cross_mrca(tnu, other):
  assert tnu > 0
  match = particle_match(tnu, other)
  if match:
    return match.cod
  else:
    return cross_mrcas.get(tnu)

def analyze_cross_mrcas(A, B):
  cross_mrcas = {}
  def half_analyze_cross_mrcas(checklist, other):
    def subanalyze_cross_mrcas(tnu, other):
      ind = particle_match(tnu, other)
      if ind:
        assert ind.dom == tnu
        if rel.is_variant(ind.relation, rel.eq):
          m = ind.cod
        else:
          print("** Non-eq articulation from particle_match\n  %s" +
                art.express(ind))
          assert False
      else:
        m = None
        for child in get_children(tnu):
          m2 = subanalyze_cross_mrcas(child, other)
          if m2:
            m = cl.mrca(m, m2)
      if m:
        assert cl.get_checklist(tnu) != cl.get_checklist(m)
        cross_mrcas[tnu] = m
        return m
      else:
        return None
    for root in cl.get_roots(checklist):
      subanalyze_cross_mrcas(root, other)
  half_analyze_cross_mrcas(A, B)
  half_analyze_cross_mrcas(B, A)
  return cross_mrcas

cross_mrcas = {}

# ---------- Particles

# A particle is mutual articulation deriving from sameness of 'intrinsic'
# node properties: name, id, rank, parent, etc.

def particle_match(tnu, other):
  return particles.get(tnu)

def find_particles(here, other):
  particles = {}
  count = [0]
  def log(tnu, message):
    if count[0] < 10:
      print("# %s: %s" % (cl.get_unique(tnu), message))
      count[0] += 1
  def subanalyze(tnu, other):
    log(tnu, "subanalyze")
    if cl.get_accepted(tnu):
      print("** Child %s of %s has an accepted name" %
            (cl.get_unique(tnu), cl.get_unique(cl.get_parent(tnu))))
      return False
    found_match = False
    for inf in get_children(tnu):
      log(tnu, "child %s" % inf)
      if subanalyze(inf, other):
        found_match = True
    if found_match:    # Some descendant is an particle
      return True
    candidate = best_match(matches_to_accepted(tnu, other))
    if candidate:
      rematch = best_match(matches_to_accepted(candidate.cod, here))
      if rematch:
        if rematch.cod == tnu:
          if cl.get_accepted(candidate.cod):
            print("** Candidate is synonymlike: %s" % cl.get_unique(candidate.cod))
          particles[tnu] = candidate    # here -> other
          particles[candidate.cod] = art.reverse(candidate)  # other -> here
          return True
        else:
          # This situation probably reflects a split!
          log(tnu,
              "Round trip fail:\n  %s\n  %s\n" %
              (art.express(candidate),
               art.express(art.reverse(rematch))))
      else:
        log(tnu, "No rematch")
    return False
  log(0, "top")
  for root in cl.get_roots(here):
    log(root, "root")
    subanalyze(root, other)
  return particles

# ---------- Matches to accepted nodes

# Three components: synonym-or-self o intensional o synonym-of-or-self
# from_accepted_articulations o intensional_matches o [to_accepted_articulation]

# Matches where the target is accepted.

def matches_to_accepted(tnu, other):
  return weak_matches_to_accepted(tnu, other)

def weak_matches_to_accepted(tnu, other):
  matches = strong_matches_to_accepted(tnu, other)
  syns = synonyms_locally(tnu)
  assert not syns or syns[0].dom == tnu
  matches = matches + [art.compose(syn, match)
                       for syn in syns
                       for match in strong_matches_to_accepted(syn.cod, other)
                       if art.composable(syn, match)]
  assert is_matches(matches)
  return collapse_matches(matches)

def strong_matches_to_accepted(tnu, other):
  # (a) Go intensional to accepted
  # (b) Go intensional from synonym to synonym, then up to accepted
  # (c) Go from synonym to accepted, then intensional to accepted

  matches = intensional_matches(tnu, other)
  assert is_matches(matches)
  matches = [match_to_accepted(match) for match in matches]
  assert is_matches(matches)
  if cl.get_accepted(tnu):
    accept = has_accepted_locally(tnu)
    assert accept.dom == tnu
    more = [art.compose(accept, match) for match in intensional_matches(accept.cod, other)]
    matches = matches + more
  return collapse_matches(matches)

# Convert match a->b to a->b' where b' is accepted (maybe b=b')

def match_to_accepted(m):
  if cl.get_accepted(m.cod):
    return art.compose(m, has_accepted_locally(m.cod))
  else:
    return m


# Intensional matches by name (no synonym following)
# This ought to be cached I think?

def intensional_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  hits = cl.get_tnus_with_value(other,
                                cl.canonical_name_field,
                                cl.get_name(node))
  matches = [bridge(node, hit, rel.same_name) for hit in hits]
  assert is_matches(matches) and True
  id_hit = cl.get_tnu_with_id(other, cl.get_tnu_id(node))
  if id_hit and not id_hit in hits:
    matches = [bridge(node, id_hit, rel.same_id)] + matches
  assert is_matches(matches)
  rank_matches = [bridge(node, match.cod, rel.same_rank)
                  for match in matches
                  if cl.get_nominal_rank(match.cod) == cl.get_nominal_rank(node)]
  return sort_by_badness(collapse_matches(matches + rank_matches))    # Should be cached.

# ---------- Within-checklist articulations

# Synonym-or-self = to-accepted
# Matches that involve going from an accepted name to a synonym are weak

# An accepted name can have many synonyms

def synonyms_locally(tnu):
  if cl.get_accepted(tnu):
    return []
  else:
    hases = [has_accepted_locally(syn) for syn in cl.get_synonyms(tnu)]
    return collapse_matches([art.reverse(ar) for ar in hases if ar])

# A synonym has only one accepted name

def has_accepted_locally(maybe_syn):   # goes from synonym to accepted
  assert maybe_syn > 0
  if cl.get_accepted(maybe_syn):
    accepted = cl.get_accepted(maybe_syn)
    if accepted:
      if accepted != maybe_syn:
        status = rel.synonym_relation(maybe_syn_status(maybe_syn))
        return art.art(maybe_syn, accepted, status)
      else:
        print("** Shouldn't happen", cl.get_unique(maybe_syn))
        return art.identity(maybe_syn)
    else:
      print("** Synonym %s has no accepted name" % cl.get_unique(maybe_syn))
      return None
  else:
    return None

def maybe_syn_status(synonym):
  return cl.get_value(synonym, cl.nomenclatural_status_field) or \
         cl.get_value(synonym, cl.taxonomic_status_field) or \
         "synonym"

# ---------- Best match among a set with common domain

# Assumes already sorted by art.badness

def best_match(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  arts = best_matches(sort_by_badness(collapse_matches(arts)))
  b = arts[0]
  if len(arts) == 1: return b
  print("** Multiple least-bad matches. Need to find more tie-breakers.")
  print("   %s -> %s" %
        (cl.get_unique(b.dom),
         [cl.get_unique(a.cod) for a in arts]))
  print(', '.join(["%s" % ar.relation.name for ar in arts]))
  return None

# There can be multiple best matches

def best_matches(sorted_matches):
  if len(sorted_matches) == 0:
    return []
  else:
    badness = art.badness(sorted_matches[0])
    matches = []
    for match in sorted_matches:
      if art.badness(match) == badness:
        matches.append(match)
      else:
        break
    return matches

def sort_by_badness(arts):
  return sorted(arts, key=art.badness)

# ---------- Utility: collapsing a set of matches

# Reduce a set of articulations so that all the codomains are
# different

def collapse_matches(arts):
  if len(arts) <= 1: return arts
  arts = sorted(arts, key=art.conjoin_sort_key)
  previous = None
  matches = []
  for ar in arts:
    assert ar
    if previous and art.conjoinable(previous, ar):
      previous = art.conjoin(previous, ar)
    else:
      if previous:
        matches.append(previous)
      previous = ar
  if previous:
    matches.append(previous)
  assert len(matches) <= len(arts)
  return matches

def bridge(dom, cod, re):
  assert cl.get_checklist(dom) != cl.get_checklist(cod)
  return art.art(dom, cod, re)

# ---------- For debugging

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

def cache_it(fun, tnu, other, dict):
  not_cached = "not cached"
  probe = dict.get(tnu, not_cached)
  if probe != not_cached:
    return probe
  else:
    result = fun(tnu, other)
    dict[x] = result
    return result

# ----

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', help='A checklist')
  parser.add_argument('B', help='B checklist')
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  args = parser.parse_args()
  main(args.A, args.B, args.out)
