import sys

import checklist as cl
import relation as rel
import articulation as art
import intension
import dribble
from intension import choose_best_match

# For each B-record, we choose an articulation with the closest
# A-record that it matches (preferably but not necessarily an '='
# articulation).

def align(B, A):

  # Precompute all best matches
  best = intension.best_intensional_match_map(B, A)

  # Extensional analysis yields <= relationships between hierarchies
  # (written as the 'matches' relation ~)
  xmrcas = infer_partners(best, A, B)
  dribble.log("# Number of cross-mrcas: %s" % len(xmrcas))

  # Turn tipward best matches into = or < articulations as appropriate
  proposal = intension.intensional_proposal(best, A, B)

  # Add extensional matches to a draft that already has intensional matches
  the_alignment = propose_alignment(proposal, best, xmrcas)
  return (the_alignment, xmrcas)

# Side-affects proposal

def propose_alignment(proposal, best, xmrcas):
  for x in xmrcas:
    if not proposal.get(x):
      alignment_step(x, best, xmrcas, proposal)
  return proposal

# Look at every set of mutually ~-related nodes in xmrcas
# Store potential = and < relationships into the proposal

def alignment_step(x, best, xmrcas, proposal):
  def luup(x, y):
    if x == cl.forest_tnu or y == cl.forest_tnu: return
    assert cl.get_checklist(x) != cl.get_checklist(y)
    if in_chain(x, y0) and in_chain(y, x0):
      bar = best.get(x)
      if bar and find_in_chain(bar.cod, y, x0):
        if bar.cod == y:
          art.proclaim(proposal, art.change_relation(bar, rel.eq, "extensional"))
          luup(cl.get_parent(x), cl.get_parent(y))
        else:
          art.proclaim(proposal, art.extensional(y, x, rel.lt, "refines", "refined by"))
          luup(x, cl.get_parent(y))
      else:
        bar = best.get(y)
        if bar and find_in_chain(bar.cod, x, y0):
          if bar.cod == y:
            art.proclaim(proposal, bar)
            luup(cl.get_parent(x), cl.get_parent(y))
          else:
            art.proclaim(proposal, art.extensional(x, y, rel.lt, "refines*", "refined by*"))
            luup(cl.get_parent(x), y)
        else:
          # neither x nor y matches by name
          art.proclaim(proposal, art.extensional(x, y, rel.eq, "similar="))
          luup(cl.get_parent(x), cl.get_parent(y))

  # See is b is in the chain (matching nodes in lineage)
  def find_in_chain(b, y, x0):
    if y == cl.forest_tnu:
      return False
    elif in_chain(y, x0):
      if b == y:
        return True
      else:
        return find_in_chain(b, cl.get_parent(y), x0)
    else:
      return False

  def in_chain(y, x0):
    ar = xmrcas.get(y)
    return ar and ar.cod == x0

  ar = xmrcas.get(x)          # x ~ y0
  if ar:
    y0 = ar.cod
    x0 = xmrcas.get(y0).cod
    if xmrcas.get(x0).cod != y0:     # x < y0
      art.proclaim(proposal, art.extensional(x, y0, rel.lt, "inferred"))
    else:                       # x could be <, =, > y0.  figure it out
      if x0 < y0:
        luup(x0, y0)
      else:
        pass    # pick it up later (or earlier)

# ---------- Cross-MRCAs (partners) x <= y

def infer_partners(best, A, B):
  xmrcas = {}
  def half_infer_partners(checklist, other):
    def subinfer_partners(x, other):
      y = None
      for child in cl.get_children(x):
        child_ar = subinfer_partners(child, other) # an articulation
        if child_ar != None:
          child_y = child_ar.cod
          if y == None:
            y = child_y
          else:
            y = cl.mrca(y, child_y)
      if y != None:
        ar = art.extensional(x, y, rel.matches, "cross-mrca")
      else:
        ar = get_mutual(best, x)
      if ar:
        assert cl.get_checklist(ar.cod) != cl.get_checklist(x)
        if dribble.watch(x):
          dribble.log("# Cross-mrca: %s" % (art.express(ar)))
        xmrcas[x] = ar
      return ar             # in B
    for root in cl.get_roots(checklist):
      subinfer_partners(root, other)
  half_infer_partners(A, B)
  half_infer_partners(B, A)
  return xmrcas

def get_mutual(best, x):
  bar = best.get(x)
  if bar:
    back = best.get(bar.cod)
    if back and back.cod == x:
      return bar
  return None
