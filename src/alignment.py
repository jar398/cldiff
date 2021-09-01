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
  proposal = propose_alignment(proposal, best, xmrcas)
  return (proposal, xmrcas)

def propose_alignment(proposal, best, xmrcas):
  for x in xmrcas:
    alignment_step(x, best, xmrcas, proposal)
  return proposal

# This side-affects the proposal.

# Look at every set of mutually ~-related nodes in xmrcas
# Store potential = and < relationships into the proposal

def alignment_step(x, best, xmrcas, proposal):
  ar = xmrcas.get(x)          # x ~ y0
  if not ar: return
  y0 = ar.cod
  x0 = xmrcas.get(y0).cod

  # Three cases:
  #  1. Bottom of a ladder, x = x0
  #  2. On ladder, but not at bottom: x > x0
  #  3. Not on a ladder: x < x0

  if x != x0:
    if xmrcas.get(x0).cod != y0:     # x0 < x < y0
      art.proclaim(proposal, art.extensional(x, y0, rel.lt, "inferred"))
    return

  if x0 > y0: return

  # x and y start at x0 and y0 and ratchet toward root of tree
  def luup(x, y):
    if x == cl.forest_tnu or y == cl.forest_tnu: return
    assert cl.get_checklist(x) != cl.get_checklist(y)

    if xmrcas.get(y).cod != x0:
      # if xmrcas.get(x).cod != y0: return      # Processed all of both sides of ladder
      # x is in ladder but y isn't, so x < y
      # art.proclaim(proposal, art.extensional(x, y, rel.lt, "refines", "refined by"))
      # luup(cl.get_parent(x), y)
      return

    elif xmrcas.get(x).cod != y0:
      # y is in ladder but x isn't, so y < x
      # art.proclaim(proposal, art.extensional(y, x, rel.lt, "refines", "refined by"))
      # luup(x, cl.get_parent(y))
      return

    else:
      xpro = proposal.get(x)
      ypro = proposal.get(y)
      if False:
        if xpro: assert xpro.cod == y
        if ypro and ypro.cod != x:
          print("## x0 %s, x %s, y0 %s, y %s,\n##  ypro %s" %
                (cl.get_unique(x0), cl.get_unique(x),
                 cl.get_unique(y0), cl.get_unique(y),
                 art.express(ypro)))
        if ypro: assert ypro.cod == x

      # Careful, we might have x < y or y < x in the proposal already!
      if xpro and xpro.relation == rel.lt:
        luup(cl.get_parent(x), y)
      elif ypro and ypro.relation == rel.lt:
        luup(x, cl.get_parent(y))
      else:
        xar = get_mutual(best, x)

        if xar and xar.cod == y:
          art.proclaim_eq(proposal, art.change_relation(xar, rel.eq, "extensional"))
          luup(cl.get_parent(x), cl.get_parent(y))

        elif xar and find_in_ladder(xar.cod, y):
          # x = y1 > y, so y < x
          # art.proclaim(proposal, art.extensional(y, x, rel.lt, "refines", "refined by"))
          luup(x, cl.get_parent(y))

        else:
          yar = get_mutual(best, y)

          if yar and find_in_ladder(yar.cod, x):
            # x < x1 = y, so x < y
            # art.proclaim(proposal, art.extensional(x, y, rel.lt, "refines", "refined by"))
            luup(cl.get_parent(x), y)

          else:
            # Neither x nor y matches intensionally.  Set up a blind date.
            art.proclaim_eq(proposal, art.extensional(x, y, rel.eq, "similar="))
            luup(cl.get_parent(x), cl.get_parent(y))

  # See if b is in the ladder (matching nodes in lineage)
  def find_in_ladder(b, y):
    if b == y:
      return True
    elif y == cl.forest_tnu:
      return False
    else:
      ar = xmrcas.get(y)
      if ar and ar.cod == x0:
        return find_in_ladder(b, cl.get_parent(y))
      else:
        return False

  luup(x0, y0)

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
      if y != None:             # x <= y
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
