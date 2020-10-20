#!/bin/env python3

debug = False

import sys, csv
import argparse

import checklist as cl
import relation as rel
import articulation as art
import eulerx
import alignment
import diff
import merge
import dribble

# A is lower priority, B is higher

def main(c1, c1_tag, c2, c2_tag, out, format):
  global dribble_file
  A = cl.read_checklist(c1, c1_tag + ".", "low-checklist")
  B = cl.read_checklist(c2, c2_tag + ".", "high-checklist")
  print ("Node counts:", len(A.get_all_nodes()), len(B.get_all_nodes()))
  with open("dribble.txt", "w") as dribfile:
    dribble.dribble_file = dribfile
    # Map each B to a corresponding A
    (al, xmrcas) = alignment.align(B, A)
    # Where to xmrcas come from?
    write_report(A, B, al, xmrcas, format, out)
    dribble.dribble_file = sys.stdout

def write_report(A, B, al, xmrcas, format, outpath):
  if outpath == "-":
    really_write_report(A, B, al, xmrcas, format, sys.stdout)
  else:
    with open(outpath, "w") as outfile:
      print ("Preparing:", outpath)
      really_write_report(A, B, al, xmrcas, format, outfile)

def really_write_report(A, B, al, xmrcas, format, outfile):
  if format == "eulerx":
    eulerx.dump_alignment(al, outfile)
  else:
    (parents, roots) = merge.merge_checklists(A, B, al, xmrcas)
    print ("# Number of roots in merge: %s" % len(roots))
    print ("# Number of non-roots in merge: %s" % len(parents))
    report(A, B, al, roots, parents, outfile)

# Default (simplified) report format

def report(A, B, al, roots, parents, outfile):
  writer = csv.writer(outfile)
  writer.writerow(["indent", "operation", "dom", "dom id", "relation", "cod id", "cod", "unchanged", "changed_props", "reason", "rank"])
  children = cl.invert_dict(parents)
  all_props = set.intersection(set(A.properties), set(B.properties))
  changed = find_changed_subtrees(roots, children, all_props)
  def process(mnode, indent):

    ch = None
    status = changed.get(mnode)

    (x, y) = mnode
    childs = children.get(mnode, [])
    re = None
    dif = None
    reason = None
    if x and y:
      op = "SHARED"
      re = al.get(x).relation.name  # ~ ≈ or =
      comparison = diff.differences(x, y, all_props)
      if not diff.same(comparison):
        props = diff.unpack(comparison)
        dif = ("; ".join(map(lambda x:x.pet_name, props)))
      px = cl.get_parent(x)
      qx = al.get(px)
      if qx and qx.cod != cl.get_parent(y):
        op += (" (moved from %s)" % cl.get_node_id(px))
      reason = art.reason(al[x])

      if not status:
        if len(childs) > 0:
          ch = "subtree="
        else:
          ch = "tip"

    elif x:
      op = "A ONLY"
      ar = al.get(x)
      if ar:
        re = ar.relation.name
        reason = art.reason(ar)
        # Equivalence, usually, but sometimes not
        y = ar.cod

        if rel.is_variant(ar.relation, rel.eq):
          op += " (merge)"
        elif rel.is_variant(ar.relation, rel.conflict):
          op += " (conflict)"
        elif rel.is_variant(ar.relation, rel.lt):
          op += " (loss of resolution)"
    else:                       # y
      op = "B ONLY"
      ar = al.get(y)
      if ar:
        re = ar.relation.revname
        reason = art.reason(ar)
        x = ar.cod
        if rel.is_variant(ar.relation, rel.eq):
          op += " (split)"
        elif rel.is_variant(ar.relation, rel.conflict):
          op += " (reorganization)"
        elif rel.is_variant(ar.relation, rel.lt):
          op += " (increased resolution)"

    report_one_articulation(op, ch, dif, x, re, y, reason, writer, indent)
    jndent = indent + "__"
    if status:
      for child in childs:
        process(child, jndent)
  for root in roots:
    process(root, "")

def report_one_articulation(op, ch, dif, x, re, y, reason, writer, indent):
  ux = None
  ix = None
  uy = None
  iy = None
  if x:
    ux = cl.get_unique(x)
    ix = cl.get_node_id(x)
    rank = cl.get_nominal_rank(x)
  if y:
    uy = cl.get_unique(y)
    iy = cl.get_node_id(y)
    rank = cl.get_nominal_rank(y)
  writer.writerow([indent, op,
                   ux, ix, re,
                   iy, uy, ch, dif, reason, rank])

# --------------------
# utilities

# Returns table with True for merged nodes all of whose descendants are
# unchanged

def find_changed_subtrees(roots, children, all_props):
  status = {}
  def process(node):
    node_changed = False
    (x, y) = node
    if not x or not y:
      node_changed = True
    else:
      comparison = diff.differences(x, y, all_props)
      if not diff.same(comparison):
        node_changed = True
    descendant_changed = False
    for child in children.get(node, []):
      if process(child):
        descendant_changed = True
    status[node] = descendant_changed       # Cache it
    return descendant_changed or node_changed
  for root in roots:
    status[root] = process(root)
  count = 0
  for key in status:
    if status[key]: count += 1
  print("# Changed status: %s, changed: %s" % (len(status), count))
  return status

# --------------------

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('low', help='lower priority checklist')
  parser.add_argument('high', help='higher priority checklist')
  parser.add_argument('--low-tag', default="A")
  parser.add_argument('--high-tag', default="B")
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  parser.add_argument('--format', help='report format', default='ad-hoc')
  args = parser.parse_args()
  main(args.low, args.low_tag, args.high, args.high_tag,
       args.out, args.format)

