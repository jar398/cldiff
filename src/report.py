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
  dribpath = out + ".log"
  with open(dribpath, "w") as dribfile:
    dribble.dribble_file = dribfile
    dribble.log ("\nLogging to %s" % (dribpath,))
    A = cl.read_checklist(c1, c1_tag + ".", "low-checklist")
    B = cl.read_checklist(c2, c2_tag + ".", "high-checklist")
    dribble.log ("Node counts: %s %s" % (len(A.get_all_nodes()), len(B.get_all_nodes())))
    # Map each B to a corresponding A
    dribble.log ("Aligning ...")
    (al, xmrcas) = alignment.align(B, A)
    dribble.log("  ... finished aligning; %s articulations\n" %
                len(al))
    # Where do xmrcas come from?
    write_report(A, B, al, xmrcas, format, out)
    dribble.dribble_file = None

def write_report(A, B, al, xmrcas, format, outpath):
  if outpath == "-":
    really_write_report(A, B, al, xmrcas, format, sys.stdout)
  else:
    with open(outpath, "w") as outfile:
      really_write_report(A, B, al, xmrcas, format, outfile)

def really_write_report(A, B, al, xmrcas, format, outfile):
  if format == "eulerx":
    eulerx.dump_alignment(al, outfile)
  else:
    (parents, roots) = merge.merge_checklists(A, B, al)
    dribble.log ("Merged.  %s roots in merge, %s nodes with parents" %
                 (len(roots), len(parents)))
    report(A, B, al, roots, parents, outfile)
    report_on_collisions(A, B, al)

def assign_ids(parents, roots, children):
  id_table = {}
  def process(node):
    id_table[node] = len(id_table) + 1
    for child in children.get(node, []):
      process(child)
  for root in roots:
    process(root)
  return id_table

canonical_name = cl.field("canonicalName")

def report_on_collisions(A, B, al):
  index = cl.index_by_value(A, canonical_name)
  for name in index:
    A_nodes = index[name]
    if len(A_nodes) == 1:
      B_nodes = cl.get_nodes_with_value(B, canonical_name, name)
      if B_nodes and len(B_nodes) == 1:
        A_node = A_nodes[0]
        B_node = B_nodes[0]
        if cl.is_accepted(A_node) and cl.is_accepted(B_node):
          ar1 = al.get(A_node)
          ar2 = al.get(B_node)
          ar1_bad = (ar1 and
                     rel.is_variant(ar1.relation, rel.eq) and
                     ar1.cod != B_node)
          ar2_bad = (ar2 and
                     rel.is_variant(ar2.relation, rel.eq) and
                     ar2.cod != A_node)
          if ar1_bad or ar2_bad:
            dribble.log("# %s names different taxa in the two checklists" % name)
            if ar1_bad: dribble.log("  %s [%s]" % (art.express(ar1), art.reason(ar1)))
            if ar2_bad: dribble.log("  %s [%s]" % (art.express(ar2), art.reason(ar2)))

# Default (simplified) report format

def report(A, B, al, roots, parents, outfile):
  writer = csv.writer(outfile)
  write_header(writer)
  children = cl.invert_dict(parents)
  all_props = set.intersection(set(A.properties), set(B.properties))
  any_descendant_differs = find_changed_subtrees(roots, children, all_props)
  id_table = assign_ids(parents, roots, children)

  def taxon_report(mnode, indent):
    nodiff = None
    different = mnode in any_descendant_differs

    id = id_table[mnode]
    (x, y) = mnode
    re = None
    dif = None
    z = None
    note = None
    if x and y:
      op = "SHARED"
      ar = al.get(x)
      comparison = diff.differences(x, y, all_props)
      if not diff.same(comparison):
        props = diff.unpack(comparison)
        dif = ("; ".join(map(lambda x:x.pet_name, props)))
      px = cl.get_parent(x)
      qx = al.get(px)
      if qx and qx.cod != cl.get_parent(y):
        #note = ("moved from %s" % cl.get_node_id(px))
        note = "moved"

      if not different:
        childs = children.get(mnode, [])
        if len(childs) > 0:
          nodiff = "subtree="
        else:
          nodiff = "shared tip"

    elif x:
      op = "A ONLY"
      ar = al.get(x)
      if ar:
        z = ar.cod
        # Equivalence, usually, but sometimes not

        if rel.is_variant(ar.relation, rel.eq):
          note = "lump"
        elif rel.is_variant(ar.relation, rel.conflict):
          note = "conflict"
        elif rel.is_variant(ar.relation, rel.lt):
          note = "loss of resolution"
    else:                       # y
      op = "B ONLY"
      ar = al.get(y)
      if ar:
        z = ar.cod
        if rel.is_variant(ar.relation, rel.eq):
          note = "split"
        elif rel.is_variant(ar.relation, rel.conflict):
          note = "reorganization"
        elif rel.is_variant(ar.relation, rel.lt):
          note = "increased resolution"

    report_one_articulation(id, op, nodiff, dif, x, y, z, ar, note, writer, indent)
    return different

  def process(mnode, indent):
    different = taxon_report(mnode, indent)
    jndent = indent + "—"    # em dash
    if different:
      for child in children.get(mnode, []):
        process(child, jndent)
  for root in roots:
    process(root, "")

def report_one_articulation(id, op, nodiff, dif, x, y, z, ar, note, writer, indent):
  (ix, ux, rankx) = node_data(x)
  (iy, uy, ranky) = node_data(y)
  (iz, uz, _) = node_data(z)
  rank = rankx or ranky
  relation = ar.relation.name if ar else None
  reason = art.reason(ar) if ar else None
  writer.writerow([indent, id, #op,
                   rank, 
                   ix, ux, 
                   iy, uy, relation, iz, uz,
                   note, reason, dif, nodiff])

def write_header(writer):
  writer.writerow(["indent", "taxonID", #"operation",
                   "rank", 
                   "A id", "A name",
                   "B id", "B name", "relation", "other id", "other name",
                   "note", "reason", "changed_props", "unchanged"])

def node_data(node):
  if node:
    return (cl.get_node_id(node), cl.get_unique(node), cl.get_nominal_rank(node))
  else:
    return (None, None, None)

# --------------------
# utilities

# Returns table with True for merged nodes all of whose descendants are
# unchanged

def find_changed_subtrees(roots, children, all_props):
  any_descendant_differs = {}
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
    if descendant_changed:
      any_descendant_differs[node] = True
    return descendant_changed or node_changed
  for root in roots:
    c = process(root)
    if c: any_descendant_differs[root] = c
  dribble.log("# %s nodes in merge have some change in their descendants" %
              (len(any_descendant_differs)))
  return any_descendant_differs

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

