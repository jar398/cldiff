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

# A is lower priority, B is higher

def main(c1, c1_tag, c2, c2_tag, out, format):
  A = cl.read_checklist(c1, c1_tag + ".", "left-checklist")
  B = cl.read_checklist(c2, c2_tag + ".", "right-checklist")
  print ("Taxname counts:", len(cl.get_all_taxnames(A)), len(cl.get_all_taxnames(B)))
  write_report(A, B, format, out)

def write_report(A, B, format, outpath):
  if outpath == "-":
    really_write_report(A, B, format, sys.stdout)
  else:
    with open(outpath, "w") as outfile:
      print ("Preparing:", outpath)
      really_write_report(A, B, format, outfile)

def really_write_report(A, B, format, outfile):
  global grafts
  # Map each B to a corresponding A
  (al, xmrcas) = alignment.align(B, A)
  # Where to xmrcas come from?

  if format == "eulerx":
    eulerx.dump_alignment(al, outfile)
  elif format == "classic":
    grafts = merge.analyze_unmatched(B, A, al)
    print ("# Number of grafts: %s\n" % len(grafts))
    report_to_io(A, al, outfile)
  else:
    (parents, roots) = merge.merge_checklists(A, B, al, xmrcas)
    print ("# Number of roots in merge: %s" % len(roots))
    print ("# Number of non-roots in merge: %s" % len(parents))
    report(A, B, al, roots, parents, outfile)

# Default (simplified) report format

def report(A, B, al, roots, parents, outfile):
  inv = invert_alignment(al)
  writer = csv.writer(outfile)
  writer.writerow(["indent", "operation", "unchanged", "dom", "dom id", "relation", "cod id", "cod"])
  children = cl.invert_dict(parents)
  changed = find_changed_merged_subtrees(roots, children)
  def process(mnode, indent):
    (x, y) = mnode
    childs = children.get(mnode, [])
    re = None
    if x and y:
      op = "KEEP"
      # TBD: Add info about changed properties and parent
    elif x:
      op = "DELETE A"
      ar = al.get(x)
      if ar:
        re = ar.relation.name
        # Equivalence, usually, but sometimes not
        y = ar.cod

        if rel.is_variant(ar.relation, rel.eq):
          op += " (merge)"
        elif rel.is_variant(ar.relation, rel.conflict):
          op += " (conflict)"
        elif rel.is_variant(ar.relation, rel.lt):
          op += " (loss of resolution)"
    else:                       # y
      op = "ADD B"
      ar = al.get(y)
      if ar:
        re = ar.relation.revname
        x = ar.cod
        if rel.is_variant(ar.relation, rel.eq):
          op += " (split)"
        elif rel.is_variant(ar.relation, rel.conflict):
          op += " (reorganization)"
        elif rel.is_variant(ar.relation, rel.lt):
          op += " (increased resolution)"

    status = changed.get(mnode)
    ch = None
    if not status:
      if len(childs) == 0:
        ch = "tip"
      else:
        ch = "subtree="

    report_one_articulation(op, ch, x, re, y, writer, indent)
    jndent = indent + "__"
    if ch:
      for child in childs:
        process(child, jndent)
  for root in roots:
    process(root, "")

def report_one_articulation(op, ch, x, re, y, writer, indent):
  if x and y:
    writer.writerow([indent, op, ch,
                     cl.get_unique(x),
                     cl.get_taxname_id(x),
                     re,
                     cl.get_taxname_id(y),
                     cl.get_unique(y)])
  elif x:
      writer.writerow([indent, op, ch,
                       cl.get_unique(x),
                       cl.get_taxname_id(x),
                       re, None, None])
  else:
    writer.writerow([indent, op, ch,
                     None, None, re,
                     cl.get_taxname_id(y),
                     cl.get_unique(y)])

# --------------------
# utilities

def invert_alignment(alignment):
  inv = {}
  for ar in alignment.values():
    rev = art.reverse(ar)
    if rev.dom in inv:
      inv[rev.dom].append(rev)
    else:
      inv[rev.dom] = [rev]
  return inv

# Returns table with True for merged nodes all of whose descendants are
# unchanged

def find_changed_merged_subtrees(roots, children):
  status = {}
  def process(node):
    node_changed = False
    (x, y) = node
    if not x or not y:
      node_changed = True
    descendant_changed = False
    for child in children.get(node, []):
      if process(child):
        descendant_changed = True
    status[node] = descendant_changed       # Cache it
    return descendant_changed or node_changed
  for root in roots:
    process(root)
  return status


# Returns table with True for nodes all of whose descendants are
# unchanged

def find_changed_subtrees(A, inv):
  status = {}
  def process(node):
    node_changed = True
    matches = inv.get(node)
    if matches and len(matches) == 1:
      match = matches[0]
      if rel.is_variant(match.relation, rel.eq):
        comparison = diff.differences(match.dom, match.cod)
        if diff.same(comparison):
          node_changed = False
    descendant_changed = False
    for child in alignment.get_children(node):
      if process(child):
        descendant_changed = True
    status[node] = descendant_changed       # Cache it
    return descendant_changed or node_changed
  for root in cl.get_roots(A):
    process(root)
  return status

# --------------------

# Fancier report format, buggy

def report_to_io(A, al, outfile):
  global descendant_changed

  # For each node, find all nodes that align to it
  inverse_alignment = invert_alignment(al)
  descendant_changed = find_changed_subtrees(A, invert_alignment)

  writer = csv.writer(outfile)
  write_header(writer)
  sink = make_sink(writer, None)
  for root in cl.get_roots(A):
    subreport(root, al, inverse_alignment, sink, "")
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
    proclaim(parent, indent, "...", None, None, None,
             "%s unchanged children" % len(rows))

def write_header(writer):
  writer.writerow(["tag", "left name", "left id", "rcc5", "right id", "right name", "justification"])

no_change_tags = ["NO CHANGE", "CHANGED ID", "..."]

# Report generation.  Recursion over the A checklist.
# TBD:
#  1. Separate nodes that match this one into = vs. other
#  2. Show all = matches preceding children and grafts

def subreport(node, al, inv, sink, indent):
  friends = report_on_matches(node, al, inv, sink, indent)
  subreport_friends(node, friends, al, inv, sink, indent)

# report_on_matches:
# Return value is a list of articulations that will get reported on
# separate lines, as "friends".
# Node is in A, matches are in B...
# Either (1) assume B is already populated, or
#        (2) explicitly add things B has that A doesn't

def report_on_matches(node, al, inv, sink, indent):

  # Which B-nodes have this A-node as their best match?

  matches = art.sort_matches(inv.get(node) or [])

  # Of these matches, which one, if any, is mutual?

  equivalents = [m for m in matches if rel.is_variant(m.relation, rel.eq)]
  relateds = [m for m in matches if not rel.is_variant(m.relation, rel.eq)]
  if len(equivalents) > 0:
    for option in equivalents:
      report_on_match(option, al, sink, indent)
    return relateds             # report indented ...?
  else:
    # relateds are the non-= matches
    if len(relateds) == 0:
      proclaim(sink, indent, "DELETE",
               node, "", None, None)
    else:
      # Relateds are non-=
      # Error or something if there are multiple equal B-nodes?
      b = relateds[0]
      friend = b.cod if b != None else None
      proclaim(sink, indent, "LUMP",
               node, "", friend, None)
    return []

# Report on A-children and other non-=-related B-nodes (mainly grafts)

def subreport_friends(node, friends, al, inv, sink, indent):
  sink = subsink(sink)
  def for_seq(node):
    b = al[node]
    return b.cod if b else node
  def sort_key(triple):
    (B_node, which, arg) = triple
    return cl.get_sequence_number(B_node)

  if descendant_changed.get(node):
    ch = [(for_seq(child), 0, child) for child in alignment.get_children(node)]
  else:
    ch = []

  agenda = ch + \
    [(friend.cod, 1, friend) for friend in friends] + \
    [(B_node, 2, B_node) for B_node in get_graftees(node)]
  indent = indent + "__"
  for (B_node, which, arg) in \
     sorted(agenda, key=sort_key):
    if which == 0:              # child
      subreport(arg, al, inv, sink, indent)
    elif which == 1:            # friend
      report_on_match(arg, al, sink, indent)
    elif which == 2:            # graft / addition
      proclaim(sink,
               indent,
               "ADD",    # tbd: add "-SUBTREE" when not tip
               "", "", arg, "")
  drain(sink)

def get_graftees(A_node):
  return (grafts.get(A_node) or [])     # grafts is global

# optionsp means ... it's not mutual

def report_on_match(match, al, sink, indent):
  tag = tag_for_match(match, al)
  proclaim(sink, indent, tag,
           match.dom,
           match.relation.name,
           match.cod,
           art.get_comment(match))

def proclaim(sink, indent, tag, dom, re, cod, remark):
  if not tag in no_change_tags: sink[0] = True
  proclaim_row((indent, tag, dom, re, cod, remark), sink)

def proclaim_row(row, sink):
  (changed, rows, writer, parent) = sink
  if parent:
    rows.append(row)
  else:
    (indent, tag, dom, re, cod, remark) = row
    (dom_name, dom_id) = get_taxname_info(dom)
    (cod_name, cod_id) = get_taxname_info(cod)
    writer.writerow([indent + tag, dom_name, dom_id, re, cod_id, cod_name, remark])

def get_taxname_info(tnu):
  if tnu:
    prefix = cl.get_checklist(tnu).prefix
    status = cl.get_taxonomic_status(tnu)
    modifiers = "?" if status == "synonym" else ""
    return (cl.get_name(tnu),
            prefix + cl.get_taxname_id(tnu) + modifiers)
  else:
    return (None, None)

def tag_for_match(match, al):
  mutualp = alignment.is_mutual(match, al)
  tag = "?"
  if rel.is_variant(match.relation, rel.eq):
    if mutualp:
      if descendant_changed.get(match.dom):
        tag = "KEEP"
      else:
        if len(alignment.get_children(match.dom)) == 0:
          tag = "KEEP-TIP"
        else:
          tag = "KEEP-SUBTREE"
      comparison = diff.differences(match.dom, match.cod)
      if not diff.same(comparison):
        props = diff.unpack(comparison)
        tag = "%s (change %s)" % (tag,
                                  "; ".join(map(lambda x:x.pet_name, props)))
    else:
      # Make a list of changes: parent, name, etc.
      if cl.is_accepted(match.cod):
        tag = "SPLIT"
      else:
        tag = "ADD-SYNONYM"
  elif rel.is_variant(match.relation, rel.lt):
    tag = "ELIDE"
  elif rel.is_variant(match.relation, rel.gt):
    tag = "INSERT"
  elif rel.is_variant(match.relation, rel.conflict):
    tag = "BREAK" if mutualp else "REFORM"
  elif rel.is_variant(match.relation, rel.disjoint):
    tag = "MUTEX"    # shouldn't happen ...
  return tag


# X aligns to Y.  Does parent X align to parent Y?

# def parent_changed(match):
#   parent = cl.get_parent(match.dom)
#   coparent = cl.get_parent(match.cod)
#   if not parent and not coparent:
#     return False
#   if not parent or not coparent:
#     return True
#   other = cl.get_checklist(coparent)
#   match = good_candidate(parent, ?other?)
#   if not match: return True
#   # if relation is not = : return True
#   return match.cod != coparent

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('left', help='A checklist')    # Lower priority
  parser.add_argument('right', help='B checklist')   # Higher priority
  parser.add_argument('--left-tag', default="A")
  parser.add_argument('--right-tag', default="B")
  parser.add_argument('--share_ids', default=False)
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  parser.add_argument('--format', help='report format', default='ad-hoc')
  args = parser.parse_args()
  alignment.shared_idspace = args.share_ids
  main(args.left, args.left_tag, args.right, args.right_tag, args.out, args.format)

