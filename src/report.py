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

# A is lower priority

def start(A, B):
  global inverse_good_candidates
  global grafts
  global the_alignment
  global change_status
  # Map each B to a corresponding A
  (the_alignment, grafts) = alignment.align(B, A)
  # Find all As that align to any given B
  inverse_good_candidates = invert_dict_by_cod(the_alignment)
  change_status = find_unchanged_subtrees(A, inverse_good_candidates)

# Partners list for reporting (A->B articulations)

def good_candidate(tnu, other):  # for children sort order
  return alignment.choose_best_match(good_candidates(tnu, other))

def good_candidates(tnu, other):
  matches = inverse_good_candidates.get(tnu) or []
  if len(matches) == 0:
    if debug:
      print("# matching: B-node", cl.get_unique(tnu))
    investigate = alignment.articulate(tnu, other, ...)
    if investigate: matches = [investigate]
  return [match for match in matches if not cl.get_accepted(match.cod)]

def invert_dict_by_cod(alignment):
  inv = {}
  for (node, ar) in alignment.items():
    if ar:      # ar: B -> A
      ar = art.reverse(ar)
      if ar.dom in inv:
        inv[ar.dom].append(ar)
      else:
        inv[ar.dom] = [ar]
  return inv

def find_unchanged_subtrees(A, inverse_good_candidates):
  status = {}
  def process(node):
    changed = True
    matches = inverse_good_candidates.get(node)
    if matches and len(matches) == 1:
      match = matches[0]
      if rel.is_variant(match.relation, rel.eq):
        comparison = diff.differences(match.dom, match.cod)
        if diff.same(comparison):
          changed = False
    chch = False
    for child in alignment.get_children(node):
      if process(child):
        chch = True
    status[node] = chch       # Cache it
    return chch or changed
  for root in cl.get_roots(A):
    process(root)
  return status

# --------------------

def main(c1, c1_tag, c2, c2_tag, out, format):
  A = cl.read_checklist(c1, c1_tag + ".", "left-checklist")
  B = cl.read_checklist(c2, c2_tag + ".", "right-checklist")
  print ("Taxname counts:", len(cl.get_all_taxnames(A)), len(cl.get_all_taxnames(B)))
  write_report(A, B, format, out)

def write_report(A, B, format, outpath):
  start(A, B)    # sets the_alignment etc.
  if outpath == "-":
    really_write_report(A, B, format, sys.stdout)
  else:
    with open(outpath, "w") as outfile:
      print ("Writing:", outpath)
      really_write_report(A, B, format, outfile)

def really_write_report(A, B, format, outfile):
  if format == "eulerx":
    eulerx.dump_alignment(the_alignment, outfile)
  else:
    report_to_io(A, B, outfile)

def report_to_io(A, B, outfile):
  writer = csv.writer(outfile)
  write_header(writer)
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
    proclaim(parent, indent, "...", None, None, None,
             "%s unchanged children" % len(rows))

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
    (dom_name, dom_id) = get_taxname_info(dom)
    (cod_name, cod_id) = get_taxname_info(cod)
    writer.writerow([indent + tag, dom_name, dom_id, re, cod_id, cod_name, remark])

def write_header(writer):
  writer.writerow(["tag", "left name", "left id", "rcc5", "right id", "right name", "justification"])

def get_taxname_info(tnu):
  if tnu:
    prefix = cl.get_checklist(tnu).prefix
    status = cl.get_taxonomic_status(tnu)
    modifiers = "?" if status == "synonym" else ""
    return (cl.get_name(tnu),
            prefix + cl.get_taxname_id(tnu) + modifiers)
  else:
    return (None, None)

# Report generation
# TBD:
#  1. Separate nodes that match this one into = vs. other
#  2. Show all = matches preceding children and grafts

def subreport(node, B, sink, indent):
  A = cl.get_checklist(node)
  assert A != B
  friends = report_on_matches(node, B, sink, indent)

  sink = subsink(sink)
  def for_seq(node):
    b = good_candidate(node, B)
    return b.cod if b else node
  def sort_key(triple):
    (B_node, which, arg) = triple
    return cl.get_sequence_number(B_node)

  if change_status[node] == False:
    ch = []
  else:
    ch = [(for_seq(child), 0, child) for child in alignment.get_children(node)]

  agenda = ch + \
    [(friend.cod, 1, friend) for friend in friends] + \
    [(B_node, 2, B_node) for B_node in get_graftees(node)]
  indent = indent + "__"
  for (B_node, which, arg) in \
     sorted(agenda, key=sort_key):
    if which == 0:              # child
      subreport(arg, B, sink, indent)
    elif which == 1:            # friend
      report_on_match(arg, True, sink, indent)
    elif which == 2:            # graft / addition
      proclaim(sink,
               indent,
               "ADD",    # tbd: add "-SUBTREE" when not tip
               "", "", arg, "")
  drain(sink)

def report_on_matches(node, B, sink, indent):
  matches = good_candidates(node, B)      # cod is accepted
  options = [m for m in matches if rel.is_variant(m.relation, rel.eq)]
  friends = [m for m in matches if not rel.is_variant(m.relation, rel.eq)]
  if len(options) > 0:
    best = alignment.choose_best_match(options)
    if best:
      options = [best]
      friends = [ok for ok in options if ok != best] + friends
  if len(options) == 0:
    proclaim(sink, indent, "DELETE",
                   node,
                   "",
                   "",
                   "%s right nodes match this left node" % len(matches))
  elif len(options) == 1:
    report_on_match(options[0], False, sink, indent)
  else:
    proclaim(sink, indent, "MULTIPLE", node,
             "?", "",
             "%s right nodes match this left node" % len(matches))
    for option in options:
      report_on_match(option, True, sink, indent)
  return friends

def get_graftees(A_node):
  return (grafts.get(A_node) or [])     # grafts is global

# optionsp means ... it's part of a multiple-match set

def report_on_match(match, optionsp, sink, indent):
  tag = tag_for_match(match, optionsp)
  proclaim(sink, indent, tag,
           None if optionsp else match.dom,
           rel.rcc5_name(match.relation),
           match.cod,
           art.get_comment(match))

def tag_for_match(match, optionsp):
  tag = "?"
  if rel.is_variant(match.relation, rel.eq):
    if optionsp:
      if cl.get_accepted(match.cod):
        tag = "ADD-SYNONYM"
      else:
        tag = "OPTION"
    else:
      if change_status[match.dom] == False:
        if len(alignment.get_children(match.dom)) == 0:
          tag = "KEEP-TIP"
        else:
          tag = "KEEP-SUBTREE"
      else:
        tag = "KEEP"
      comparison = diff.differences(match.dom, match.cod)
      if not diff.same(comparison):
        props = diff.unpack(comparison)
        tag = "%s (change %s)" % (tag,
                                  "; ".join(map(lambda x:x.pet_name, props)))
      # Make a list of changes: parent, name, etc.
  elif rel.is_variant(match.relation, rel.lt):
    tag = "ELIDE"
  elif rel.is_variant(match.relation, rel.gt):
    tag = "INSERT"
  elif rel.is_variant(match.relation, rel.conflict):
    tag = "REFORM" if optionsp else "BREAK"
  elif rel.is_variant(match.relation, rel.disjoint):
    tag = "MUTEX"    # shouldn't happen ...
  return tag

# X aligns to Y.  Does parent X align to parent Y?

def parent_changed(match):
  parent = cl.get_parent(match.dom)
  coparent = cl.get_parent(match.cod)
  if not parent and not coparent:
    return False
  if not parent or not coparent:
    return True
  other = cl.get_checklist(coparent)
  match = good_candidate(parent, other)
  if not match: return True
  # if relation is not = : return True
  return match.cod != coparent

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

