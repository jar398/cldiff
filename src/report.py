debug = False

import sys, csv
import argparse

import checklist as cl
import relation as rel
import articulation as art
import eulerx
import diff

def main(c1, c1_tag, c2, c2_tag, out, format):
  A = cl.read_checklist(c1, c1_tag + ".", "checklist 1")
  B = cl.read_checklist(c2, c2_tag + ".", "checklist 2")
  print ("TNU counts:", len(cl.get_all_tnus(A)), len(cl.get_all_tnus(B)))
  diff.start(A, B)
  write_report(A, B, format, out)

def write_report(A, B, format, outpath):
  if outpath == "-":
    really_write_report(A, B, format, sys.stdout)
  else:
    with open(outpath, "w") as outfile:
      print ("Writing:", outpath)
      really_write_report(A, B, format, outfile)

def really_write_report(A, B, format, outfile):
  if format == "eulerx":
    eulerx.dump_alignment(diff.finish_alignment(B, A), outfile)
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
    (dom_name, dom_id) = get_tnu_info(dom)
    (cod_name, cod_id) = get_tnu_info(cod)
    writer.writerow([indent + tag, dom_name, dom_id, re, cod_id, cod_name, remark])

def write_header(writer):
  writer.writerow(["tag", "left name", "left id", "rcc5", "right id", "right name", "justification"])

def get_tnu_info(tnu):
  if tnu:
    prefix = cl.get_checklist(tnu).prefix
    status = cl.get_taxonomic_status(tnu)
    modifiers = "?" if status == "synonym" else ""
    return (cl.get_name(tnu),
            prefix + cl.get_tnu_id(tnu) + modifiers)
  else:
    return (None, None)

# Report generation

def subreport(node, B, sink, indent):
  A = cl.get_checklist(node)
  assert A != B
  multiple = report_on_matches(node, B, sink, indent)
  sink = subsink(sink)
  def for_seq(node):
    b = diff.good_candidate(node, B)
    return b.cod if b else node
  def sort_key(triple):
    (B_node, which, arg) = triple
    return cl.get_sequence_number(B_node)
  agenda = \
    [(for_seq(child), 0, child) for child in diff.get_children(node)] + \
    [(option.cod, 1, option) for option in multiple] + \
    [(B_node, 2, B_node) for B_node in diff.get_graftees(node)]
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
                       arg,
                       "")
  drain(sink)

def report_on_matches(node, B, sink, indent):
  matches = diff.good_candidates(node, B)      # cod is accepted
  if len(matches) == 0:
    proclaim(sink, indent, "REMOVE",
                     node,
                     "",
                     "",
                     "%s right nodes match this left node" % len(matches))
    return []
  elif len(matches) == 1:
    report_on_match(matches[0], False, sink, indent)
    return []
  else:
    proclaim(sink, indent, "MULTIPLE",
                     node,
                     "?",
                     "",
                     "%s right nodes match this left node" % len(matches))
    return matches

def report_on_match(match, splitp, sink, indent):
  tag = tag_for_match(match, splitp)
  proclaim(sink, indent, tag,
                   None if splitp else match.dom,
                   rel.rcc5_name(match.relation),
                   match.cod,
                   art.get_comment(match))

def tag_for_match(match, splitp):
  tag = "?"
  if rel.is_variant(match.relation, rel.eq):
    if splitp:
      if cl.get_accepted(match.cod):
        tag = "ADD SYNONYM"
      else:
        tag = "OPTION"
    elif diff.parent_changed(match):
      tag = "MOVE"
    else:
      changes = []
      
      if match.differences != 0:
        tag = "CHANGED %o" % match.differences
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


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('left', help='A checklist')
  parser.add_argument('right', help='B checklist')
  parser.add_argument('--left-tag', default="A")
  parser.add_argument('--right-tag', default="B")
  parser.add_argument('--idspace', default=False)
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  parser.add_argument('--format', help='report format', default='ad-hoc')
  args = parser.parse_args()
  diff.shared_idspace = args.idspace
  main(args.left, args.left_tag, args.right, args.right_tag, args.out, args.format)

