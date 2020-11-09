# Difference file in the form of commands that transform A into B (e.g. for EOL)

import csv
import re

import checklist as cl
import changes
import dribble
import property
import relation as rel
import articulation as art

def write_diff_set(A, B, al, keyprop, outpath):
  with open(outpath, "w") as outfile:
    print("Preparing diff")
    writer = csv.writer(outfile)
    write_row("operation", "taxonID", "property", "oldvalue", "newvalue", writer)

    def process_B_node(y):
      yid = merged_id(y, al, keyprop)
      x = eq_partner(y, al)
      q = cl.get_parent(y)
      if x:
        xid = merged_id(x, al, keyprop)
        (drop, change, add) = changes.differences_in_record(x, y)
        if drop != 0:
          write_changes(drop, "drop", xid, x, y, keyprop, writer)
        if change:    # some property changed going from x to y:
          write_changes(change, "change", xid, x, y, keyprop, writer)
        if add:
          write_changes(add, "add", xid, x, y, keyprop, writer)
        # Deal with topology changes (insertions)
        if q != cl.forest_tnu:
          p = cl.get_parent(x)
          pid = merged_id(p, al, keyprop)
          qid = merged_id(q, al, keyprop)
          if pid != qid:
            write_row("move", xid, None, pid, qid, writer)
      else:    # there's a y, but no x
        # Either the page id is 'repurposed' or the node is new
        if get_node_with_key(A, keyprop, get_key(y, keyprop)):
          # See "repurpose" below
          pass
        else:
          write_row("new", None, None, None, yid, writer)
      for child in cl.get_children(y): process_B_node(child)
    for y in cl.get_roots(B): process_B_node(y)
    def process_A_node(x):
      childs = cl.get_children(x)
      if not eq_partner(x, al):
        xid = merged_id(x, al, keyprop)
        z = get_node_with_key(B, keyprop, get_key(x, keyprop))
        if z:
          # Would be useful to show the relationship between the two...
          write_row("repurpose_id", xid, None, None, qualified_id(z, keyprop), writer)
        else:
          # Show what should replace it
          ar = al.get(x)
          zid = None; relate = None
          if ar:
            relate = ar.relation.name
            z = ar.cod
            zid = merged_id(z, al, keyprop)
          # Also show art.express(al.get(x)) ??
          write_row("deprecate", xid, relate, None, zid, writer)
      for child in childs: process_A_node(child)
    for x in cl.get_roots(A): process_A_node(x)
  print("Wrote %s" % outpath)

# N.b. xid == get_page_id(x), possibly prefixed

def write_row(op, mid, prop, before, after, writer):
  # Ultimately we don't need the 'before' value.  Just for debugging.
  writer.writerow([mid, op, prop, before, after])

def write_changes(mask, op, mid, x, y, keyprop, writer):
  assert mask > 0
  props = changes.unpack1(mask)
  for prop in props:
    u = cl.get_value(x, prop)
    v = cl.get_value(y, prop)
    if prop == cl.taxon_rank and u == "subspecies" and v == "infraspecies":
      pass
    elif prop == cl.taxonomic_status and u == "accepted" and v == "valid":
      pass
    elif prop == keyprop and u:
      write_row("do not " + op, mid, prop.pet_name, u, v, writer)
    elif (prop == cl.scientific_name and
          u and v and scipat.search(u) and not scipat.search(v)):
      # Lots of these, just drop them silently
      # dribble.log("Losing scientific name: %s -> %s" % (u, v))
      pass
    else:
      write_row(op, mid, prop.pet_name, u, v, writer)

scipat = re.compile(" [12][0-9][0-9][0-9]\\b")

# Id is primary key of A node, if there is one, otherwise nonconflicting
# key in B, otherwise is made up on the fly

def merged_id(x, al, keyprop):
  if x and x != cl.forest_tnu:
    x_page = get_key(x, keyprop)
    if x_page:
      # See if there's a conflict in the opposite checklist
      other = get_other(x, al)
      if (get_node_with_key(other, keyprop, x_page) == 
          eq_partner(x, al)):
        return x_page
      else:
        return cl.get_checklist(x).prefix + x_page
    else:
      return cl.get_checklist(x).prefix[0:-1] + "#" + cl.get_taxon_id(x)

def qualified_id(x, keyprop):
  x_page = get_key(x, keyprop)
  if x_page:
    return cl.get_checklist(x).prefix + x_page
  else:
    return cl.get_checklist(x).prefix[0:-1] + "#" + cl.get_taxon_id(x)

def get_other(x, al):
  while x != cl.forest_tnu: 
    ar = al.get(x)
    if ar:
      return cl.get_checklist(ar.cod)
    x = cl.get_parent(x)
  return None   # hope it doesn't happen

def get_key(node, keyprop):
  if keyprop == None: return None
  return cl.get_value(node, keyprop)

def get_node_with_key(checklist, keyprop, id):
  if keyprop == None: return None
  others = cl.get_nodes_with_value(checklist, keyprop, id)
  if len(others) == 1:
    return others[0]
  elif len(others) == 0:
    return None
  else:
    print("** This seems like a bug.  Key %s is ambiguous." % id)
    return None

def eq_partner(node, al):
  ar = al.get(node)
  if ar and ar.relation == rel.eq:
    return ar.cod
  else:
    return None

