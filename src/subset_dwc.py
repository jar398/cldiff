"""
 Makes a subset of a checklist based on one subtree of a taxonmoy.

 python3 subset_dwc.py [--taxonomy tax_dwc] source_dwc out_dwc

 Assumption: every accepted record has a taxonID
"""

debug = False

import sys, os, csv, argparse

def main(checklist, tax_path, root_id, outpath):
  topo = read_topology(tax_path)
  all = closure(topo, root_id)
  write_subset(checklist, root_id, all, topo, outpath)

def write_subset(checklist, root_id, all, topo, outpath):
  print("Writing subset to %s" % outpath, flush=True)

  (delimiter, quotechar, mode) = csv_parameters(checklist)
  with open(checklist, "r") as infile:
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    head = next(reader)

    tid_column = head.index("taxonID") 
    aid_column = head.index("acceptedNameUsageID")
    pid_column = head.index("parentNameUsageID")
    sid_column = head.index("taxonomicStatus")

    with open(outpath, "w") as outfile:
      (delimiter, quotechar, mode) = csv_parameters(outpath)
      writer = csv.writer(outfile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
      writer.writerow(head)
      for row in reader:
        row = clean(row, tid_column, pid_column, aid_column, sid_column, topo)
        tid = row[tid_column]
        if tid in all:
          writer.writerow(row)

# Transitive closure of accepted records

def closure(topo, root_id):
  print("Computing transitive closure starting from %s" % root_id, flush=True)
  all = {}
  empty = []
  def descend(id):
    if not id in all:
      all[id] = True
      if id in topo:
        (children, synonyms, _) = topo[id]
        for child in children:
          descend(child)
        for syn in synonyms:
          descend(syn)
  descend(root_id)
  print ("  Nodes in transitive closure: %s" % len(all))
  return all

def read_topology(tax_path):
  # Keyed by taxon id
  topo = {}
  (delimiter, quotechar, mode) = csv_parameters(tax_path)
  counter = 0
  with open(tax_path, "r") as infile:
    print("Scanning %s to obtain topology" % tax_path, flush=True)
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    head = next(reader)

    tid_column = head.index("taxonID") 
    pid_column = head.index("parentNameUsageID")
    aid_column = head.index("acceptedNameUsageID")
    sid_column = head.index("taxonomicStatus")

    if tid_column == None:      # usually 0
      print("** No taxonID column found")
    if pid_column == None:
      print("** No taxonID column found")
    if aid_column == None:
      print("** No acceptedNameUsageID column found")
    if sid_column == None:
      print("** No taxonomicStatus column found")

    for row in reader:
      counter += 1
      tid = row[tid_column]
      parent_id = row[pid_column]
      accepted_id = row[aid_column]
      status = row[sid_column]
      is_syn = is_synonym_status(status)
      get_topo_record(tid, topo)[2] = is_syn
      if is_syn:
        if accepted_id != '':
          (_, syns, _) = get_topo_record(accepted_id, topo)
          syns.append(tid)
      else:
        if parent_id != '':
          (children, _, _) = get_topo_record(parent_id, topo)
          children.append(tid)
    print("  %s nodes of which %s have children and/or synonyms" %
          (counter, len(topo)))

  return topo

def get_topo_record(tid, topo):
  record = topo.get(tid)
  if not record:
    record = [[], [], None]
    topo[tid] = record
  return record

def clean(row, tid_column, pid_column, aid_column, sid_column, topo):
  tid = row[tid_column]
  status = row[sid_column]
  accepted_id = row[aid_column]
  parent_id = row[pid_column]

  # Clean up this record

  if parent_id == tid:
    parent_id = ''
  if accepted_id == tid:
    accepted_id = ''

  # Compare with checklist.validate

  record = topo.get(tid)
  if record:
    (children, synonyms, is_syn) = record
    assert is_syn == is_synonym_status(status)
  else:
    children = []
    synonyms = []
    is_syn = is_synonym_status(status)

  if is_syn:
    if len(children) > 0:
      print("** Synonym %s (%s) has children" % (tid, status))
    if len(synonyms) > 0:
      print("** Synonym %s (%s) has synonyms" % (tid, status))
    if accepted_id == '':
      print("** Synonym %s (%s) has no accepted id" % (tid, status))
      if parent_id != '':
        print("** Synonym %s (%s) has a parent %s" % (tid, status, parent_id))
    else:
      (_, _, is_syn_) = topo.get(accepted_id)
      if is_syn_ == None:
        print("** Synonym %s accepted id %s -> nowhere" % (tid, accepted_id))
      elif is_syn_:
        print("** Synonym %s accepted id %s -> synonym" % (tid, accepted_id))
      if parent_id != '':
        # This occurs a lot in GBIF, usually (always?) points to accepted's parent
        # print("** Synonym %s has a parent %s" % (tid, parent_id))
        parent_id = ''
  else:
    if accepted_id != '':
      print("** Accepted node %s (%s) has accepted id %s" % (tid, status, accepted_id))
      accepted_id = ''
    if parent_id == '':
      pass    # this is OK, it would be a root
    else:
      (_, _, is_syn_) = topo.get(parent_id)
      if is_syn_ == None:
        print("** Accepted %s parent id %s -> nowhere" % (tid, parent_id))
      elif is_syn_:
        print("** Accepted %s has parent id %s -> synonym" % (tid, parent_id))
    for id in children:
      if id in topo:
        (_, _, is_syn) = topo[id]
        if is_syn:
          print("** Accepted %s (%s) has child %s which is a synonym" % (tid, status, id))
    for id in synonyms:
      if id in topo:
        (_, _, is_syn) = topo[id]
        if not is_syn:
          print("** Accepted %s (%s) has synonym %s which is accepted" % (tid, status, id))

  row[pid_column] = parent_id
  row[aid_column] = accepted_id

  return row

def is_synonym_status(status):
  return ("synonym" in status) or (status == "misapplied")

def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

# main(checklist, taxonomy, root_id, outfile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--taxonomy', help="""taxonomy from which to extract
                            hierarchy; defaults to source""")
  parser.add_argument('source', help='taxonomy or checklist from which to extract subset')
  parser.add_argument('id', help="taxon id of subset's root")
  parser.add_argument('--out', help='where to store the subset')
  args = parser.parse_args()
  main(args.source, args.taxonomy or args.source, args.id, args.out)
