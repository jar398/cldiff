"""
 Makes a subset of a checklist based on one subtree of a taxonmoy.

 python3 subset_dwc.py [--taxonomy tax_dwc] source_dwc out_dwc

 Assumption: every accepted record has a taxonID
"""

debug = False

import sys, os, csv, argparse

def main(checklist, tax_path, root_id, outpath):
  (topo, synonyms) = read_topology(tax_path)
  print ("Nodes with children: %s" % len(topo))
  print ("Nodes with synonyms: %s" % len(synonyms))
  all = closure(topo, synonyms, root_id)
  print ("Nodes in subset: %s" % len(all))
  write_subset(checklist, root_id, all, topo, outpath)

def write_subset(checklist, root_id, all, topo, outpath):

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
        if row[sid_column] != "doubtful":
          row = clean(row, tid_column, pid_column, aid_column, sid_column, topo)
          tid = row[tid_column]
          if tid in all:
            writer.writerow(row)

# Transitive closure of accepted records

def closure(topo, synonyms, root_id):
  all = {}
  empty = []
  def descend(id):
    if id in all:
      pass
    else:
      all[id] = True
      for child in topo.get(id, empty):
        descend(child)
      for syn in synonyms.get(id, empty):
        descend(syn)
  descend(root_id)
  return all

def read_topology(tax_path):
  children = {}
  synonyms = {}
  (delimiter, quotechar, mode) = csv_parameters(tax_path)
  with open(tax_path, "r") as infile:
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    head = next(reader)

    tid_column = head.index("taxonID") 
    pid_column = head.index("parentNameUsageID")
    aid_column = head.index("acceptedNameUsageID")

    if tid_column == None:      # usually 0
      print("** No taxonID column found")
    if pid_column == None:
      print("** No taxonID column found")
    if aid_column == None:
      print("** No acceptedNameUsageID column found")

    counter = 0
    for row in reader:
      counter += 1
      tid = row[tid_column]
      parent_id = row[pid_column]
      accepted_id = row[aid_column]

      if parent_id != '':
        childs = children.get(parent_id, None)
        if childs:
          childs.append(tid)
        else:
          children[parent_id] = [tid]
      if accepted_id != '':
        syns = synonyms.get(accepted_id, None)
        if syns:
          syns.append(tid)
        else:
          synonyms[accepted_id] = [tid]
    print("Nodes: %s" % counter)

  return (children, synonyms)

def clean(row, tid_column, pid_column, aid_column, sid_column, topo):
  tid = row[tid_column]
  parent_id = row[pid_column]
  accepted_id = row[aid_column]
  status = row[sid_column]

  # Clean up this record
  if parent_id == tid:
    row[pid_column] = ''
    parent_id = ''
  if accepted_id == tid:
    row[aid_column] = ''
    accepted_id = ''

  if parent_id != '' and accepted_id != '':
    # if it has children, it has to be accepted
    if status == "synonym" and not topo.get(tid):
      row[pid_column] = ''
      parent_id = ''
    else:
      row[aid_column] = ''
      accepted_id = ''
  elif parent_id != '' and accepted_id == '':
    if status == "synonym":
      print("** Taxonomic status is synonym, but no accepted node: %s" % tid)
      row[status] = "error"

  return row

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
