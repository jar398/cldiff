"""
 Makes a subset of a checklist based on one subtree of a taxonmoy.

 python3 subset_dwc.py checklist root_id taxonomy out_dir
   'checklist' and 'taxonomy' are dwca unzip directories (can be the same)
      (must use the same taxon ids!)
   'root_id' is the taxonID from which to start
   'out_dir' will contain the resulting .tsv file

"""

import sys, os, csv

def main(checkdir, taxdir, root_id, outdir):
  assert os.path.exists(checkdir)
  assert os.path.exists(taxdir)
  assert os.path.exists(outdir)
  topo = read_topology(taxdir)
  print ("Nodes with children:", len(topo))
  all = closure(topo, root_id)
  print ("Nodes in subset:", len(all))
  write_subset(checkdir, root_id, all, outdir)

def write_subset(checkdir, root_id, all, outdir):
  with open(get_tnu_path(checkdir), "r") as infile:
    reader = csv.reader(infile, delimiter="\t", quotechar='\a', quoting=csv.QUOTE_NONE)
    head = next(reader)
    tid_column = head.index("taxonID") 
    aid_column = None
    if "acceptedNameUsageID" in head:
      aid_column = head.index("acceptedNameUsageID")

    with open(os.path.join(outdir, "taxon.tsv"), "w") as outfile:
      writer = csv.writer(outfile, delimiter="\t", quotechar='\a', quoting=csv.QUOTE_NONE)
      writer.writerow(head)
      for row in reader:
        tid = row[tid_column]
        if tid in all: 
          writer.writerow(row)
        elif aid_column:
          aid = row[aid_column]
          if aid in all:
            writer.writerow(row)

def closure(topo, root_id):
  all = {}
  empty = []
  def descend(id):
    all[id] = True
    for child in topo.get(id, empty):
      descend(child)
  descend(root_id)
  return all

def read_topology(taxdir):
  children = {}
  with open(get_tnu_path(taxdir), "r") as infile:
    reader = csv.reader(infile, delimiter="\t", quotechar='\a', quoting=csv.QUOTE_NONE)
    head = next(reader)
    tid_column = head.index("taxonID") 
    pid_column = head.index("parentNameUsageID")
    for row in reader:
      child_id = row[tid_column]
      parent_id = row[pid_column]
      childs = children.get(parent_id, None)
      if childs:
        childs.append(child_id)
      else:
        children[parent_id] = [child_id]
  return children

def get_tnu_path(dwca_dir):
  for name in ["taxon.tsv",
               "Taxon.tsv",
               "taxon.txt",
               "Taxon.txt"]:
    path = os.path.join(dwca_dir, name)
    if os.path.exists(path):
      return path
  raise valueError("cannot find TNU file", dwca_dir)

main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
