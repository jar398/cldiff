"""
 Converts an NCBI dump to a DwCA TNU (taxon) file.
 Fold the scientific name (-> canonicalName) and authority (if it
 extends the scientific name) into the taxon record.

 python3 ncbi_to_dwca.py from to
   'from' is a directory containing .dmp files from ncbi
   'to' is a directory that will contain files for DwCA contents

 E.g.
  python3 ncbi_to_dwca.py ~/ncbi/2015-01-01/dump ~/ncbi/2015-01-01/dwca
"""

import sys, os, csv, argparse

def main(indir, outpath):
  assert os.path.exists(indir)
  accepteds = read_accepteds(os.path.join(indir, "nodes.dmp"))
  names = read_names(os.path.join(indir, "names.dmp"))
  merged = read_merged(os.path.join(indir, "merged.dmp"))
  (synonyms, scinames, authorities) = collate_names(names, accepteds)
  emit_dwc(accepteds, synonyms, scinames, authorities, merged, outpath)

def write_row(writer,
              taxonID, ncbi_id, parentNameUsageID, taxonRank,
              acceptedNameUsageID, scientificName, canonicalName,
              taxonomicStatus, nomenclaturalStatus):
  writer.writerow([taxonID, ncbi_id, parentNameUsageID, taxonRank,
                   acceptedNameUsageID, scientificName, canonicalName,
                   taxonomicStatus, nomenclaturalStatus])

def emit_dwc(accepteds, synonyms, scinames, authorities, merged, outpath):
  outdir = os.path.basename(outpath)
  if not os.path.isdir(outdir): os.mkdir(outdir)
  (delimiter, quotechar, mode) = csv_parameters(outpath)
  print ("Writing", outpath)
  with open(outpath, "w") as outfile:
    writer = csv.writer(outfile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    write_row(writer,
              "taxonID", "NCBI Taxonomy ID", "parentNameUsageID", "taxonRank",
              "acceptedNameUsageID", "scientificName", "canonicalName",
              "taxonomicStatus", "nomenclaturalStatus")
    for (id, parent_id, rank) in accepteds:
      sci = scinames.get(id, None)
      write_row(writer,
                id, id, parent_id, rank,
                None, authorities.get(id, None), sci,
                "accepted", None)
      if sci and "BOLD:" in sci:
        z = sci.split("BOLD:")
        write_row(writer,
                  id + ".BOLD", None, None, None,
                  id, None, "BOLD:" + z[-1],
                  "synonym", "BOLD id")
    for (id, text, kind, spin) in synonyms:
      if "BOLD:" in text:
        z = text.split("BOLD:")
        write_row(writer,
                  id + ".BOLD", None, None, None,
                  id, None, "BOLD:" + z[-1],
                  "synonym", "BOLD id")
      # synonym is a taxonomic status, not a nomenclatural status
      elif kind == "synonym": kind = None
      node_id = id + "." + str(spin)
      write_row(writer,
                node_id, None, None, None,
                id, None, text,
                "synonym", kind)
    for (old_id, new_id) in merged:
      canonical = "%s merged into %s" % (old_id, new_id)
      write_row(writer,
                old_id, old_id, None, None,
                new_id, None, canonical,
                "synonym", "merged id")

# Input: list of (id, text, kind, spin) from names.txt file
# Output: list of (id, text, kind, spin); 
#         dict: id -> text [canonical names];
#         dict: id -> text [authorities]

def collate_names(names, accepteds):
  keep = []
  scinames = {}
  for row in names:
    (id, text, kind, spin) = row
    if kind == "scientific name":
      scinames[id] = text
    else:
      keep.append(row)
  # Remove the canonical names, keep the rest
  names2 = keep
  keep = None #GC
  print (len(scinames), "canonicalNames (NCBI scientific names)")
  synonyms = []
  authorities = {}
  for row in names2:
    (id, text, kind, spin) = row
    if kind == "authority":
      probe = scinames.get(id, None)
      if probe and text.startswith(probe):
        authorities[id] = text
      else:
        synonyms.append(row)
    else:
      synonyms.append(row)
  # Remove authorities (= scientific names), keep the rest
  print (len(authorities), "scientificNames (NCBI authorities)")
  return (synonyms, scinames, authorities)

def read_accepteds(nodes_path):
  accepteds = []
  # Read the nodes file
  with open(nodes_path, "r") as infile:
    for row in csv.reader(infile,
                          delimiter="\t",
                          quotechar="\a",
                          quoting=csv.QUOTE_NONE):
      # tax_id, |, parent tax_id, |, rank, ... other stuff we don't use ...
      rank = row[4]
      if rank == "clade" or rank == "no rank":
        rank = None
      accepteds.append((row[0], row[2], rank))
  print (len(accepteds), "accepteds")
  return accepteds

def read_names(names_path):
  names = []
  # Read the names file
  with open(names_path, "r") as infile:
    # Depends on names being grouped by taxa
    previous_id = None
    spin = -1
    for row in csv.reader(infile,
                          delimiter="\t",
                          quotechar="\a",
                          quoting=csv.QUOTE_NONE):
      # tax_id, |, text, |, <unused>, |, name class, ... other stuff we don't use ...
      id = row[0]
      if str(id) != previous_id:
        spin = 1
        previous_id = id
      else:
        spin += 1
      names.append((id, row[2], row[6], spin))
  print (len(names), "names")
  return names

def read_merged(merged_path):
  merged = []
  # Read the merged file
  with open(merged_path, "r") as infile:
    for row in csv.reader(infile,
                          delimiter="\t",
                          quotechar="\a",
                          quoting=csv.QUOTE_NONE):
      # old_tax_id, |, new_tax_id
      merged.append((row[0], row[2]))
  print (len(merged), "merged")
  return merged

def csv_parameters(path):
  if path.endswith(".csv"):
    print("CSV")
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    print("TSV")
    return ("\t", "\a", csv.QUOTE_NONE)

# When invoked from command line:

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dump', help='directory containing taxdump files')
  parser.add_argument('--out', help='where to store the DwC version')
  args = parser.parse_args()
  main(args.dump, args.out)
