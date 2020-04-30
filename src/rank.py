# Thanks to Open Tree

# Height counts the number of groupings, bottom up, starting at 0 =
# fewest groups (everything is in a different group).

# Depth counts the number of bifurcations, starting at 0 = fewest
# bifurcations (everything is in one group).

# Approximately 50 currently named ranks

rank_configuration = [
  [[],
   ["individual"]],  # added by jar for fun

  [["no rank",       # NCBI
    "population"],     
   ["cluster",      # added by jar, for SILVA ?
    "subform",
    "forma",
    "form"]],           # 2016 GBIF

  [["subvariety",
    "varietas"],
   ["variety",
    "natio"]],          # worms

  [["subspecies",
    "infraspecificname"],
   ["species",
    "species subgroup",
    "species group"]],

  [["subgenus"],
   ["genus",
    "subtribe",
    "tribe",
    "supertribe"]],     # worms

  [["subfamily"],
   ["family",
    "superfamily",
    "subsection",     # worms
    "section"]],        # worms

  [["parvorder",
    "infraorder",
    "suborder"],
   ["order",
    "superorder",
    "subcohort",      # NCBI
    "cohort"]],         # NCBI Polyneoptera

  [["subterclass",    # worms Colobognatha
    "infraclass",
    "subclass"],
   ["class",
    "superclass"]],

  [["subdivision",     # worms
    "infraphylum",     # worms
    "subphylum"],
   ["phylum",
    "superphylum"]],

  [["infrakingdom",     # worms
    "division",         # h2007
    "subkingdom"],
   ["kingdom",
    "superkingdom",
    "domain",
    "root"]],

  [[], ["anything"]]    # added by JAR for fun
  ]

def process_ranks(groups):
  rank_to_height = {}
  height_to_name = {}
  max_depth = 0
  def associate(name, i):
    rank_to_height[name] = i
    height_to_name[i] = name
  linnean_increment = 100000
  foo = 10
  sublinnean_increment = int(linnean_increment / foo)
  height = 0
  for (colinnea, linnea) in groups:
    if len(colinnea) > foo/2:
      print("** Too many subgroup levels")
    if len(linnea) > foo/2:
      print("** Too many supergroup levels")
    print("height(%s) = %s" % (linnea[0], height))
    down = height
    for name in colinnea:
      down -= sublinnean_increment
      associate(name, down)
    up = height
    for name in linnea:
      associate(name, up)
      max_height = up
      up += sublinnean_increment
    height += linnean_increment
  return (rank_to_height, height_to_name, max_height)

(rank_to_height_table, height_to_name_table, max_height) = \
  process_ranks(rank_configuration)

max_depth = max_height

def rank_to_height(name):
  return rank_to_height_table.get(name)

def height_to_depth(height):
  return max_height - height

def depth_to_height(height):
  return max_depth - height

def rank_to_depth(name):
  height = rank_to_height_table.get(name)
  if height:
    return height_to_depth(height)
  return None

for rank in ["subspecies", "individual", "anything"]:
  print("height %s -> %s" % (rank, rank_to_height(rank)))
  print("depth %s -> %s" % (rank, rank_to_depth(rank)))

individual_height = 0
anything_depth = 0

def next_higher(height):        # more phylumward
  return height + 1

root_height = rank_to_height("root")
