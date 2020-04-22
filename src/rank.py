# Thanks to Open Tree

rank_names = ["domain",
		          "superkingdom",
		          "kingdom",
		          "subkingdom",
              "division",       # h2007
		          "infrakingdom",		# worms
		          "superphylum",
		          "phylum",
		          "subphylum",
		          "infraphylum",    # worms
		          "subdivision",    # worms
		          "superclass",
		          "class",
		          "subclass",
		          "infraclass",
              "subterclass",    # worms Colobognatha
              "cohort",         # NCBI Polyneoptera
              "subcohort",      # NCBI
		          "superorder",
		          "order",
		          "suborder",
		          "infraorder",
		          "parvorder",
		          "section",				# worms
		          "subsection",			# worms
		          "superfamily",
		          "family",
		          "subfamily",
		          "supertribe",			# worms
		          "tribe",
		          "subtribe",
		          "genus",
		          "subgenus",
		          "species group",
		          "species subgroup",
		          "species",
		          "infraspecificname",
		          "subspecies",
              "natio",          # worms
		          "variety",
		          "varietas",
		          "subvariety",
		          "form",           # 2016 GBIF
		          "forma",
		          "subform",
              "cluster",
              ]

name_to_rank = {}
rank_to_name = {}

i = 100
for rankname in rank_names:
  name_to_rank[rankname] = i
  rank_to_name[i] = rankname
  i += 100

