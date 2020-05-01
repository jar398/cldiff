

n_points = 10
n_sets = 1 << n_points

n_eq = 0
n_lt = 0
n_gt = 0
n_disjoint = 0
n_conflict = 0

n_total = 0
for set1 in range(0,n_sets):
  for set2 in range(0,n_sets):
    n_total += 1
    if set1 == set2:
      n_eq += 1
    if set1 & set2 == 0:
      n_disjoint += 1
    if set1 & set2 == set1:
      n_lt += 1
    if set1 & set2 == set2:
      n_gt += 1
    else:
      n_conflict += 1

print("=:  %s" % (n_eq / n_total))
print("<:  %s" % (n_lt / n_total))
print(">:  %s" % (n_gt / n_total))
print("><: %s" % (n_conflict / n_total))
print("!:  %s" % (n_disjoint / n_total))

