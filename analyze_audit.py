import json, collections
with open('audit_gold_standard_results.json') as f:
    data = json.load(f)
upgrade = data['summary']['needs_upgrade_papers']
ic = collections.Counter()
for p in upgrade:
    for i in p['issues']:
        ic[i] += 1
print("Issue frequency across all 60 deficient papers:")
for issue, count in ic.most_common():
    print(f"  {count:3d}x  {issue}")
