# fair-substring-mining
## Compilation Guide
### Baseline Method
```
g++ -o fair_substring_baseline baseline.cpp
```

### WAQ Method
```
g++ -o fair_substring_waq WAQbasedFairSubstring.cpp GeneralizedSuffixTree.cpp Node.cpp ActivePoint.cpp
```

### LCA Method
```
g++ -o fair_substring_lca LCAbasedFairSubstring.cpp GeneralizedSuffixTree.cpp rmq-offline.cpp Node.cpp ActivePoint.cpp
```
