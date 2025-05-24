# fair-substring-mining
## Compilation Guide
### Baseline Methods
```
g++ -o fair_substring_baseline1 baseline1.cpp
g++ -o fair_substring_baseline2 baseline2.cpp
```

### WAQ Method
```
g++ -o fair_substring_waq WAQbasedFairSubstring.cpp GeneralizedSuffixTree.cpp Node.cpp ActivePoint.cpp
```

### LCA Method
```
g++ -o fair_substring_lca LCAbasedFairSubstring.cpp GeneralizedSuffixTree.cpp rmq-offline.cpp Node.cpp ActivePoint.cpp
```
