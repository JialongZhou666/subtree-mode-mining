# fair-substring-mining
## Compilation Guide
### Baseline Methods
#### Baseline 1
```
g++ -o fair_substring_baseline1 baseline1.cpp
```

#### Baseline 2
```
g++ -o fair_substring_baseline2 baseline2.cpp
```

### WAQ-based Method
```
g++ -o fair_substring_waq WAQbasedFairSubstring.cpp GeneralizedSuffixTree.cpp Node.cpp ActivePoint.cpp
```

### LCA-based Method
```
g++ -o fair_substring_lca LCAbasedFairSubstring.cpp GeneralizedSuffixTree.cpp rmq-offline.cpp Node.cpp ActivePoint.cpp
```
