#ifndef WAQ_BASED_FAIR_SUBSTRING_H
#define WAQ_BASED_FAIR_SUBSTRING_H

#include "GeneralizedSuffixTree.h"
#include <vector>
#include <string>
#include <map>

class WAQbasedFairSubstring {
private:
    std::vector<std::string> inputStrings;
    
    std::map<int, std::vector<GeneralizedSuffixTree*> > levelTrees;

public:
    WAQbasedFairSubstring(const std::vector<std::string>& strings);
    ~WAQbasedFairSubstring();
    
    size_t calculateTreeCountForLevel(int level) const;
    
    size_t getGroupSizeForLevel(int level) const;
    
    int getMaxLevelCount() const;

    void processLevel(int level);
    
    void processAllLevels();
    
    size_t getTreeCount(int level) const;
    
    GeneralizedSuffixTree* getTree(int level, size_t index);

    void cleanupLevelsBelowCurrent(int currentLevel);
    
    void traverseHighestLevelTree();
};

#endif // WAQ_BASED_FAIR_SUBSTRING_H