#include "WAQbasedFairSubstring.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

WAQbasedFairSubstring::WAQbasedFairSubstring(const vector<string>& strings) 
    : inputStrings(strings) {
}

WAQbasedFairSubstring::~WAQbasedFairSubstring() {
    for (auto& levelPair : levelTrees) {
        for (auto tree : levelPair.second) {
            delete tree;
        }
    }
}

size_t WAQbasedFairSubstring::calculateTreeCountForLevel(int level) const {
    size_t groupSize = getGroupSizeForLevel(level);
    return inputStrings.size() / groupSize;
}

size_t WAQbasedFairSubstring::getGroupSizeForLevel(int level) const {
    return 1 << (level + 1);
}

int WAQbasedFairSubstring::getMaxLevelCount() const {
    if (inputStrings.empty()) return 0;
    
    return static_cast<int>(log2(inputStrings.size()));
}

void WAQbasedFairSubstring::cleanupLevelsBelowCurrent(int currentLevel) {
    for (auto it = levelTrees.begin(); it != levelTrees.end(); ) {
        if (it->first < currentLevel) {
            cout << "Cleaning up level " << (it->first + 1) << " trees..." << endl;
            
            for (auto tree : it->second) {
                delete tree;
            }
            
            it = levelTrees.erase(it);
        } else {
            ++it;
        }
    }
}

void WAQbasedFairSubstring::processLevel(int level) {
    if (level < 0) {
        cout << "Error: Invalid level " << level << endl;
        return;
    }
    
    if (level > 0 && levelTrees.find(level-1) == levelTrees.end()) {
        cout << "Error: Level " << (level-1) << " not processed yet. Process it first." << endl;
        return;
    }
    
    size_t groupSize = getGroupSizeForLevel(level);
    
    if (inputStrings.size() < groupSize) {
        cout << "Error: Not enough strings to process level " << level 
             << ". Need at least " << groupSize << " strings." << endl;
        return;
    }
    
    int treeCount = inputStrings.size() / groupSize;
    cout << "Processing " << treeCount << " string groups for level " << (level+1) << "..." << endl;
    
    vector<GeneralizedSuffixTree*> currentLevelTrees;
    
    // Process first level - pair of strings
    if (level == 0) {
        for (size_t i = 0; i < inputStrings.size(); i += 2) {
            if (i + 1 >= inputStrings.size()) break;
            
            cout << "Building level 1 GST for strings " << i << " and " << (i+1) << "..." << endl;
            
            GeneralizedSuffixTree* tree = new GeneralizedSuffixTree(
                inputStrings[i], 
                inputStrings[i+1]
            );
            
            currentLevelTrees.push_back(tree);
            
            cout << "Finished building level 1 GST for pair " << (i/2) << endl;
        }
    }
    // Process higher levels - using vector constructor and passing subtrees
    else {
        vector<GeneralizedSuffixTree*>& previousLevelTrees = levelTrees[level-1];
        
        for (size_t i = 0; i < inputStrings.size(); i += groupSize) {
            if (i + groupSize - 1 >= inputStrings.size()) {
                cout << "Skipping incomplete group at index " << i << endl;
                break;
            }
            
            cout << "Building level " << (level+1) << " GST for strings " 
                 << i << " to " << (i + groupSize - 1) << "..." << endl;
            
            vector<string> groupStrings;
            for (size_t j = 0; j < groupSize && (i + j) < inputStrings.size(); j++) {
                groupStrings.push_back(inputStrings[i + j]);
            }
            
            int childGroupSize = groupSize/2;
            size_t leftChildIndex = (i/childGroupSize)/2;
            size_t rightChildIndex = leftChildIndex + 1;
            
            if (leftChildIndex >= previousLevelTrees.size() || 
                rightChildIndex >= previousLevelTrees.size()) {
                cout << "Error: Invalid child tree indices: " << leftChildIndex 
                     << " and " << rightChildIndex << endl;
                continue;
            }
            
            GeneralizedSuffixTree* leftSubtree = previousLevelTrees[leftChildIndex];
            GeneralizedSuffixTree* rightSubtree = previousLevelTrees[rightChildIndex];
            
            cout << "Using child trees from indices " << leftChildIndex << " and " 
                 << rightChildIndex << " from level " << level << endl;
            
            GeneralizedSuffixTree* tree = new GeneralizedSuffixTree(groupStrings, leftSubtree, rightSubtree);
            
            currentLevelTrees.push_back(tree);
            
            cout << "Finished building level " << (level+1) << " GST for group " << (i/groupSize) << endl;
        }
    }
    
    levelTrees[level] = currentLevelTrees;
    
    cout << "Total number of level " << (level+1) << " GSTs built: " << currentLevelTrees.size() << endl;
}

void WAQbasedFairSubstring::traverseHighestLevelTree() {
    int highestLevel = -1;
    for (const auto& levelPair : levelTrees) {
        highestLevel = max(highestLevel, levelPair.first);
    }
    
    if (highestLevel == -1 || levelTrees[highestLevel].empty()) {
        cout << "No trees available to traverse." << endl;
        return;
    }
    
    cout << "\n--- Traversing highest level tree (Level " << (highestLevel+1) << ") ---\n" << endl;
    
    GeneralizedSuffixTree* highestTree = levelTrees[highestLevel][0];
    highestTree->traverseDFS();
}

void WAQbasedFairSubstring::processAllLevels() {
    int maxLevel = getMaxLevelCount();
    
    cout << "Processing up to " << (maxLevel+1) << " levels based on input size..." << endl;
    
    for (int level = 0; level <= maxLevel; level++) {
        size_t groupSize = getGroupSizeForLevel(level);
        
        if (inputStrings.size() < groupSize) {
            cout << "Stopping at level " << (level) << " as there are not enough strings for level " 
                 << (level+1) << endl;
            break;
        }
        
        cout << "\n--- Processing Level " << (level+1) << " ---\n" << endl;
        processLevel(level);
        
        if (level > 0) {
            cout << "\n--- Memory Cleanup after Level " << (level+1) << " ---\n" << endl;
            cleanupLevelsBelowCurrent(level);
        }
    }
}

size_t WAQbasedFairSubstring::getTreeCount(int level) const {
    auto it = levelTrees.find(level);
    if (it != levelTrees.end()) {
        return it->second.size();
    }
    return 0;
}

GeneralizedSuffixTree* WAQbasedFairSubstring::getTree(int level, size_t index) {
    auto it = levelTrees.find(level);
    if (it != levelTrees.end() && index < it->second.size()) {
        return it->second[index];
    }
    return nullptr;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }
    
    string inputFileName = argv[1];
    string outputFileName = argv[2];
    
    ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        cout << "Error: Could not open input file " << inputFileName << endl;
        return 1;
    }
    
    ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        cout << "Error: Could not open output file " << outputFileName << endl;
        return 1;
    }
    
    streambuf* coutBuffer = cout.rdbuf();
    cout.rdbuf(outputFile.rdbuf());
    
    vector<string> inputStrings;
    string line;
    
    cout << "Reading strings from file: " << inputFileName << endl;
    
    while (getline(inputFile, line) && !line.empty()) {
        inputStrings.push_back(line);
    }
    
    cout << "Read " << inputStrings.size() << " strings" << endl;
    
    if (!inputStrings.empty()) {
        size_t totalLength = 0;
        size_t minLength = inputStrings[0].length();
        size_t maxLength = inputStrings[0].length();
        
        for (const string& s : inputStrings) {
            totalLength += s.length();
            minLength = min(minLength, s.length());
            maxLength = max(maxLength, s.length());
        }
        
        cout << "String statistics:" << endl;
        cout << "- Average length: " << (double)totalLength / inputStrings.size() << endl;
        cout << "- Minimum length: " << minLength << endl;
        cout << "- Maximum length: " << maxLength << endl;
        cout << "- Total characters: " << totalLength << endl;
    }
    
    if (inputStrings.size() < 2) {
        cout << "Need at least 2 strings to process." << endl;
        cout.rdbuf(coutBuffer);
        return 1;
    }
    
    WAQbasedFairSubstring processor(inputStrings);
    
    processor.processAllLevels();
    processor.traverseHighestLevelTree();
    
    cout << "\n=== Summary ===" << endl;
    int maxLevel = processor.getMaxLevelCount();
    for (int level = 0; level <= maxLevel; level++) {
        size_t treeCount = processor.getTreeCount(level);
        if (treeCount > 0) {
            cout << "Level " << (level+1) << " trees: " << treeCount << endl;
        }
    }
    
    cout.rdbuf(coutBuffer);
    
    cout << "Processing completed. Results written to " << outputFileName << endl;
    cout << "Processed " << inputStrings.size() << " strings" << endl;
    
    return 0;
}