#include "LCAbasedFairSubstring.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

LCAbasedFairSubstring::LCAbasedFairSubstring(const vector<string>& strings) 
    : inputStrings(strings) {
}

LCAbasedFairSubstring::~LCAbasedFairSubstring() {
    for (auto& levelPair : levelTrees) {
        for (auto tree : levelPair.second) {
            delete tree;
        }
    }
}

size_t LCAbasedFairSubstring::calculateTreeCountForLevel(int level) const {
    size_t groupSize = getGroupSizeForLevel(level);
    return inputStrings.size() / groupSize;
}

size_t LCAbasedFairSubstring::getGroupSizeForLevel(int level) const {
    return 1 << (level + 1);
}

int LCAbasedFairSubstring::getMaxLevelCount() const {
    if (inputStrings.empty()) return 0;
    
    return static_cast<int>(log2(inputStrings.size()));
}

void LCAbasedFairSubstring::cleanupLevelsBelowCurrent(int currentLevel) {
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

void LCAbasedFairSubstring::processLevel(int level) {  
    size_t groupSize = getGroupSizeForLevel(level);
    
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
            cout << "Building level " << (level+1) << " GST for strings " 
                 << i << " to " << (i + groupSize - 1) << "..." << endl;
            
            vector<string> groupStrings;
            for (size_t j = 0; j < groupSize && (i + j) < inputStrings.size(); j++) {
                groupStrings.push_back(inputStrings[i + j]);
            }
            
            int childGroupSize = groupSize/2;
            size_t leftChildIndex = (i/childGroupSize)/2;
            size_t rightChildIndex = leftChildIndex + 1;
            
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

void LCAbasedFairSubstring::processAllLevels() {
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

size_t LCAbasedFairSubstring::getTreeCount(int level) const {
    auto it = levelTrees.find(level);
    if (it != levelTrees.end()) {
        return it->second.size();
    }
    return 0;
}

GeneralizedSuffixTree* LCAbasedFairSubstring::getTree(int level, size_t index) {
    auto it = levelTrees.find(level);
    if (it != levelTrees.end() && index < it->second.size()) {
        return it->second[index];
    }
    return nullptr;
}

void LCAbasedFairSubstring::traverseHighestLevelTree() {
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
        
    }
    
    if (inputStrings.size() < 2) {
        cout << "Need at least 2 strings to process." << endl;
        cout.rdbuf(coutBuffer);
        return 1;
    }
    
    LCAbasedFairSubstring processor(inputStrings);
    
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