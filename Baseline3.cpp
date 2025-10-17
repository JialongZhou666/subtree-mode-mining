#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <stack>
#include <map>
#include <thread>
#include <atomic>
#include <mutex>
#include <fstream>
#include <sstream>
#include <cmath>

// Forward declarations for divsufsort
extern "C" {
    int32_t divsufsort(const uint8_t *T, int32_t *SA, int32_t n);
    int64_t divsufsort64(const uint8_t *T, int64_t *SA, int64_t n);
}

using namespace std;

using INT = int64_t;

void computeLCP(const string& text, const vector<INT>& sa, const vector<INT>& isa, vector<INT>& lcp) {
    INT n = text.length();
    lcp.assign(n, 0);
    INT k = 0;
    for (INT i = 0; i < n; ++i) {
        if (isa[i] == n - 1) {
            k = 0;
            continue;
        }
        INT j = sa[isa[i] + 1];
        while (i + k < n && j + k < n && text[i + k] == text[j + k]) {
            k++;
        }
        lcp[isa[i]] = k;
        if (k > 0) {
            k--;
        }
    }
}

class LCP_RMQ {
private:
    vector<vector<INT>> st;
    vector<INT> log_table;
    INT n;

public:
    LCP_RMQ() : n(0) {}

    LCP_RMQ(const vector<INT>& arr) {
        n = arr.size();
        if (n == 0) return;
        log_table.resize(n + 1);
        log_table[1] = 0;
        for (INT i = 2; i <= n; i++) {
            log_table[i] = log_table[i / 2] + 1;
        }

        INT k = log_table[n];
        st.resize(n, vector<INT>(k + 1));
        for (INT i = 0; i < n; i++) {
            st[i][0] = arr[i];
        }

        for (INT j = 1; j <= k; j++) {
            for (INT i = 0; i + (1 << j) <= n; i++) {
                st[i][j] = min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    INT query(INT l, INT r) {
        if (l > r) swap(l, r);
        if (l >= n || r >= n || l < 0 || r < 0) return 0;
        INT j = log_table[r - l + 1];
        return min(st[l][j], st[r - (1 << j) + 1][j]);
    }
};

class EnhancedSuffixArray {
public:
    vector<string> originalStrings;
    string combinedText;
    vector<INT> SA;
    vector<INT> ISA;
    vector<INT> LCP;
    vector<INT> max_freq;
    vector<INT> min_freq;
    vector<INT> textPosToStringId;
    map<pair<INT, INT>, INT> leafCache;
    LCP_RMQ rmq;

    EnhancedSuffixArray(const string& singleString) {
        originalStrings.push_back(singleString);
        buildAndProcess(1);
    }

    EnhancedSuffixArray(const string& str1, const string& str2) {
        originalStrings.push_back(str1);
        originalStrings.push_back(str2);
        buildAndProcess(2);
    }

    EnhancedSuffixArray(const vector<string>& inputStrings, EnhancedSuffixArray* leftSubtree, EnhancedSuffixArray* rightSubtree)
        : originalStrings(inputStrings) {
        
        buildCombinedText();
        INT n = combinedText.length();

        SA.resize(n);
        ISA.resize(n);
        divsufsort64(reinterpret_cast<const uint8_t*>(combinedText.c_str()), SA.data(), n);
        for(INT i=0; i<n; ++i) ISA[SA[i]] = i;

        computeLCP(combinedText, SA, ISA, LCP);
        rmq = LCP_RMQ(LCP);

        updateFrequencies_Merge(leftSubtree, rightSubtree);
    }

private:
    char getTerminatorChar(size_t index) {
        return (char)(index + 1);
    }

    void buildCombinedText() {
        stringstream ss;
        textPosToStringId.clear();
        for (size_t i = 0; i < originalStrings.size(); ++i) {
            ss << originalStrings[i];
            ss << getTerminatorChar(i);
            for(size_t j = 0; j < originalStrings[i].length() + 1; ++j) {
                textPosToStringId.push_back(i);
            }
        }
        combinedText = ss.str();
    }
    
    void buildAndProcess(int num_groups) {
        buildCombinedText();
        INT n = combinedText.length();

        SA.resize(n);
        ISA.resize(n);
        divsufsort64(reinterpret_cast<const uint8_t*>(combinedText.c_str()), SA.data(), n);
        for(INT i=0; i<n; ++i) ISA[SA[i]] = i;

        computeLCP(combinedText, SA, ISA, LCP);
        rmq = LCP_RMQ(LCP);

        if (num_groups == 1) updateFrequencies_Single();
        else if (num_groups == 2) updateFrequencies_TwoGroups();
    }

    void updateFrequencies_Single() {
        INT n = SA.size();
        max_freq.assign(n, 0);
        min_freq.assign(n, 0);
        
        vector<INT> child(n, -1);
        stack<INT> s;
        s.push(0);

        for (INT i = 1; i < n; ++i) {
            INT last_child = -1;
            while (LCP[i] < (s.top() == 0 ? 0 : LCP[s.top()])) {
                last_child = s.top();
                s.pop();
                if (child[s.top()] != -1) processNode_Base(child[s.top()]);
                processNode_Base(last_child);
            }
            if (LCP[i] > (s.top() == 0 ? 0 : LCP[s.top()])) {
                child[i] = last_child;
                s.push(i);
            }
        }
        while(!s.empty()){
            if (!s.empty() && child[s.top()] != -1) processNode_Base(child[s.top()]);
            s.pop();
        }
    }
    
    void updateFrequencies_TwoGroups() {
        INT n = SA.size();
        max_freq.assign(n, 0);
        min_freq.assign(n, 0);
        
        vector<INT> child(n, -1);
        stack<INT> s;
        s.push(0);

        for (INT i = 1; i < n; ++i) {
            INT last_child = -1;
            while (LCP[i] < (s.top() == 0 ? 0 : LCP[s.top()])) {
                last_child = s.top();
                s.pop();
                if (child[s.top()] != -1) processNode_Base(child[s.top()]);
                processNode_Base(last_child);
            }
            if (LCP[i] > (s.top() == 0 ? 0 : LCP[s.top()])) {
                child[i] = last_child;
                s.push(i);
            }
        }
        while(!s.empty()){
             if (!s.empty() && child[s.top()] != -1) processNode_Base(child[s.top()]);
            s.pop();
        }
    }

    void processNode_Base(INT sa_rank) {
         max_freq[sa_rank] = 1; 
         min_freq[sa_rank] = 1;
    }

    void buildLeafCache() {
        if (!leafCache.empty()) return;
        INT current_pos = 0;
        for (size_t i = 0; i < originalStrings.size(); ++i) {
            for (size_t j = 0; j < originalStrings[i].length(); ++j) {
                leafCache[{(INT)i, (INT)j}] = ISA[current_pos + j];
            }
            current_pos += originalStrings[i].length() + 1;
        }
    }

    void updateFrequencies_Merge(EnhancedSuffixArray* left, EnhancedSuffixArray* right) {
        left->buildLeafCache();
        right->buildLeafCache();
        
        INT n = SA.size();
        max_freq.assign(n, 0);
        min_freq.assign(n, 0);
    }
};

class ESAbasedFairSubstring {
private:
    vector<string> inputStrings;
    map<int, vector<EnhancedSuffixArray*>> levelESAs;

    void cleanupSpecificLevel(int level) {
        auto it = levelESAs.find(level);
        if (it != levelESAs.end()) {
            for (auto esa : it->second) {
                if (esa != nullptr) delete esa;
            }
            levelESAs.erase(it);
        }
    }
    
    vector<string> mergeStringLists(const vector<string>& left, const vector<string>& right) {
        vector<string> combined = left;
        combined.insert(combined.end(), right.begin(), right.end());
        return combined;
    }

    void processHigherLevelParallel(int level, vector<EnhancedSuffixArray*>& prevLevel, vector<EnhancedSuffixArray*>& currentLevel) {
        int completePairs = prevLevel.size() / 2;
        bool hasOdd = (prevLevel.size() % 2 == 1);
        int totalGroups = completePairs + (hasOdd ? 1 : 0);
        
        vector<bool> treeTransferred(prevLevel.size(), false);

        #pragma omp parallel for
        for (int group = 0; group < totalGroups; ++group) {
            if (group < completePairs) {
                size_t leftIdx = group * 2;
                size_t rightIdx = leftIdx + 1;
                EnhancedSuffixArray* leftSub = prevLevel[leftIdx];
                EnhancedSuffixArray* rightSub = prevLevel[rightIdx];
                if (leftSub && rightSub) {
                    vector<string> groupStrings = mergeStringLists(leftSub->originalStrings, rightSub->originalStrings);
                    currentLevel[group] = new EnhancedSuffixArray(groupStrings, leftSub, rightSub);
                }
            } else {
                size_t oddIdx = prevLevel.size() - 1;
                #pragma omp critical
                {
                    if (!treeTransferred[oddIdx]) {
                        currentLevel[group] = prevLevel[oddIdx];
                        treeTransferred[oddIdx] = true;
                        prevLevel[oddIdx] = nullptr;
                    }
                }
            }
        }
    }

    void processLevel(int level, int num_threads) {
        if (level == 0) {
            size_t totalStrings = inputStrings.size();
            int completePairs = totalStrings / 2;
            bool hasOdd = (totalStrings % 2 == 1);
            int totalTasks = completePairs + (hasOdd ? 1 : 0);
            
            vector<EnhancedSuffixArray*> currentLevel(totalTasks, nullptr);

            #pragma omp parallel for num_threads(num_threads)
            for (int task = 0; task < totalTasks; ++task) {
                if (task < completePairs) {
                    size_t leftIdx = task * 2;
                    size_t rightIdx = leftIdx + 1;
                    currentLevel[task] = new EnhancedSuffixArray(inputStrings[leftIdx], inputStrings[rightIdx]);
                } else {
                    currentLevel[task] = new EnhancedSuffixArray(inputStrings.back());
                }
            }
            levelESAs[level] = currentLevel;
        } else {
            vector<EnhancedSuffixArray*>& prevLevel = levelESAs[level - 1];
            int completePairs = prevLevel.size() / 2;
            bool hasOdd = (prevLevel.size() % 2 == 1);
            int totalGroups = completePairs + (hasOdd ? 1 : 0);
            vector<EnhancedSuffixArray*> currentLevel(totalGroups, nullptr);
            
            processHigherLevelParallel(level, prevLevel, currentLevel);
            
            levelESAs[level] = currentLevel;
        }

        if (level > 0) {
            cleanupSpecificLevel(level - 1);
        }
    }

public:
    ESAbasedFairSubstring(const vector<string>& strings) : inputStrings(strings) {}
    
    ~ESAbasedFairSubstring() {
        for (auto& levelPair : levelESAs) {
            for (auto esa : levelPair.second) {
                if (esa != nullptr) delete esa;
            }
        }
        levelESAs.clear();
    }

    void processAllLevels(int num_threads) {
        int theoreticalMaxLevel = 0;
        if (!inputStrings.empty()) {
            theoreticalMaxLevel = ceil(log2(inputStrings.size()));
        }
        
        for (int currentLevel = 0; currentLevel <= theoreticalMaxLevel; ++currentLevel) {
            processLevel(currentLevel, num_threads);
            auto it = levelESAs.find(currentLevel);
            if (it != levelESAs.end() && it->second.size() <= 1) {
                break;
            }
        }
    }
    
    EnhancedSuffixArray* getFinalESA() {
        if (levelESAs.empty()) return nullptr;
        int highestLevel = levelESAs.rbegin()->first;
        if (levelESAs[highestLevel].empty()) return nullptr;
        return levelESAs[highestLevel][0];
    }
    
    vector<string> findSubstringsWithinThreshold(int threshold) {
        vector<string> results;
        EnhancedSuffixArray* finalESA = getFinalESA();
        if (finalESA) {
            // Logic to traverse finalESA intervals would go here.
        }
        return results;
    }
};

vector<string> readStringsFromFile(const string& filename) {
    ifstream inFile(filename);
    if (!inFile.is_open()) {
        return {};
    }

    vector<string> strings;
    string currentLine;

    while (getline(inFile, currentLine)) {
        if (currentLine.empty()) continue;
        
        string processedLine;
        processedLine.reserve(currentLine.length());
        for (char c : currentLine) {
             if (c >= 32 && c <= 126) {
                processedLine += c;
            }
        }

        if (!processedLine.empty()) {
            strings.push_back(move(processedLine));
        }
    }
    strings.shrink_to_fit();
    return strings;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        return 1;
    }

    string inputFileName = argv[1];
    string outputFileName = argv[2];
    
    int num_threads = thread::hardware_concurrency();
    if (num_threads == 0) {
        num_threads = 1;
    }

    auto strings = readStringsFromFile(inputFileName);
    if (strings.empty()) {
        return 1;
    }

    ESAbasedFairSubstring* processor = new ESAbasedFairSubstring(strings);
    processor->processAllLevels(num_threads);

    int threshold = 1; 
    vector<string> fairSubstrings = processor->findSubstringsWithinThreshold(threshold);
    
    // Writing results to output file
    ofstream outFile(outputFileName);
    if(outFile.is_open()){
        for(const auto& str : fairSubstrings){
            outFile << str << "\n";
        }
    }

    delete processor;

    return 0;
}