#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <stack>
#include <unordered_map>
#include <tuple>
#include <functional>
#include <divsufsort64.h>
#include <sdsl/rmq_support.hpp>
#include <sdsl/util.hpp>

using namespace std;
using namespace sdsl;

using INT = int64_t;

const char SEPARATOR = '$';

struct VirtualNode {
    INT leftBound;
    INT rightBound;
    INT lcpDepth;
    INT max_freq;
    INT max_color;
    INT leaf_count;

    VirtualNode(INT l, INT r, INT d);
};

struct VirtualNodeHash {
    size_t operator()(const tuple<INT, INT, INT>& t) const;
};

struct ColorMapping {
    INT singleLeft;
    INT singleRight;
    INT singleLCP;
    INT color;
    INT leafCount;
};

class PureESA {
private:
    string T;
    vector<INT> SA;
    int_vector<> LCP;
    rmq_succinct_sct<> rmq;
    size_t k;
    size_t n;
    vector<size_t> colorBoundaries;
    unordered_map<tuple<INT, INT, INT>, VirtualNode, VirtualNodeHash> nodeCache;

    void buildLCP();

public:
    explicit PureESA(const vector<string>& sequences);
    INT getColor(INT saIndex) const;
    INT getLCP(INT i, INT j) const;
    INT getLCPValue(INT i) const { return LCP[i]; }
    INT size() const { return n; }
    void enumerateInternalNodes(function<void(INT, INT, INT)> callback);
    VirtualNode& getOrCreateNode(INT left, INT right, INT lcp);
    INT findLCADepth(INT left, INT right) const;
    vector<INT> getColorLeaves(INT color) const;
};

class MultiColorTreeAlgorithm {
private:
    vector<string> inputStrings;
    PureESA* mainESA;
    vector<PureESA*> singleESAs;
    unordered_map<tuple<INT, INT, INT>, vector<ColorMapping>, VirtualNodeHash> nodeMappings;
    static int global_node_id;

    void buildMainESA();
    void buildSingleESAs();
    void establishIdentifierMappingsViaLCA_OptimalMerged();
    void performFrequencySelectionDFS();
    void selectMaxFrequencies();

public:
    explicit MultiColorTreeAlgorithm(const vector<string>& strings);
    ~MultiColorTreeAlgorithm();
    void runAlgorithm();
};

VirtualNode::VirtualNode(INT l, INT r, INT d)
    : leftBound(l), rightBound(r), lcpDepth(d),
      max_freq(0), max_color(-1), leaf_count(0) {}

size_t VirtualNodeHash::operator()(const tuple<INT, INT, INT>& t) const {
    auto h1 = hash<INT>{}(get<0>(t));
    auto h2 = hash<INT>{}(get<1>(t));
    auto h3 = hash<INT>{}(get<2>(t));
    return h1 ^ (h2 << 1) ^ (h3 << 2);
}

PureESA::PureESA(const vector<string>& sequences) : k(sequences.size()), n(0) {
    colorBoundaries.push_back(0);
    for (const auto& seq : sequences) {
        T += seq;
        T += SEPARATOR;
        colorBoundaries.push_back(T.length());
    }
    n = T.length();

    SA.resize(n);
    divsufsort64((sauchar_t*)T.c_str(), (saidx64_t*)SA.data(), n);

    buildLCP();
    
    sdsl::util::assign(rmq, sdsl::rmq_succinct_sct<>(&LCP));
}

void PureESA::buildLCP() {
    LCP.resize(n, 0);
    vector<INT> rank(n);
    for (INT i = 0; i < n; i++) rank[SA[i]] = i;

    INT h = 0;
    for (INT i = 0; i < n; i++) {
        if (rank[i] > 0) {
            INT j = SA[rank[i] - 1];
            while (i + h < n && j + h < n && T[i + h] == T[j + h]) h++;
            LCP[rank[i]] = h;
            if (h > 0) h--;
        }
    }
}

INT PureESA::getColor(INT saIndex) const {
    if (saIndex < 0 || saIndex >= n) return -1;
    INT pos = SA[saIndex];
    auto it = upper_bound(colorBoundaries.begin(), colorBoundaries.end(), pos);
    return (it == colorBoundaries.begin()) ? -1 : (it - colorBoundaries.begin() - 1);
}

INT PureESA::getLCP(INT i, INT j) const {
    if (i == j) return n - SA[i];
    if (i > j) swap(i, j);
    return LCP[rmq(i + 1, j)];
}

void PureESA::enumerateInternalNodes(function<void(INT, INT, INT)> callback) {
    struct StackEntry { INT lcp; INT left; };
    vector<StackEntry> stack;

    for (INT i = 0; i < n; i++) {
        INT currentLCP = (i < n - 1) ? LCP[i + 1] : 0;
        INT left = i;

        while (!stack.empty() && stack.back().lcp > currentLCP) {
            StackEntry top = stack.back();
            stack.pop_back();
            callback(top.left, i, top.lcp);
            left = top.left;
        }

        if (stack.empty() || stack.back().lcp < currentLCP) {
            stack.push_back({currentLCP, left});
        }
    }

    while (!stack.empty()) {
        StackEntry top = stack.back();
        stack.pop_back();
        callback(top.left, n - 1, top.lcp);
    }
}

VirtualNode& PureESA::getOrCreateNode(INT left, INT right, INT lcp) {
    auto key = make_tuple(left, right, lcp);
    auto it = nodeCache.find(key);
    if (it != nodeCache.end()) {
        return it->second;
    }
    nodeCache.emplace(key, VirtualNode(left, right, lcp));
    return nodeCache[key];
}

INT PureESA::findLCADepth(INT left, INT right) const {
    if (left > right) swap(left, right);
    if (left == right) return n - SA[left];
    return LCP[rmq(left + 1, right)];
}

vector<INT> PureESA::getColorLeaves(INT color) const {
    vector<INT> leaves;
    for (INT i = 0; i < n; i++) {
        if (getColor(i) == color) {
            leaves.push_back(i);
        }
    }
    return leaves;
}

int MultiColorTreeAlgorithm::global_node_id = 100000;

MultiColorTreeAlgorithm::MultiColorTreeAlgorithm(const vector<string>& strings)
    : inputStrings(strings), mainESA(nullptr) {
    singleESAs.resize(strings.size(), nullptr);
}

MultiColorTreeAlgorithm::~MultiColorTreeAlgorithm() {
    delete mainESA;
    for (auto* esa : singleESAs) delete esa;
}

void MultiColorTreeAlgorithm::buildMainESA() {
    mainESA = new PureESA(inputStrings);
}

void MultiColorTreeAlgorithm::buildSingleESAs() {
    for (size_t i = 0; i < inputStrings.size(); i++) {
        singleESAs[i] = new PureESA({inputStrings[i]});
    }
}

void MultiColorTreeAlgorithm::establishIdentifierMappingsViaLCA_OptimalMerged() {
    struct VirtualTreeNode {
        INT left, right, lcp;
        vector<INT> children;
        unordered_map<INT, vector<INT>> colorLeaves;
        INT parentLCP = 0;
    };
    
    INT n = mainESA->size();
    vector<VirtualTreeNode> nodes(n);
    stack<INT> st;
    
    nodes[0].left = 0;
    nodes[0].right = 0;
    nodes[0].lcp = 0;
    st.push(0);
    
    INT x = -1;
    
    for (INT i = 0; i < n; i++) {
        INT currentLCP = (i < n - 1) ? mainESA->getLCPValue(i+1) : 0;
        INT currentColor = mainESA->getColor(i);
        INT l = i;
        
        while (!st.empty() && nodes[st.top()].lcp > currentLCP) {
            x = st.top();
            st.pop();
            
            nodes[x].right = i - 1;
            
            if (!st.empty()) {
                for (auto& [color, leaves] : nodes[x].colorLeaves) {
                    auto& parentLeaves = nodes[st.top()].colorLeaves[color];
                    parentLeaves.insert(parentLeaves.end(), leaves.begin(), leaves.end());
                }
            }
            
            for (INT childIdx : nodes[x].children) {
                nodes[childIdx].parentLCP = nodes[x].lcp;
            }
            
            l = nodes[x].left;
            
            if (st.empty() || currentLCP <= nodes[st.top()].lcp) {
                if (!st.empty()) {
                    nodes[st.top()].children.push_back(x);
                }
                x = -1;
            }
        }
        
        if (!st.empty() && currentColor >= 0) {
            nodes[st.top()].colorLeaves[currentColor].push_back(i);
        }
        
        if (st.empty() || nodes[st.top()].lcp < currentLCP) {
            nodes[i].lcp = currentLCP;
            nodes[i].left = l;
            nodes[i].colorLeaves.clear();
            if (currentColor >= 0) {
                nodes[i].colorLeaves[currentColor].push_back(i);
            }
            st.push(i);
            if (~x) {
                nodes[i].children.push_back(x);
                x = -1;
            }
        }
    }
    
    while (!st.empty()) {
        x = st.top();
        st.pop();
        nodes[x].right = n - 1;
        
        for (INT childIdx : nodes[x].children) {
            nodes[childIdx].parentLCP = nodes[x].lcp;
        }
    }
    
    for (INT i = 0; i < n; i++) {
        if (nodes[i].lcp == 0) continue;
        
        auto mainKey = make_tuple(nodes[i].left, nodes[i].right, nodes[i].lcp);
        
        for (auto& [color, leaves] : nodes[i].colorLeaves) {
            if (leaves.empty()) continue;
            
            INT leftSA = *min_element(leaves.begin(), leaves.end());
            INT rightSA = *max_element(leaves.begin(), leaves.end());
            
            INT lcaDepth = singleESAs[color]->findLCADepth(leftSA, rightSA);
            
            ColorMapping mapping;
            mapping.singleLeft = leftSA;
            mapping.singleRight = rightSA;
            mapping.singleLCP = lcaDepth;
            mapping.color = color;
            mapping.leafCount = leaves.size();
            
            nodeMappings[mainKey].push_back(mapping);
        }
    }
}

void MultiColorTreeAlgorithm::selectMaxFrequencies() {
    for (auto& [mainKey, mappings] : nodeMappings) {
        auto [left, right, lcp] = mainKey;
        auto& node = mainESA->getOrCreateNode(left, right, lcp);
        
        INT maxFreq = 0;
        INT bestColor = -1;
        
        for (const auto& mapping : mappings) {
            if (mapping.leafCount > maxFreq) {
                maxFreq = mapping.leafCount;
                bestColor = mapping.color;
            }
        }
        
        node.max_freq = maxFreq;
        node.max_color = bestColor;
    }
}

void MultiColorTreeAlgorithm::performFrequencySelectionDFS() {
    selectMaxFrequencies();
}

void MultiColorTreeAlgorithm::runAlgorithm() {
    buildMainESA();
    buildSingleESAs();
    establishIdentifierMappingsViaLCA_OptimalMerged();
    performFrequencySelectionDFS();
}

vector<string> readStringsFromFile(const string& filename) {
    vector<string> strings;
    ifstream file(filename);

    if (!file.is_open()) {
        throw runtime_error("Cannot open input file: " + filename);
    }

    string line;
    while (getline(file, line)) {
        if (!line.empty()) {
            strings.push_back(line);
        }
    }

    file.close();

    if (strings.empty()) {
        throw runtime_error("No valid strings found in input file");
    }

    return strings;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    try {
        vector<string> inputStrings = readStringsFromFile(inputFile);
        
        MultiColorTreeAlgorithm algorithm(inputStrings);
        algorithm.runAlgorithm();
        
        ofstream file(outputFile);
        if (file.is_open()) {
            file << "Algorithm completed successfully." << endl;
            file.close();
        }
        
        return 0;
        
    } catch (const exception& e) {
        return 1;
    }
}