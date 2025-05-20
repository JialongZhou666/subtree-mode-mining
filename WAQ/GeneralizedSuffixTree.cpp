#include "GeneralizedSuffixTree.h"

#include <queue>
#include <fstream>
#include <iostream>
#include <cstring>
#include <set>
#include <algorithm>
#include <unordered_set>

using namespace std;

GeneralizedSuffixTree::GeneralizedSuffixTree(string s) {
    strings.push_back(s);
    initialize(strings);
    clearStringInfo();
    addNextString();
    
    calculateNodeDepths();
}

GeneralizedSuffixTree::GeneralizedSuffixTree(string s1, string s2) {
    vector<string> inputStrings;
    inputStrings.push_back(s1);
    inputStrings.push_back(s2);

    initialize(inputStrings);
    
    for (int i = 0; i < (int)inputStrings.size(); i++) {
        clearStringInfo();
        addNextString();
    }
    cout << "Added 2 strings [OK]" << endl;

    updateNodeInfoDFS_Lower(rootIdx, 0);
    cout << "Node info updated [OK]" << endl;
}

GeneralizedSuffixTree::GeneralizedSuffixTree(vector<string>& strings, GeneralizedSuffixTree* leftSubtree, GeneralizedSuffixTree* rightSubtree) {
    initialize(strings);

    for (int i = 0; i < (int) strings.size(); i++) {
        if (i % 1000 == 0)
            cout << "Adding string " << currentStringIdx << "\r";
        clearStringInfo();
        addNextString();
    }

    cout << "Adding string " << currentStringIdx - 1 << " [OK]" << endl;
    
    calculateNodeDepths();

    updateFrequenciesDFS_Higher(leftSubtree, rightSubtree);
}

GeneralizedSuffixTree::~GeneralizedSuffixTree() {
    nodes.clear();
    strings.clear();
}

int GeneralizedSuffixTree::addString(string& s) {
    strings.push_back(s);

    clearStringInfo();
    addNextString();

    calculateNodeDepths();

    return currentStringIdx - 1;
}

char GeneralizedSuffixTree::getLastChar(Node& n) {
    return strings[n.stringIdx][n.labelEndIdx];
}

string GeneralizedSuffixTree::getEdgeString(Node& node) {
    return getEdgeStringWithTerm(node, true);
}

string GeneralizedSuffixTree::getEdgeStringWithTerm(Node& node, bool withTerm) {
    if (node.labelStartIdx == node.labelEndIdx)
        return strings[node.stringIdx].substr(node.labelStartIdx, (withTerm ? 1 : 0));
    else
        return strings[node.stringIdx].substr(node.labelStartIdx, (node.labelEndIdx - node.labelStartIdx) + (withTerm ? 1 : 0));
}

void GeneralizedSuffixTree::initialize(vector<string>& ss) {
    strings = ss;

    currentNodeID = 0;
    rootIdx = addNode(0, 0, 0);

    currentStringIdx = 0;
}

void GeneralizedSuffixTree::clearStringInfo() {
    activePoint.reset();
    currentStringLength = 0;
}

void GeneralizedSuffixTree::addNextString() {
    strings[currentStringIdx].push_back('$');
    currentStringLength = (short) strings[currentStringIdx].length();

    for (currentCharIdx = 0; currentCharIdx < currentStringLength; currentCharIdx++)
        extend();

    currentStringIdx++;
}

void GeneralizedSuffixTree::extend() {

    char c = strings[currentStringIdx][currentCharIdx];

    lastInsertedNode = 0;
    activePoint.remainder++;

    while (activePoint.remainder > 0) {

        if (activePoint.idx == 0)
            activePoint.edge = currentCharIdx;

        if (nodes[activePoint.node].children.count(getActiveEdge()) == 0) {
            int newLeaf = addLeaf(currentStringIdx, currentCharIdx, currentStringLength - 1);
            nodes[activePoint.node].children[getActiveEdge()] = newLeaf;
            setSuffixLink(activePoint.node);
        } else {

            if (activePoint.node == rootIdx) {
                activePoint.edge = currentCharIdx - activePoint.remainder + 1;
                activePoint.idx = activePoint.remainder - 1;
            }

            int childID = nodes[activePoint.node].children[getActiveEdge()];
            Node& child = nodes[childID];

            if (walkDown(childID))
                continue;

            if (strings[child.stringIdx][child.labelStartIdx + activePoint.idx] == c) {
                if (c != '$') {
                    activePoint.idx++;
                    setSuffixLink(activePoint.node);
                    break;
                } else
                    child.addSuffix(currentStringIdx, currentCharIdx - (activePoint.remainder - 1));

            } else {
                int newInnerNode = addNode(child.stringIdx, child.labelStartIdx, child.labelStartIdx + activePoint.idx - 1);
                nodes[activePoint.node].children[getActiveEdge()] = newInnerNode;
                int newLeaf = addLeaf(currentStringIdx, currentCharIdx, currentStringLength - 1);
                nodes[newInnerNode].children[c] = newLeaf;
                child.labelStartIdx += activePoint.idx;
                nodes[newInnerNode].children[strings[child.stringIdx][child.labelStartIdx]] = childID;
                setSuffixLink(newInnerNode);
            }
        }

        activePoint.remainder--;

        if (activePoint.node == rootIdx && activePoint.idx > 0) {
            activePoint.idx--;
            activePoint.edge = currentCharIdx - activePoint.remainder + 1;
        } else {
            activePoint.node = nodes[activePoint.node].suffixLink > 0 ? nodes[activePoint.node].suffixLink : rootIdx;

        }
    }
}

char GeneralizedSuffixTree::getActiveEdge() {
    return strings[currentStringIdx][activePoint.edge];
}

Node GeneralizedSuffixTree::getActiveNode() {
    return nodes[activePoint.node];
}

bool GeneralizedSuffixTree::walkDown(int node) {
    short labelLength = nodes[node].getLabelLength();

    if (activePoint.idx >= labelLength) {
        activePoint.edge += labelLength;
        activePoint.idx -= labelLength;
        activePoint.node = node;
        return true;
    }
    return false;
}

int GeneralizedSuffixTree::addLeaf(int stringIdx, short labelStartIdx, short labelEndIdx) {
    nodes[currentNodeID] = *(new Node(stringIdx, labelStartIdx, labelEndIdx));
    nodes[currentNodeID].addSuffix(stringIdx, currentCharIdx - (activePoint.remainder - 1));
    return currentNodeID++;
}

int GeneralizedSuffixTree::addNode(int stringIdx, short labelStartIdx, short labelEndIdx) {
    nodes[currentNodeID] = *(new Node(stringIdx, labelStartIdx, labelEndIdx));
    return currentNodeID++;
}

void GeneralizedSuffixTree::setSuffixLink(int node) {
    if (lastInsertedNode > 0)
        nodes[lastInsertedNode].suffixLink = node;

    lastInsertedNode = node;
}

pair<int, int> GeneralizedSuffixTree::updateNodeInfoDFS_Lower(int nodeIdx, int parentDepth) {
    Node& node = nodes[nodeIdx];
    
    node.depth = parentDepth;
    
    if (node.children.empty()) {
        int string_left = (node.stringIdx == 0) ? 1 : 0;
        int string_right = (node.stringIdx == 0) ? 0 : 1;
        
        node.max_freq = max(string_left, string_right);
        node.min_freq = min(string_left, string_right);
        
        return make_pair(string_left, string_right);
    }
    
    int string_left = 0;
    int string_right = 0;
    
    for (auto& child_pair : node.children) {
        int childIdx = child_pair.second;
        Node& child = nodes[childIdx];
        
        int edgeLength = child.labelEndIdx - child.labelStartIdx + 1;
        
        pair<int, int> child_freqs = updateNodeInfoDFS_Lower(childIdx, parentDepth + edgeLength);
        
        string_left += child_freqs.first;
        string_right += child_freqs.second;
    }
    
    node.max_freq = max(string_left, string_right);
    node.min_freq = min(string_left, string_right);
    
    return make_pair(string_left, string_right);
}

void GeneralizedSuffixTree::updateFrequenciesDFS_Higher(GeneralizedSuffixTree* leftSubtree, GeneralizedSuffixTree* rightSubtree) {
    vector<pair<int, int>> nodeQueries;
    
    for (const auto& nodePair : nodes) {
        int nodeId = nodePair.first;
        const Node& node = nodePair.second;
        nodeQueries.push_back(make_pair(nodeId, node.depth));
    }

    vector<int> leftResults(nodeQueries.size(), -1);
    vector<int> rightResults(nodeQueries.size(), -1);
    
    leftSubtree->processWeightedAncestorQueries(nodeQueries, leftResults);
    
    rightSubtree->processWeightedAncestorQueries(nodeQueries, rightResults);
    
    for (size_t i = 0; i < nodeQueries.size(); i++) {
        int nodeId = nodeQueries[i].first;
        Node& currentNode = nodes[nodeId];
        
        int leftMinFreq = 0;
        int leftMaxFreq = 0;
        int rightMinFreq = 0;
        int rightMaxFreq = 0;
        
        int leftNodeId = leftResults[i];
        if (leftNodeId != -1 && leftSubtree->nodes.count(leftNodeId) > 0) {
            const Node& leftNode = leftSubtree->nodes[leftNodeId];
            leftMinFreq = leftNode.min_freq;
            leftMaxFreq = leftNode.max_freq;
        }
        
        int rightNodeId = rightResults[i];
        if (rightNodeId != -1 && rightSubtree->nodes.count(rightNodeId) > 0) {
            const Node& rightNode = rightSubtree->nodes[rightNodeId];
            rightMinFreq = rightNode.min_freq;
            rightMaxFreq = rightNode.max_freq;
        }
        
        currentNode.max_freq = max(leftMaxFreq, rightMaxFreq);
        
        if (leftNodeId != -1 && rightNodeId != -1 && leftMinFreq > 0 && rightMinFreq > 0) {
            currentNode.min_freq = min(leftMinFreq, rightMinFreq);
        } else {
            currentNode.min_freq = 0;
        }
    }
    
    cout << "Node frequencies updated from subtrees [OK]" << endl;
}

void GeneralizedSuffixTree::calculateNodeDepths() {
    unordered_set<int> visited;
    
    function<void(int, int)> dfs = [&](int nodeIdx, int parentDepth) {
        if (nodes.find(nodeIdx) == nodes.end()) {
            cerr << "Warning: Node " << nodeIdx << " not found!" << endl;
            return;
        }
        
        if (visited.find(nodeIdx) != visited.end()) {
            cerr << "Warning: Cycle detected at node " << nodeIdx << endl;
            return;
        }
        
        visited.insert(nodeIdx);
        
        Node& node = nodes[nodeIdx];
        
        node.depth = parentDepth;
        
        for (auto& child_pair : node.children) {
            int childIdx = child_pair.second;
            
            if (nodes.find(childIdx) == nodes.end()) {
                cerr << "Warning: Child node " << childIdx << " not found!" << endl;
                continue;
            }
            
            Node& child = nodes[childIdx];
            
            int edgeLength = child.labelEndIdx - child.labelStartIdx + 1;
            
            dfs(childIdx, parentDepth + edgeLength);
        }
    };
    
    dfs(rootIdx, 0);
    
    cout << "Node depths calculated [OK]" << endl;
}

void GeneralizedSuffixTree::processWeightedAncestorQueries(
    const vector<pair<int, int>>& queries, vector<int>& results) {
    
    set<int> uniqueDepths;
    for (const auto& nodePair : nodes) {
        uniqueDepths.insert(nodePair.second.depth);
    }
    
    for (const auto& query : queries) {
        uniqueDepths.insert(query.second);
    }
    
    vector<int> depths(uniqueDepths.begin(), uniqueDepths.end());
    sort(depths.begin(), depths.end(), greater<int>());
    
    unordered_map<int, int> uf;
    for (const auto& nodePair : nodes) {
        uf[nodePair.first] = -1;
    }
    
    auto find = [&uf](int x) -> int {
        if (uf.find(x) == uf.end()) {
            return -1;
        }
        
        int r = x;
        while (uf[r] >= 0) r = uf[r];
        
        while (x != r && uf[x] >= 0) {
            int tmp = uf[x];
            uf[x] = r;
            x = tmp;
        }
        
        return r;
    };
    
    auto join = [&uf, &find](int x, int y) {
        x = find(x);
        y = find(y);
        
        if (x == -1 || y == -1 || x == y) return;
        
        if (uf[x] < uf[y]) {
            uf[x] += uf[y];
            uf[y] = x;
        } else {
            uf[y] += uf[x];
            uf[x] = y;
        }
    };
    
    for (int depth : depths) {
        for (const auto& nodePair : nodes) {
            int nodeIdx = nodePair.first;
            const Node& node = nodePair.second;
            
            if (node.depth == depth) {
                for (const auto& child_pair : node.children) {
                    join(nodeIdx, child_pair.second);
                }
            }
        }
        
        for (size_t i = 0; i < queries.size(); i++) {
            if (queries[i].second == depth) {
                results[i] = find(queries[i].first);
            }
        }
    }
    
    cout << "Weighted ancestor queries processed [OK]" << endl;
}

void GeneralizedSuffixTree::traverseDFS() {
    unordered_set<int> visited;
    
    int totalNodes = 0;
    int fairSubstrings = 0;
    
    function<void(int, int, string)> dfs = [&](int nodeIdx, int depth, string currentString) {
        if (visited.find(nodeIdx) != visited.end()) {
            return;
        }
        
        visited.insert(nodeIdx);
        totalNodes++;
        
        Node& node = nodes[nodeIdx];
        
        string edgeString;
        if (nodeIdx != rootIdx) {
            if (node.labelStartIdx <= node.labelEndIdx && 
                node.stringIdx < strings.size()) {
                edgeString = strings[node.stringIdx].substr(
                    node.labelStartIdx, 
                    node.labelEndIdx - node.labelStartIdx + 1);
                
                if (!edgeString.empty() && edgeString.back() == '$') {
                    edgeString.pop_back();
                }
                
                currentString += edgeString;
            }
        }
        
        if (node.min_freq > 0 && (node.max_freq - node.min_freq < 1)) {
            fairSubstrings++;
            if (!currentString.empty()) {
                cout << "Fair substring: \"" << currentString << "\" (min_freq: " << node.min_freq 
                     << ", max_freq: " << node.max_freq << ")" << endl;
            }
        }
        
        for (auto& child_pair : node.children) {
            dfs(child_pair.second, depth + 1, currentString);
        }
    };
    
    dfs(rootIdx, 0, "");
}