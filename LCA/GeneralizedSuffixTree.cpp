#include "GeneralizedSuffixTree.h"
#include "rmq-offline.h"

#include <queue>
#include <fstream>
#include <iostream>
#include <cstring>
#include <unordered_set>

using namespace std;

GeneralizedSuffixTree::GeneralizedSuffixTree(string s) {
    strings.push_back(s);
    initialize(strings);
    clearStringInfo();
    addNextString();
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
    
    updateFrequenciesDFS_Lower();
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

    updateFrequenciesDFS_Higher(leftSubtree, rightSubtree);
}

GeneralizedSuffixTree::~GeneralizedSuffixTree() {
}

int GeneralizedSuffixTree::addString(string& s) {
    strings.push_back(s);
    clearStringInfo();
    addNextString();
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

void GeneralizedSuffixTree::updateFrequenciesDFS_Lower() {
    updateNodeFrequencies_Lower(rootIdx);
    cout << "Frequencies updated [OK]" << endl;
}

pair<int, int> GeneralizedSuffixTree::updateNodeFrequencies_Lower(int nodeIdx) {
    Node& node = nodes[nodeIdx];
    
    if (node.children.empty()) {
        node.string_left = (node.stringIdx == 0) ? 1 : 0;
        node.string_right = (node.stringIdx == 0) ? 0 : 1;
        
        node.max_freq = max(node.string_left, node.string_right);
        node.min_freq = min(node.string_left, node.string_right);
        
        return make_pair(node.string_left, node.string_right);
    }
    
    node.string_left = 0;
    node.string_right = 0;
    
    for (auto& child_pair : node.children) {
        pair<int, int> child_freqs = updateNodeFrequencies_Lower(child_pair.second);
        node.string_left += child_freqs.first;
        node.string_right += child_freqs.second;
    }
    
    node.max_freq = max(node.string_left, node.string_right);
    node.min_freq = min(node.string_left, node.string_right);
    
    return make_pair(node.string_left, node.string_right);
}

void GeneralizedSuffixTree::updateFrequenciesDFS_Higher(GeneralizedSuffixTree* leftSubtree, GeneralizedSuffixTree* rightSubtree) {
    int totalStrings = (int)strings.size();
    int leftGroupSize = totalStrings / 2;
    
    map<int, Node*> leftGroupLeftmostMap;
    map<int, Node*> leftGroupRightmostMap;
    map<int, Node*> rightGroupLeftmostMap;
    map<int, Node*> rightGroupRightmostMap;
    
    collectLeafNodeInfo(rootIdx, leftGroupSize, leftGroupLeftmostMap, leftGroupRightmostMap, 
                        rightGroupLeftmostMap, rightGroupRightmostMap);
    
    vector<pair<int, int>> leftNodePairs;
    vector<int> leftLCAResults;

    for (auto& entry : leftGroupLeftmostMap) {
        int nodeIdx = entry.first;
        leftNodePairs.push_back(make_pair(nodeIdx, nodeIdx));
    }

    if (!leftNodePairs.empty()) {
        leftLCAResults.resize(leftNodePairs.size());
        batchLCA(leftNodePairs, leftLCAResults);
    }

    vector<pair<int, int>> rightNodePairs;
    vector<int> rightLCAResults;

    for (auto& entry : rightGroupLeftmostMap) {
        int nodeIdx = entry.first;
        rightNodePairs.push_back(make_pair(nodeIdx, nodeIdx));
    }

    if (!rightNodePairs.empty()) {
        rightLCAResults.resize(rightNodePairs.size());
        batchLCA(rightNodePairs, rightLCAResults);
    }
    
    for (auto& nodePair : nodes) {
        Node& node = nodePair.second;
        node.string_left = 0;
        node.string_right = 0;
    }
    
    int leftIndex = 0;
    for (auto& entry : leftGroupLeftmostMap) {
        int nodeIdx = entry.first;
        int lcaResult = leftLCAResults[leftIndex++];
        
        if (lcaResult != -1 && leftSubtree && leftSubtree->nodes.count(lcaResult) > 0) {
            const Node& leftNode = leftSubtree->nodes[lcaResult];
            nodes[nodeIdx].string_left = leftNode.string_left + leftNode.string_right;
        }
    }
    
    int rightIndex = 0;
    for (auto& entry : rightGroupLeftmostMap) {
        int nodeIdx = entry.first;
        int lcaResult = rightLCAResults[rightIndex++];
        
        if (lcaResult != -1 && rightSubtree && rightSubtree->nodes.count(lcaResult) > 0) {
            const Node& rightNode = rightSubtree->nodes[lcaResult];
            nodes[nodeIdx].string_right = rightNode.string_left + rightNode.string_right;
        }
    }
    
    for (auto& nodePair : nodes) {
        Node& node = nodePair.second;
        node.max_freq = max(node.string_left, node.string_right);
        node.min_freq = min(node.string_left, node.string_right);
    }
    
    cout << "Node frequencies updated based on LCA information [OK]" << endl;
}

pair<vector<Node*>, vector<Node*> > GeneralizedSuffixTree::collectLeafNodeInfo(
    int nodeIdx, int leftGroupSize,
    map<int, Node*>& leftGroupLeftmostMap,
    map<int, Node*>& leftGroupRightmostMap,
    map<int, Node*>& rightGroupLeftmostMap,
    map<int, Node*>& rightGroupRightmostMap) {
    
    Node& node = nodes[nodeIdx];
    
    vector<Node*> leftGroupInfo;
    vector<Node*> rightGroupInfo;
    
    if (node.children.empty()) {
        bool isLeftGroup = (node.stringIdx < leftGroupSize);
        
        if (isLeftGroup) {
            leftGroupInfo.push_back(&node);
            leftGroupInfo.push_back(&node);
            
            leftGroupLeftmostMap[nodeIdx] = &node;
            leftGroupRightmostMap[nodeIdx] = &node;
        } else {
            rightGroupInfo.push_back(&node);
            rightGroupInfo.push_back(&node);
            
            rightGroupLeftmostMap[nodeIdx] = &node;
            rightGroupRightmostMap[nodeIdx] = &node;
        }
        
        return make_pair(leftGroupInfo, rightGroupInfo);
    }
    
    Node* leftGroupLeftmost = nullptr;
    Node* leftGroupRightmost = nullptr;
    Node* rightGroupLeftmost = nullptr;
    Node* rightGroupRightmost = nullptr;
    
    for (auto& child_pair : node.children) {
        pair<vector<Node*>, vector<Node*> > childInfo = collectLeafNodeInfo(
            child_pair.second, leftGroupSize, 
            leftGroupLeftmostMap, leftGroupRightmostMap, 
            rightGroupLeftmostMap, rightGroupRightmostMap);
        
        vector<Node*>& childLeftGroupInfo = childInfo.first;
        vector<Node*>& childRightGroupInfo = childInfo.second;
        
        if (!childLeftGroupInfo.empty()) {
            Node* childLeftmost = childLeftGroupInfo[0];
            Node* childRightmost = childLeftGroupInfo[1];
            
            if (!leftGroupLeftmost || 
                childLeftmost->stringIdx < leftGroupLeftmost->stringIdx ||
                (childLeftmost->stringIdx == leftGroupLeftmost->stringIdx && 
                 childLeftmost->labelStartIdx < leftGroupLeftmost->labelStartIdx)) {
                leftGroupLeftmost = childLeftmost;
            }
            
            if (!leftGroupRightmost || 
                childRightmost->stringIdx > leftGroupRightmost->stringIdx ||
                (childRightmost->stringIdx == leftGroupRightmost->stringIdx && 
                 childRightmost->labelStartIdx > leftGroupRightmost->labelStartIdx)) {
                leftGroupRightmost = childRightmost;
            }
        }
        
        if (!childRightGroupInfo.empty()) {
            Node* childLeftmost = childRightGroupInfo[0];
            Node* childRightmost = childRightGroupInfo[1];
            
            if (!rightGroupLeftmost || 
                childLeftmost->stringIdx < rightGroupLeftmost->stringIdx ||
                (childLeftmost->stringIdx == rightGroupLeftmost->stringIdx && 
                 childLeftmost->labelStartIdx < rightGroupLeftmost->labelStartIdx)) {
                rightGroupLeftmost = childLeftmost;
            }
            
            if (!rightGroupRightmost || 
                childRightmost->stringIdx > rightGroupRightmost->stringIdx ||
                (childRightmost->stringIdx == rightGroupRightmost->stringIdx && 
                 childRightmost->labelStartIdx > rightGroupRightmost->labelStartIdx)) {
                rightGroupRightmost = childRightmost;
            }
        }
    }
    
    if (leftGroupLeftmost && leftGroupRightmost) {
        leftGroupLeftmostMap[nodeIdx] = leftGroupLeftmost;
        leftGroupRightmostMap[nodeIdx] = leftGroupRightmost;
        
        leftGroupInfo.push_back(leftGroupLeftmost);
        leftGroupInfo.push_back(leftGroupRightmost);
    }
    
    if (rightGroupLeftmost && rightGroupRightmost) {
        rightGroupLeftmostMap[nodeIdx] = rightGroupLeftmost;
        rightGroupRightmostMap[nodeIdx] = rightGroupRightmost;
        
        rightGroupInfo.push_back(rightGroupLeftmost);
        rightGroupInfo.push_back(rightGroupRightmost);
    }
    
    return make_pair(leftGroupInfo, rightGroupInfo);
}

void GeneralizedSuffixTree::generateELR(vector<int>& E, vector<int>& L, vector<int>& R) {
    int n = nodes.size();
    
    E.clear();
    L.clear();
    R.assign(n, -1);
    
    E.reserve(2*n-1);
    L.reserve(2*n-1);
    
    generateELRDFS(rootIdx, 0, E, L, R);
    
    cout << "Generated Euler tour with " << E.size() << " entries." << endl;
}

void GeneralizedSuffixTree::generateELRDFS(int nodeIdx, int depth, vector<int>& E, vector<int>& L, vector<int>& R) {
    if (R[nodeIdx] == -1) {
        R[nodeIdx] = E.size();
    }
    
    E.push_back(nodeIdx);
    L.push_back(depth);
    
    vector<pair<char, int> > sortedChildren;
    for (const auto& child : nodes[nodeIdx].children) {
        sortedChildren.push_back(child);
    }
    sort(sortedChildren.begin(), sortedChildren.end());
    
    for (size_t i = 0; i < sortedChildren.size(); i++) {
        int childIdx = sortedChildren[i].second;
        
        generateELRDFS(childIdx, depth + 1, E, L, R);
        
        if (i < sortedChildren.size() - 1) {
            E.push_back(nodeIdx);
            L.push_back(depth);
        }
    }
}

void GeneralizedSuffixTree::batchLCA(const vector<pair<int, int> >& nodePairs, 
                                     vector<int>& results) {
    vector<int> E, L, R;
    generateELR(E, L, R);
    
    int q = nodePairs.size();
    results.resize(q);
    
    Query* queries = (Query*)calloc(q, sizeof(Query));
    
    for (int i = 0; i < q; i++) {
        int u = nodePairs[i].first;
        int v = nodePairs[i].second;
        
        int pos_u = R[u];
        int pos_v = R[v];
        
        if (pos_u > pos_v) {
            swap(pos_u, pos_v);
        }
        
        queries[i].L = pos_u;
        queries[i].R = pos_v;
        queries[i].O = -1;
    }
    
    INT* L_array = (INT*)calloc(L.size(), sizeof(INT));
    for (size_t i = 0; i < L.size(); i++) {
        L_array[i] = L[i];
    }
    
    rmq_offline(L_array, L.size(), queries, q);
    
    for (int i = 0; i < q; i++) {
        int pos_lca = queries[i].O;
        results[i] = E[pos_lca];
    }
    
    free(queries);
    free(L_array);
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