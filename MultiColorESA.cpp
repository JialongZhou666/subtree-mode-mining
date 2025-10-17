/**
    SCM
    Copyright (C) 2025 Jialong Zhou.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include "MultiColorESA.h"
#include <divsufsort64.h>
#include <cstring>
#include <climits>
#include <queue>
#include <set>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <tuple>
#include <functional>

using namespace std;
using namespace chrono;

const char SEPARATOR = '$';

SCMStatistics globalSCMStats;

MultiColorESA::MultiColorESA(const vector<string>& sequences)
    : k(sequences.size()) {
    
    colorBoundaries.push_back(0);
    for (INT i = 0; i < (INT)sequences.size(); i++) {
        T += sequences[i];
        T += SEPARATOR;
        colorBoundaries.push_back(T.length());
    }
    n = T.length();
    
    SA.resize(n);
    sauchar_t* text_ptr = (sauchar_t*)T.c_str();
    divsufsort64(text_ptr, (saidx64_t*)SA.data(), n);
    
    buildLCP();
    buildMainTree();
    
    sdsl::util::assign(rmq, sdsl::rmq_succinct_sct<>(&LCP));
    
    buildMainTreeIndex();
    precomputeColorData();
    
    globalSCMStats.totalMainTreeNodes = mainTree.size();
}

void MultiColorESA::buildLCP() {
    LCP.resize(n, 0);
    vector<INT> rank(n);
    
    for (INT i = 0; i < n; i++) {
        rank[SA[i]] = i;
    }
    
    INT h = 0;
    for (INT i = 0; i < n; i++) {
        if (rank[i] > 0) {
            INT j = SA[rank[i] - 1];
            while (i + h < n && j + h < n && T[i + h] == T[j + h]) {
                h++;
            }
            LCP[rank[i]] = h;
            if (h > 0) h--;
        }
    }
}

void MultiColorESA::buildMainTree() {
    mainTree.clear();
    mainTree.reserve(2 * n);
    
    struct StackEntry {
        INT lcpValue;
        INT nodeIndex;
    };
    
    vector<StackEntry> stack;
    stack.reserve(n);
    
    for (INT i = 0; i < n; i++) {
        INT currentLCP = (i < n - 1) ? LCP[i + 1] : 0;
        INT lastChildIndex = i;
        
        while (!stack.empty() && stack.back().lcpValue > currentLCP) {
            StackEntry top = stack.back();
            stack.pop_back();
            
            INT parentIdx = stack.empty() ? -1 : stack.back().nodeIndex;
            
            mainTree[top.nodeIndex].parent = parentIdx;
            if (parentIdx != -1) {
                mainTree[parentIdx].children.push_back(top.nodeIndex);
                mainTree[parentIdx].rightBound = mainTree[top.nodeIndex].rightBound;
            }
            
            lastChildIndex = top.nodeIndex;
        }
        
        if (stack.empty() || stack.back().lcpValue < currentLCP) {
            mainTree.push_back(MainTreeNode());
            INT newNodeIdx = mainTree.size() - 1;
            mainTree[newNodeIdx].lcpDepth = currentLCP;
            mainTree[newNodeIdx].leftBound = i;
            mainTree[newNodeIdx].rightBound = i;
            mainTree[newNodeIdx].children.push_back(lastChildIndex);
            
            stack.push_back({currentLCP, newNodeIdx});
        }
    }
    
    while (!stack.empty()) {
        StackEntry top = stack.back();
        stack.pop_back();
        
        INT parentIdx = stack.empty() ? -1 : stack.back().nodeIndex;
        mainTree[top.nodeIndex].parent = parentIdx;
        if (parentIdx != -1) {
            mainTree[parentIdx].children.push_back(top.nodeIndex);
            mainTree[parentIdx].rightBound = mainTree[top.nodeIndex].rightBound;
        }
    }
}

void MultiColorESA::buildMainTreeIndex() {
    mainTreeIndex.clear();
    mainTreeIndex.reserve(mainTree.size());
    
    for (INT i = 0; i < (INT)mainTree.size(); i++) {
        const MainTreeNode& node = mainTree[i];
        auto key = make_tuple(node.leftBound, node.rightBound, node.lcpDepth);
        mainTreeIndex[key] = i;
    }
}

void MultiColorESA::precomputeColorData() {
    precomputedColorLeaves.resize(k);
    
    for (INT i = 0; i < n; i++) {
        INT pos = SA[i];
        auto it = upper_bound(colorBoundaries.begin(), colorBoundaries.end(), pos);
        if (it != colorBoundaries.begin()) {
            INT color = (it - colorBoundaries.begin() - 1);
            if (color >= 0 && color < k) {
                precomputedColorLeaves[color].push_back(i);
            }
        }
    }
    
    precomputedAdjacentLCPs.resize(k);
    
    for (INT color = 0; color < k; color++) {
        const vector<INT>& leaves = precomputedColorLeaves[color];
        INT numLeaves = leaves.size();
        
        if (numLeaves <= 1) continue;
        
        precomputedAdjacentLCPs[color].reserve(numLeaves - 1);
        
        for (INT i = 0; i < numLeaves - 1; i++) {
            INT pos1 = leaves[i];
            INT pos2 = leaves[i + 1];
            if (pos1 > pos2) swap(pos1, pos2);
            
            INT lcp = LCP[rmq(pos1 + 1, pos2)];
            precomputedAdjacentLCPs[color].push_back(lcp);
        }
    }
}

INT MultiColorESA::getColor(INT saIndex) {
    if (saIndex < 0 || saIndex >= n) return -1;
    
    INT pos = SA[saIndex];
    
    auto it = upper_bound(colorBoundaries.begin(), colorBoundaries.end(), pos);
    if (it == colorBoundaries.begin()) return -1;
    
    return (it - colorBoundaries.begin() - 1);
}

INT MultiColorESA::getLCPBetween(INT i, INT j) {
    if (i == j) return n - SA[i];
    if (i > j) swap(i, j);
    return LCP[rmq(i + 1, j)];
}

INT MultiColorESA::findMainTreeNode(INT leftSA, INT rightSA, INT targetLCP) {
    auto key = make_tuple(leftSA, rightSA, targetLCP);
    auto it = mainTreeIndex.find(key);
    if (it != mainTreeIndex.end()) {
        return it->second;
    }
    return -1;
}

void MultiColorESA::buildSingleColorTreeKasai(INT color, const vector<INT>& colorLeaves,
                                               const vector<INT>& adjacentLCPs) {
    vector<SingleColorNode>& colorTree = singleColorTrees[color];
    colorTree.clear();
    
    if (colorLeaves.empty()) return;
    
    INT numLeaves = colorLeaves.size();
    colorTree.reserve(2 * numLeaves);
    
    for (INT leaf : colorLeaves) {
        colorTree.push_back(SingleColorNode());
        INT leafIdx = colorTree.size() - 1;
        colorTree[leafIdx].leftBound = leaf;
        colorTree[leafIdx].rightBound = leaf;
    }
    
    if (numLeaves == 1) return;
    
    struct StackEntry {
        INT lcpValue;
        INT nodeIndex;
    };
    
    vector<StackEntry> stack;
    stack.reserve(numLeaves);
    
    for (INT i = 0; i < numLeaves; i++) {
        INT currentLCP = (i < numLeaves - 1) ? adjacentLCPs[i] : 0;
        INT lastChildIndex = i;
        
        while (!stack.empty() && stack.back().lcpValue > currentLCP) {
            StackEntry top = stack.back();
            stack.pop_back();
            
            INT parentIdx = stack.empty() ? -1 : stack.back().nodeIndex;
            colorTree[top.nodeIndex].parent = parentIdx;
            
            if (parentIdx != -1) {
                colorTree[parentIdx].children.push_back(top.nodeIndex);
                colorTree[parentIdx].rightBound = colorTree[top.nodeIndex].rightBound;
            }
            
            INT leftLeaf = colorTree[top.nodeIndex].leftBound;
            INT rightLeaf = colorTree[top.nodeIndex].rightBound;
            INT targetLCP = colorTree[top.nodeIndex].lcpDepth;
            
            colorTree[top.nodeIndex].mainTreeId = findMainTreeNode(leftLeaf, rightLeaf, targetLCP);
            
            lastChildIndex = top.nodeIndex;
        }
        
        if (stack.empty() || stack.back().lcpValue < currentLCP) {
            colorTree.push_back(SingleColorNode());
            INT newIdx = colorTree.size() - 1;
            colorTree[newIdx].lcpDepth = currentLCP;
            colorTree[newIdx].leftBound = colorLeaves[i];
            colorTree[newIdx].rightBound = colorLeaves[i];
            colorTree[newIdx].children.push_back(lastChildIndex);
            
            stack.push_back({currentLCP, newIdx});
        }
    }
    
    while (!stack.empty()) {
        StackEntry top = stack.back();
        stack.pop_back();
        
        INT parentIdx = stack.empty() ? -1 : stack.back().nodeIndex;
        colorTree[top.nodeIndex].parent = parentIdx;
        
        if (parentIdx != -1) {
            colorTree[parentIdx].children.push_back(top.nodeIndex);
            colorTree[parentIdx].rightBound = colorTree[top.nodeIndex].rightBound;
        }
        
        INT leftLeaf = colorTree[top.nodeIndex].leftBound;
        INT rightLeaf = colorTree[top.nodeIndex].rightBound;
        INT targetLCP = colorTree[top.nodeIndex].lcpDepth;
        colorTree[top.nodeIndex].mainTreeId = findMainTreeNode(leftLeaf, rightLeaf, targetLCP);
    }
}

void MultiColorESA::splitTree() {
    singleColorTrees.resize(k);
    
    for (INT color = 0; color < k; color++) {
        buildSingleColorTreeKasai(color, 
                                  precomputedColorLeaves[color], 
                                  precomputedAdjacentLCPs[color]);
    }
    
    INT totalColorNodes = 0;
    for (INT i = 0; i < (INT)singleColorTrees.size(); i++) {
        totalColorNodes += singleColorTrees[i].size();
    }
    globalSCMStats.totalSingleColorNodes = totalColorNodes;
}

void MultiColorESA::countColors() {
    for (INT color = 0; color < (INT)singleColorTrees.size(); color++) {
        vector<SingleColorNode>& colorTree = singleColorTrees[color];
        
        if (colorTree.empty()) continue;
        
        function<INT(INT)> dfs = [&](INT nodeIdx) -> INT {
            SingleColorNode& node = colorTree[nodeIdx];
            
            if (node.children.empty()) {
                node.count = 1;
                return 1;
            }
            
            INT totalCount = 0;
            for (INT childIdx : node.children) {
                totalCount += dfs(childIdx);
            }
            
            node.count = totalCount;
            return totalCount;
        };
        
        for (INT i = 0; i < (INT)colorTree.size(); i++) {
            if (colorTree[i].parent == -1) {
                dfs(i);
            }
        }
    }
}

void MultiColorESA::mergeCounts() {
    vector<pair<INT, INT>> results(mainTree.size(), {0, -1});
    
    vector<INT> saColors(n, -1);
    for (INT i = 0; i < n; i++) {
        INT pos = SA[i];
        auto it = upper_bound(colorBoundaries.begin(), colorBoundaries.end(), pos);
        if (it != colorBoundaries.begin()) {
            saColors[i] = (it - colorBoundaries.begin() - 1);
        }
    }
    
    vector<vector<pair<INT, INT>>> mainToColorNodes(mainTree.size());
    
    for (INT color = 0; color < (INT)singleColorTrees.size(); color++) {
        const vector<SingleColorNode>& colorTree = singleColorTrees[color];
        
        for (INT i = 0; i < (INT)colorTree.size(); i++) {
            INT mainId = colorTree[i].mainTreeId;
            INT count = colorTree[i].count;
            
            if (mainId >= 0 && mainId < (INT)mainTree.size() && count > 0) {
                mainToColorNodes[mainId].push_back({color, count});
            }
        }
    }
    
    vector<INT> internalChildCount(mainTree.size(), 0);
    
    for (INT i = 0; i < (INT)mainTree.size(); i++) {
        for (INT childIdx : mainTree[i].children) {
            if (childIdx >= 0 && childIdx < (INT)mainTree.size()) {
                internalChildCount[i]++;
            }
        }
    }
    
    queue<INT> q;
    for (INT i = 0; i < (INT)mainTree.size(); i++) {
        if (internalChildCount[i] == 0) {
            q.push(i);
        }
    }
    
    vector<bool> processed(mainTree.size(), false);
    
    while (!q.empty()) {
        INT nodeIdx = q.front();
        q.pop();
        
        if (processed[nodeIdx]) continue;
        processed[nodeIdx] = true;
        
        if (nodeIdx < 0 || nodeIdx >= (INT)mainTree.size()) continue;
        
        MainTreeNode& node = mainTree[nodeIdx];
        
        INT maxCount = 0;
        INT bestColor = -1;
        
        for (INT childIdx : node.children) {
            pair<INT, INT> childResult;
            
            if (childIdx >= 0 && childIdx < (INT)mainTree.size()) {
                childResult = results[childIdx];
            } else if (childIdx >= (INT)mainTree.size() && childIdx < n) {
                INT saIdx = childIdx;
                INT color = saColors[saIdx];
                if (color >= 0 && color < k) {
                    childResult = {1, color};
                } else {
                    childResult = {0, -1};
                }
            } else {
                childResult = {0, -1};
            }
            
            if (childResult.first > maxCount) {
                maxCount = childResult.first;
                bestColor = childResult.second;
            }
        }
        
        for (const auto& p : mainToColorNodes[nodeIdx]) {
            INT color = p.first;
            INT count = p.second;
            
            if (count > maxCount) {
                maxCount = count;
                bestColor = color;
            }
        }
        
        results[nodeIdx] = {maxCount, bestColor};
        
        INT parentIdx = node.parent;
        if (parentIdx >= 0 && parentIdx < (INT)mainTree.size()) {
            internalChildCount[parentIdx]--;
            if (internalChildCount[parentIdx] == 0) {
                q.push(parentIdx);
            }
        }
    }
}

void MultiColorESA::runSCM() {
    globalSCMStats.numColors = k;
    
    splitTree();
    countColors();
    mergeCounts();
}