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

#ifndef MULTI_COLOR_ESA_H
#define MULTI_COLOR_ESA_H

#include "defs.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>
#include <sdsl/rmq_support.hpp>

using namespace std;

struct SCMStatistics {
    INT totalMainTreeNodes = 0;
    INT totalSingleColorNodes = 0;
    INT numColors = 0;
};

extern SCMStatistics globalSCMStats;

struct MainTreeNode {
    INT lcpDepth;
    INT leftBound;
    INT rightBound;
    INT parent;
    vector<INT> children;
    
    MainTreeNode() : lcpDepth(0), leftBound(0), rightBound(0), parent(-1) {}
};

struct SingleColorNode {
    INT lcpDepth;
    INT leftBound;
    INT rightBound;
    INT parent;
    vector<INT> children;
    INT count;
    INT mainTreeId;
    
    SingleColorNode() : lcpDepth(0), leftBound(0), rightBound(0), 
                        parent(-1), count(0), mainTreeId(-1) {}
};

class MultiColorESA {
private:
    string T;
    INT n;
    INT k;
    
    vector<INT> SA;
    vector<INT> LCP;
    vector<INT> colorBoundaries;
    
    vector<MainTreeNode> mainTree;
    vector<vector<SingleColorNode>> singleColorTrees;
    
    sdsl::rmq_succinct_sct<> rmq;
    
    struct TupleHash {
        size_t operator()(const tuple<INT, INT, INT>& t) const {
            size_t h1 = hash<INT>{}(get<0>(t));
            size_t h2 = hash<INT>{}(get<1>(t));
            size_t h3 = hash<INT>{}(get<2>(t));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
    
    unordered_map<tuple<INT, INT, INT>, INT, TupleHash> mainTreeIndex;
    
    vector<vector<INT>> precomputedColorLeaves;
    vector<vector<INT>> precomputedAdjacentLCPs;
    
    void buildLCP();
    void buildMainTree();
    void buildMainTreeIndex();
    void precomputeColorData();
    
    void buildSingleColorTreeKasai(INT color, 
                                    const vector<INT>& colorLeaves,
                                    const vector<INT>& adjacentLCPs);
    
    INT getColor(INT saIndex);
    INT getLCPBetween(INT i, INT j);
    INT findMainTreeNode(INT leftSA, INT rightSA, INT targetLCP);
    
    void splitTree();
    void countColors();
    void mergeCounts();
    
public:
    MultiColorESA(const vector<string>& sequences);
    void runSCM();
};

#endif // MULTI_COLOR_ESA_H