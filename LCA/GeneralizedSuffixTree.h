#ifndef GENERALIZEDSUFFIXTREE_H_
#define GENERALIZEDSUFFIXTREE_H_

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <functional>
#include <unordered_set>

#include "ActivePoint.h"
#include "Node.h"

using namespace std;

class GeneralizedSuffixTree {

public:
    GeneralizedSuffixTree(string s);
    GeneralizedSuffixTree(string s1, string s2);
    GeneralizedSuffixTree(vector<string>& strings, GeneralizedSuffixTree* leftSubtree, GeneralizedSuffixTree* rightSubtree);
    virtual ~GeneralizedSuffixTree();

    int addString(string& s);
    void exportInDotFormat(char* fileName);

    char getLastChar(Node& n);
    string getEdgeString(Node& node);
    string getEdgeStringWithTerm(Node& node, bool withTerm);

    int lastInsertedNode = -1;
    int currentStringIdx = -1;
    short currentCharIdx = -1;

    vector<string> strings;
    int rootIdx;
    unordered_map<int, Node> nodes;

    void updateFrequenciesDFS_Lower();
    pair<int, int> updateNodeFrequencies_Lower(int nodeIdx);
    
    void updateFrequenciesDFS_Higher(GeneralizedSuffixTree* leftSubtree = nullptr, GeneralizedSuffixTree* rightSubtree = nullptr);
    void generateELR(vector<int>& E, vector<int>& L, vector<int>& R);
    
    void batchLCA(const vector<pair<int, int> >& nodePairs, vector<int>& results);
    
    void traverseDFS();

private:
    void initialize(vector<string>& strings);
    void clearStringInfo();

    void addNextString();
    void extend();

    char getActiveEdge();
    Node getActiveNode();

    bool walkDown(int node);

    int addNode(int stringIdx, short labelStartIdx, short labelEndInx);
    int addLeaf(int stringIdx, short labelStartIdx, short labelEndInx);

    void setSuffixLink(int node);

    void generateELRDFS(int nodeIdx, int depth, vector<int>& E, vector<int>& L, vector<int>& R);

    ActivePoint activePoint;
    int currentNodeID = -1;
    short currentStringLength = -1;

    pair<vector<Node*>, vector<Node*> > collectLeafNodeInfo(
        int nodeIdx, int leftGroupSize,
        map<int, Node*>& leftGroupLeftmostMap,
        map<int, Node*>& leftGroupRightmostMap,
        map<int, Node*>& rightGroupLeftmostMap,
        map<int, Node*>& rightGroupRightmostMap);
};

#endif /* GENERALIZEDSUFFIXTREE_H_ */