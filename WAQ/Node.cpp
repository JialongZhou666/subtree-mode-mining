#include "Node.h"

using namespace std;

Node::Node() {
    max_freq = 0;
    min_freq = 0;
    string_left = 0;
    string_right = 0;
    depth = 0;
}

Node::Node(int si, short lsi, short lei) {
    stringIdx = si;
    labelStartIdx = lsi;
    labelEndIdx = lei;
    max_freq = 0;
    min_freq = 0;
    string_left = 0;
    string_right = 0;
    depth = 0;
}

Node::~Node() {
}

void Node::addSuffix(int stringIdx, short charIdx){
    suffixes.push_back(make_pair(stringIdx, charIdx));
}

short Node::getLabelLength(){
    return labelEndIdx - labelStartIdx + 1;
}