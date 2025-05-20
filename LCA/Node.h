#ifndef NODE_H_
#define NODE_H_

#include <unordered_map>
#include <vector>

using namespace std;

class Node {
public:

	typedef pair<int,short> suffix;

	Node();
	Node(int stringIdx, short labelStartIdx, short labelEndInx);
	virtual ~Node();

	void addSuffix(int stringIdx, short suffixIdx);
	short getLabelLength();

	int stringIdx = -1;
	short labelStartIdx = -1;
	short labelEndIdx = -1;

	int suffixLink = -1;
	unordered_map<char,int> children;

	vector<suffix> suffixes;

	int max_freq = 0;
	int min_freq = 0;
	int string_left = 0;
	int string_right = 0;
};

#endif /* NODE_H_ */
