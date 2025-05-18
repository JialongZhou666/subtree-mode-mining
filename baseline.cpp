#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <functional>

using namespace std;

struct STvertex;
struct STedge {
    int l, r;
    STvertex* v;
    STedge(int _l = 0, int _r = 0, STvertex* _v = nullptr) : l(_l), r(_r), v(_v) {}
};

struct STvertex {
    map<unsigned char, STedge> g;
    STvertex* f;
    STvertex* parent;  
    char edge_char;    
    vector<int> freq;  
    vector<int> string_ids;
    int numer;
    
    STvertex() : f(nullptr), parent(nullptr), edge_char(0), numer(-1) {}
};

class SuffixTree {
private:
    STvertex* root;
    string combined_text;
    vector<int> str_starts;
    int string_count;
    int liscie;

    void Canonize(STedge& kraw) {
        if (kraw.l <= kraw.r) {
            STedge e = kraw.v->g[combined_text[kraw.l]];
            while (e.r - e.l <= kraw.r - kraw.l) {
                kraw.l += e.r - e.l + 1;
                kraw.v = e.v;
                if (kraw.l <= kraw.r) 
                    e = kraw.v->g[combined_text[kraw.l]];
            }
        }
    }

    bool Test_and_split(STvertex*& w, const STedge& kraw) {
        w = kraw.v;
        if (kraw.l <= kraw.r) {
            char c = combined_text[kraw.l];
            STedge e = kraw.v->g[c];
            if (combined_text[kraw.r + 1] == combined_text[e.l + kraw.r - kraw.l + 1]) 
                return true;
            
            w = new STvertex();
            w->numer = -1;
            kraw.v->g[c].r = e.l + kraw.r - kraw.l;
            kraw.v->g[c].v = w;
            w->parent = kraw.v;
            w->edge_char = c;
            e.l += kraw.r - kraw.l + 1;
            w->g[combined_text[e.l]] = e;
            e.v->parent = w;
            e.v->edge_char = combined_text[e.l];
            return false;
        }
        return kraw.v->g.find(combined_text[kraw.l]) != kraw.v->g.end();
    }

    void Update(STedge& kraw, int n) {
        STvertex* oldr = root;
        STvertex* w;
        
        while (!Test_and_split(w, kraw)) {
            STedge e;
            e.v = new STvertex();
            e.l = kraw.r + 1;
            e.r = n - 1;
            e.v->numer = liscie++;
            
            int pos = e.l;
            for (int i = str_starts.size() - 1; i >= 0; i--) {
                if (pos >= str_starts[i]) {
                    e.v->string_ids.push_back(i);
                    break;
                }
            }
            
            w->g[combined_text[kraw.r + 1]] = e;
            e.v->parent = w;
            e.v->edge_char = combined_text[kraw.r + 1];
            
            if (oldr != root) 
                oldr->f = w;
            oldr = w;
            kraw.v = kraw.v->f;
            Canonize(kraw);
        }
        if (oldr != root) 
            oldr->f = kraw.v;
    }

    string getEdgeString(const STedge& edge) {
        return combined_text.substr(edge.l, edge.r - edge.l + 1);
    }

    void getPath(STvertex* v, vector<pair<string,vector<int>>>& results) {
        if (v == root) return;
        
        string path = "";
        STvertex* curr = v;
        int lastEdgeLen = 0;
        
        while (curr != root) {
            STvertex* parent = curr->parent;
            for (auto& edge : parent->g) {
                if (edge.second.v == curr) {
                    string edgeStr = getEdgeString(edge.second);
                    if (curr == v) {
                        lastEdgeLen = edgeStr.length();
                    }
                    path = edgeStr + path;
                    break;
                }
            }
            curr = curr->parent;
        }
        
        for(int i = 0; i <= lastEdgeLen-1; i++) {
            results.push_back({path.substr(0, path.length()-i), v->freq});
        }
    }

    void dfsCalculateFreq(STvertex* v, int threshold, vector<pair<string, vector<int>>>& results) {
        if (!v) return;

        for (auto& edge : v->g) {
            dfsCalculateFreq(edge.second.v, threshold, results);
        }

        if (v != root) {
            if (v->g.empty()) {
            } else {
                fill(v->freq.begin(), v->freq.end(), 0);
                for (const auto& edge : v->g) {
                    for (int i = 0; i < string_count; i++) {
                        v->freq[i] += edge.second.v->freq[i];
                    }
                }
            }

            int maxFreq = *max_element(v->freq.begin(), v->freq.end());
            int minFreq = *min_element(v->freq.begin(), v->freq.end());
            
            if (minFreq > 0 && maxFreq - minFreq <= threshold) {
                vector<pair<string, vector<int>>> temp_results;
                getPath(v, temp_results);
                for(auto& res : temp_results) {
                    if (!res.first.empty() && res.first.find('$') == string::npos) {
                        results.push_back(res);
                    }
                }
            }
        }
    }

    void initializeFrequencies(STvertex* v) {
        if (!v) return;
        v->freq.resize(string_count, 0);
        
        if (v->g.empty()) {
            for (int i = 0; i < string_count; i++) {
                if (find(v->string_ids.begin(), v->string_ids.end(), i) != v->string_ids.end()) {
                    v->freq[i] = 1;
                }
            }
        }
        
        for (auto& edge : v->g) {
            initializeFrequencies(edge.second.v);
        }
    }

public:
    SuffixTree(const vector<string>& strings) {
        string_count = strings.size();
        str_starts.clear();
        combined_text.clear();
        int curr_pos = 0;
        
        for (int i = 0; i < strings.size(); i++) {
            str_starts.push_back(curr_pos);
            combined_text += strings[i] + "$" + to_string(i);
            curr_pos += strings[i].length() + 1 + to_string(i).length();
        }

        STvertex* top = new STvertex();
        root = new STvertex();
        
        STedge e;
        e.v = root;
        liscie = 0;
        
        for (int i = 0; i < combined_text.length(); i++) {
            e.r = -i;
            e.l = -i;
            top->g[combined_text[i]] = e;
        }
        
        root->f = top;
        
        e.l = 0;
        e.v = root;
        for (int i = 0; i < combined_text.length(); i++) {
            e.r = i - 1;
            Update(e, combined_text.length());
            e.r++;
            Canonize(e);
        }
    }

    vector<pair<string, vector<int>>> findCommonSubstrings(int threshold) {
        vector<pair<string, vector<int>>> results;
        initializeFrequencies(root);
        dfsCalculateFreq(root, threshold, results);
        
        sort(results.begin(), results.end(),
             [](const auto& a, const auto& b) {
                 return a.first.length() > b.first.length();
             });
        
        return results;
    }

    ~SuffixTree() {
        function<void(STvertex*)> deleteVertex = [&](STvertex* v) {
            if (!v) return;
            for (auto& edge : v->g) {
                deleteVertex(edge.second.v);
            }
            delete v;
        };
        deleteVertex(root);
    }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }
    
    string inputFile = argv[1];
    string outputFile = argv[2];
    
    ifstream inFile(inputFile);
    if (!inFile.is_open()) {
        cerr << "Error opening input file: " << inputFile << endl;
        return 1;
    }
    
    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cerr << "Error opening output file: " << outputFile << endl;
        return 1;
    }
    
    string line;
    vector<string> strings;
    size_t lineCount = 0;
    size_t total_length = 0;
    
    while (getline(inFile, line)) {
        if (!line.empty()) {
            strings.push_back(line);
            total_length += line.length();
            
            lineCount++;
        }
    }
    
    inFile.close();
    
    outFile << "=== Suffix Tree Analysis ===" << endl;
    
    // Algorithm execution
    SuffixTree st(strings);

    auto results = st.findCommonSubstrings(1); // Using threshold 1
    
    const int MAX_RESULTS_TO_SAVE = 10;
    outFile << "Top " << min(MAX_RESULTS_TO_SAVE, (int)results.size()) << " longest common substrings:" << endl;
    
    for (int i = 0; i < min(MAX_RESULTS_TO_SAVE, (int)results.size()); i++) {
        const auto& [substr, freq] = results[i];
        outFile << (i+1) << ". Length: " << substr.length() << ", Substring: \"" << substr << "\"" << endl;
        outFile << "   Frequencies: [";
        for (int j = 0; j < min(20, (int)freq.size()); j++) {
            outFile << freq[j];
            if (j < min(20, (int)freq.size()) - 1) outFile << ",";
        }
        if (freq.size() > 20) outFile << ",...";
        outFile << "]" << endl;
    }
    
    outFile.close();
    
    return 0;
}