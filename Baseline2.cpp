#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <functional>
#include <cmath>
#include <climits>
#include <set>

using namespace std;

class ExactMode_Sqrt {
private:
    int* array;         
    int len;            
    int s;              
    int t;              
    int delta;          
    int* A_Prime;       
    int** Qa;           
    int* QaSize;        
    
    struct BlockTable {
        int** S;
        int** S_freq;
        int** S_min;
        int** S_min_freq;
        int block_start, block_end;
        
        BlockTable() : S(nullptr), S_freq(nullptr), S_min(nullptr), S_min_freq(nullptr) {}
        ~BlockTable() {
            if (S) {
                int size = block_end - block_start + 1;
                for (int i = 0; i < size; i++) {
                    delete[] S[i];
                    delete[] S_freq[i];
                    delete[] S_min[i];
                    delete[] S_min_freq[i];
                }
                delete[] S;
                delete[] S_freq;
                delete[] S_min;
                delete[] S_min_freq;
            }
        }
    };
    
    vector<BlockTable*> block_tables;
    int blocks_per_table;
    
public:
    ExactMode_Sqrt(int* array, int len, int s) {
        this->array = array;
        this->len = len;
        this->s = s;
        this->A_Prime = nullptr;
        this->Qa = nullptr;
        this->QaSize = nullptr;
        blocks_per_table = min(s, 200);
    }
    
    ~ExactMode_Sqrt() {
        if (QaSize) delete[] QaSize;
        if (Qa) {
            for (int i = 0; i < delta; i++) {
                if (Qa[i]) delete[] Qa[i];
            }
            delete[] Qa;
        }
        for (auto* table : block_tables) {
            delete table;
        }
        if (A_Prime) delete[] A_Prime;
    }
    
    void construct() {
        buildPositionIndex();
        buildBlockTables();
    }
    
private:
    void buildPositionIndex() {
        delta = 0;
        for (int i = 0; i < len; i++) {
            delta = max(delta, array[i]);
        }
        
        A_Prime = new int[len];
        QaSize = new int[delta]();
        
        for (int i = 0; i < len; i++) {
            QaSize[array[i]-1]++;
        }
        
        Qa = new int*[delta];
        for (int i = 0; i < delta; i++) {
            if (QaSize[i] > 0) {
                Qa[i] = new int[QaSize[i]];
            } else {
                Qa[i] = nullptr;
            }
        }
        
        vector<int> index(delta, 0);
        for (int i = 0; i < len; i++) {
            A_Prime[i] = index[array[i]-1]++;
            Qa[array[i]-1][A_Prime[i]] = i + 1;
        }
        
        t = (len + s - 1) / s;
    }
    
    void buildBlockTables() {
        for (int chunk_start = 0; chunk_start < s; chunk_start += blocks_per_table) {
            int chunk_end = min(chunk_start + blocks_per_table - 1, s - 1);
            
            BlockTable* table = new BlockTable();
            table->block_start = chunk_start;
            table->block_end = chunk_end;
            
            int chunk_size = chunk_end - chunk_start + 1;
            
            table->S = new int*[chunk_size];
            table->S_freq = new int*[chunk_size];
            table->S_min = new int*[chunk_size];
            table->S_min_freq = new int*[chunk_size];
            
            for (int i = 0; i < chunk_size; i++) {
                table->S[i] = new int[s]();
                table->S_freq[i] = new int[s]();
                table->S_min[i] = new int[s]();
                table->S_min_freq[i] = new int[s]();
                fill(table->S_min_freq[i], table->S_min_freq[i] + s, INT_MAX);
            }
            
            computeBlockTable(table, chunk_start, chunk_end);
            block_tables.push_back(table);
        }
    }
    
    void computeBlockTable(BlockTable* table, int chunk_start, int chunk_end) {
        vector<int> countArray(delta + 1);
        
        for (int i = chunk_start; i <= chunk_end; i++) {
            for (int j = i; j < s; j++) {
                fill(countArray.begin(), countArray.end(), 0);
                
                int max_freq = 0, mode = 0;
                int min_freq = INT_MAX, min_element = 0;
                
                int start_pos = i * t;
                int end_pos = min((j + 1) * t - 1, len - 1);
                
                for (int pos = start_pos; pos <= end_pos; pos++) {
                    countArray[array[pos]]++;
                    
                    if (countArray[array[pos]] > max_freq) {
                        max_freq = countArray[array[pos]];
                        mode = array[pos];
                    }
                }
                
                for (int k = 1; k <= delta; k++) {
                    if (countArray[k] > 0 && countArray[k] < min_freq) {
                        min_freq = countArray[k];
                        min_element = k;
                    }
                }
                
                int local_i = i - chunk_start;
                table->S[local_i][j] = mode;
                table->S_freq[local_i][j] = max_freq;
                table->S_min[local_i][j] = min_element;
                table->S_min_freq[local_i][j] = min_freq;
            }
        }
    }
    
    tuple<int, int, int, int> getPrecomputedValues(int bi, int bj) {
        if (bi > bj) return make_tuple(0, 0, 0, INT_MAX);
        
        BlockTable* table = nullptr;
        int local_i = -1;
        
        for (auto* t : block_tables) {
            if (bi >= t->block_start && bi <= t->block_end) {
                table = t;
                local_i = bi - t->block_start;
                break;
            }
        }
        
        if (!table || local_i < 0) {
            return make_tuple(0, 0, 0, INT_MAX);
        }
        
        return make_tuple(
            table->S[local_i][bj],
            table->S_freq[local_i][bj],
            table->S_min[local_i][bj],
            table->S_min_freq[local_i][bj]
        );
    }
    
public:
    int getElementFreq(int element, int start_index, int end_index) {
        if (element < 1 || element > delta || !Qa[element-1]) return 0;
        
        int* positions = Qa[element - 1];
        int size = QaSize[element - 1];
        
        int left = 0, right = size - 1;
        int first_valid = -1, last_valid = -1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            if (positions[mid] >= start_index) {
                first_valid = mid;
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        
        if (first_valid == -1) return 0;
        
        left = first_valid;
        right = size - 1;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (positions[mid] <= end_index) {
                last_valid = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        
        return (last_valid >= first_valid) ? (last_valid - first_valid + 1) : 0;
    }
    
    pair<int, int> findMode(int start_index, int end_index) {
        start_index--; end_index--;
        
        if (start_index > end_index || start_index < 0 || end_index >= len) {
            return make_pair(0, 0);
        }
        
        int bi = start_index / t;
        int bj = end_index / t;
        
        int c = 0, freq_c = 0;
        
        if (bi + 1 <= bj - 1) {
            auto [mode, mode_freq, min_elem, min_freq] = getPrecomputedValues(bi + 1, bj - 1);
            c = mode;
            freq_c = mode_freq;
        }
        
        set<int> candidates;
        
        int prefix_end = min(bi * t + t - 1, end_index);
        for (int i = start_index; i <= prefix_end; i++) {
            candidates.insert(array[i]);
        }
        
        int suffix_start = max(bj * t, start_index);
        if (suffix_start <= end_index) {
            for (int i = suffix_start; i <= end_index; i++) {
                candidates.insert(array[i]);
            }
        }
        
        for (int element : candidates) {
            int total_freq = getElementFreq(element, start_index + 1, end_index + 1);
            if (total_freq > freq_c) {
                c = element;
                freq_c = total_freq;
            }
        }
        
        return make_pair(c, freq_c);
    }
    
    pair<int, int> findMinFreq(int start_index, int end_index) {
        start_index--; end_index--;
        
        if (start_index > end_index || start_index < 0 || end_index >= len) {
            return make_pair(0, INT_MAX);
        }
        
        int bi = start_index / t;
        int bj = end_index / t;
        
        int min_element = 0, min_freq = INT_MAX;
        
        if (bi + 1 <= bj - 1) {
            auto [mode, mode_freq, min_elem, min_f] = getPrecomputedValues(bi + 1, bj - 1);
            min_element = min_elem;
            min_freq = min_f;
        }
        
        set<int> candidates;
        
        int prefix_end = min(bi * t + t - 1, end_index);
        for (int i = start_index; i <= prefix_end; i++) {
            candidates.insert(array[i]);
        }
        
        int suffix_start = max(bj * t, start_index);
        if (suffix_start <= end_index) {
            for (int i = suffix_start; i <= end_index; i++) {
                candidates.insert(array[i]);
            }
        }
        
        for (int element : candidates) {
            int total_freq = getElementFreq(element, start_index + 1, end_index + 1);
            if (total_freq > 0 && total_freq < min_freq) {
                min_freq = total_freq;
                min_element = element;
            }
        }
        
        return make_pair(min_element, min_freq);
    }
};

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
    vector<int> string_ids;
    int numer;
    int converted_id;  
    
    STvertex() : f(nullptr), parent(nullptr), edge_char(0), numer(-1), converted_id(-1) {}
};

class SuffixTreeToSubtreeConverter {
private:
    STvertex* root;
    string combined_text;
    vector<int> str_starts;
    int string_count;
    int liscie;
    vector<int> converted_colors;  
    vector<pair<int,int>> node_ranges;  
    int total_nodes;
    
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
    
    void convertToSubtreeProblem(STvertex* v, int& current_color_index) {
        if (!v) return;
        
        v->converted_id = total_nodes++;
        int start_index = converted_colors.size();
        
        if (v->g.empty()) {
            for (int string_id : v->string_ids) {
                converted_colors.push_back(string_id + 1);
            }
        }
        
        for (auto& edge : v->g) {
            convertToSubtreeProblem(edge.second.v, current_color_index);
        }
        
        int end_index = converted_colors.size() - 1;
        
        if (start_index <= end_index) {
            node_ranges[v->converted_id] = {start_index + 1, end_index + 1};
        } else {
            node_ranges[v->converted_id] = {-1, -1};
        }
    }
    
    string getNodePath(STvertex* v) {
        if (v == root) return "";
        
        string path = "";
        STvertex* curr = v;
        
        while (curr != root && curr->parent) {
            STvertex* parent = curr->parent;
            for (auto& edge : parent->g) {
                if (edge.second.v == curr) {
                    string edgeStr = combined_text.substr(edge.second.l, 
                                                        edge.second.r - edge.second.l + 1);
                    path = edgeStr + path;
                    break;
                }
            }
            curr = curr->parent;
        }
        
        return path;
    }

public:
    SuffixTreeToSubtreeConverter(const vector<string>& strings) {
        string_count = strings.size();
        total_nodes = 0;
        str_starts.clear();
        combined_text.clear();
        int curr_pos = 0;
        
        combined_text.reserve(strings.size() * 110);
        
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
        
        int color_index = 0;
        node_ranges.resize(combined_text.length() * 2);
        convertToSubtreeProblem(root, color_index);
        node_ranges.resize(total_nodes);
    }
    
    void solveWithExactMode(ofstream& output_file) {
        if (converted_colors.empty()) return;
        
        int* colors_array = new int[converted_colors.size()];
        for (int i = 0; i < converted_colors.size(); i++) {
            colors_array[i] = converted_colors[i];
        }
        
        int s = ceil(sqrt(converted_colors.size()));
        
        ExactMode_Sqrt* query_structure = new ExactMode_Sqrt(colors_array, converted_colors.size(), s);
        query_structure->construct();
        
        vector<tuple<string, vector<int>, int, int>> results;
        analyzeAllNodes(root, query_structure, results);
        
        sort(results.begin(), results.end(),
            [](const auto& a, const auto& b) {
                return get<0>(a).length() > get<0>(b).length();
            });
        
        output_file << "=== Suffix Tree Analysis ===" << endl;
        
        const int MAX_RESULTS_TO_SAVE = 10;
        output_file << "Top " << min(MAX_RESULTS_TO_SAVE, (int)results.size()) << " longest common substrings:" << endl;
        
        int displayed = 0;
        for (const auto& [substr, freq_array, max_freq, min_freq] : results) {
            if (substr.find('$') != string::npos || substr.empty()) continue;
            
            output_file << (displayed+1) << ". Length: " << substr.length() << ", Substring: \"" << substr << "\"" << endl;
            output_file << "   Frequencies: [";
            for (int i = 0; i < min(20, (int)freq_array.size()); i++) {
                output_file << freq_array[i];
                if (i < min(20, (int)freq_array.size()) - 1) output_file << ",";
            }
            if (freq_array.size() > 20) output_file << ",...";
            output_file << "]" << endl;
            
            if (++displayed >= MAX_RESULTS_TO_SAVE) break;
        }
        
        delete[] colors_array;
        delete query_structure;
    }
    
    void analyzeAllNodes(STvertex* v, ExactMode_Sqrt* query_structure, 
                        vector<tuple<string, vector<int>, int, int>>& results) {
        if (!v || v->converted_id < 0 || v->converted_id >= node_ranges.size()) return;
        
        auto [start, end] = node_ranges[v->converted_id];
        
        if (start > 0 && end > 0 && start <= end) {
            vector<int> frequencies(string_count, 0);
            for (int i = start; i <= end; i++) {
                int color = converted_colors[i - 1];
                if (color >= 1 && color <= string_count) {
                    frequencies[color - 1]++;
                }
            }
            
            auto mode_result = query_structure->findMode(start, end);
            auto min_freq_result = query_structure->findMinFreq(start, end);
            
            int max_freq = mode_result.second;
            int min_freq = min_freq_result.second;
            
            string path = getNodePath(v);
            results.push_back(make_tuple(path, frequencies, max_freq, min_freq));
        }
        
        for (auto& edge : v->g) {
            analyzeAllNodes(edge.second.v, query_structure, results);
        }
    }
    
    ~SuffixTreeToSubtreeConverter() {
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

vector<string> readStringsFromFile(const string& filename) {
    vector<string> strings;
    ifstream inFile(filename);
    
    if (!inFile.is_open()) {
        cerr << "Error: Cannot open file '" << filename << "'" << endl;
        return strings;
    }
    
    string line;
    while (getline(inFile, line)) {
        if (!line.empty()) {
            strings.push_back(line);
        }
    }
    
    inFile.close();
    return strings;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input_filename> <output_filename>" << endl;
        return 1;
    }
    
    string input_filename = argv[1];
    string output_filename = argv[2];
    
    ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        cerr << "Error: Cannot open output file '" << output_filename << "'" << endl;
        return 1;
    }
    
    vector<string> input_strings = readStringsFromFile(input_filename);
    
    if (input_strings.empty()) {
        cerr << "Error: No valid strings read from file!" << endl;
        return 1;
    }
    
    SuffixTreeToSubtreeConverter converter(input_strings);
    converter.solveWithExactMode(output_file);
    
    output_file.close();
    return 0;
}
