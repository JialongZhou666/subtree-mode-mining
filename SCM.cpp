/**
    SCM
    Copyright (C) 2025 Jialong Zhou.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <divsufsort.h>

#define MAXN 2000000
#define MAX_DELTA 20000
#define SEPARATOR_CHAR 1

int Delta;
char** strings;
int* string_lengths;
short* string_ids;

unsigned char* concatenated_text;
int text_length;
int* SA;
int* ISA;
int* LCP;

short* position_to_string_id;
int* position_to_string_pos;

int* sa_to_node;

typedef struct {
    int l;
    int r;
    int lcp;
    int parent;
    int* children;
    int child_count;
    int capacity;
} SuffixTreeNode;

SuffixTreeNode* suffix_tree_nodes;
int node_count;
int node_capacity;

typedef struct {
    int* leaves;
    int* originals;
    int size;
    int capacity;
} LeafList;

LeafList* L;

int* euler;
int* depth_euler;
int* first;
int euler_size;

int block_size;
int block_count;
int* block_min;
int** block_sparse;

int max_patterns;
int*** pattern_table;
int* block_pattern;

typedef struct {
    int st_node;
    int ti;
    int depth_val;
} TreeNode;

typedef struct {
    TreeNode* nodes;
    int size;
    int capacity;
} SingleColorTree;

SingleColorTree* T_color;

typedef struct {
    int st_node;
    int count;
} NodeCount;

typedef struct {
    NodeCount* counts;
    int size;
    int capacity;
} ColorCounts;

ColorCounts* leaf_counts;

typedef struct {
    int color_id;
    int st_node;
    int count;
} ColorNodePair;

typedef struct {
    ColorNodePair* pairs;
    int count;
    int capacity;
} MappedNodes;

MappedNodes* reverse_mapping;

typedef struct {
    int frequency;
    int string_id;
} MaxCount;

MaxCount* max_count;

static inline int integer_log2(int n) {
    int result = 0;
    while (n > 1) {
        n >>= 1;
        result++;
    }
    return result;
}

void readInputStrings(const char* filename, int max_strings) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        exit(1);
    }
    
    char line[10000];
    Delta = 0;
    int total_length = 0;
    
    while (fgets(line, sizeof(line), file)) {
        int len = strlen(line);
        if (len > 0 && line[len-1] == '\n') {
            len--;
            line[len] = '\0';
        }
        if (len > 0) {
            Delta++;
            total_length += len + 1;
            if (max_strings > 0 && Delta >= max_strings) {
                break;
            }
        }
    }
    
    if (Delta == 0) {
        fprintf(stderr, "Error: No strings found\n");
        exit(1);
    }
    
    if (Delta > MAX_DELTA) {
        fprintf(stderr, "Error: Delta exceeds MAX_DELTA\n");
        exit(1);
    }
    
    strings = (char**)malloc(Delta * sizeof(char*));
    string_lengths = (int*)malloc(Delta * sizeof(int));
    string_ids = (short*)malloc(Delta * sizeof(short));
    
    rewind(file);
    int idx = 0;
    while (fgets(line, sizeof(line), file) && idx < Delta) {
        int len = strlen(line);
        if (len > 0 && line[len-1] == '\n') {
            len--;
            line[len] = '\0';
        }
        if (len == 0) continue;
        
        string_lengths[idx] = len;
        string_ids[idx] = idx;
        strings[idx] = (char*)malloc((len + 1) * sizeof(char));
        
        for (int j = 0; j < len; j++) {
            unsigned char ch = (unsigned char)line[j];
            strings[idx][j] = (char)(ch + 2);
        }
        strings[idx][len] = '\0';
        idx++;
    }
    
    fclose(file);
    
    text_length = total_length;
    concatenated_text = (unsigned char*)malloc(text_length * sizeof(unsigned char));
    position_to_string_id = (short*)malloc(text_length * sizeof(short));
    position_to_string_pos = (int*)malloc(text_length * sizeof(int));
    
    int pos = 0;
    for (int i = 0; i < Delta; i++) {
        for (int j = 0; j < string_lengths[i]; j++) {
            concatenated_text[pos] = (unsigned char)strings[i][j];
            position_to_string_id[pos] = string_ids[i];
            position_to_string_pos[pos] = j;
            pos++;
        }
        concatenated_text[pos] = SEPARATOR_CHAR + i;
        position_to_string_id[pos] = -1;
        position_to_string_pos[pos] = -1;
        pos++;
    }
}

void buildSuffixArray() {
    SA = (int*)malloc(text_length * sizeof(int));
    ISA = (int*)malloc(text_length * sizeof(int));
    
    if (divsufsort((sauchar_t*)concatenated_text, (saidx_t*)SA, (saidx_t)text_length) != 0) {
        fprintf(stderr, "Error: SA computation failed\n");
        exit(1);
    }
    
    for (int i = 0; i < text_length; i++) {
        ISA[SA[i]] = i;
    }
}

void buildLCPArray() {
    LCP = (int*)malloc(text_length * sizeof(int));
    LCP[0] = 0;
    
    int j = 0;
    for (int i = 0; i < text_length; i++) {
        if (ISA[i] == 0) {
            j = 0;
            continue;
        }
        
        int prev_suffix = SA[ISA[i] - 1];
        if (i > 0 && ISA[i-1] != 0) {
            j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]] - 1 : 0;
        } else {
            j = 0;
        }
        
        while (i + j < text_length && prev_suffix + j < text_length &&
               concatenated_text[i + j] == concatenated_text[prev_suffix + j]) {
            j++;
        }
        
        LCP[ISA[i]] = j;
    }
}

void addChild(int parent_idx, int child_idx) {
    SuffixTreeNode* parent = &suffix_tree_nodes[parent_idx];
    if (parent->child_count >= parent->capacity) {
        parent->capacity = (parent->capacity == 0) ? 4 : parent->capacity * 2;
        parent->children = (int*)realloc(parent->children, parent->capacity * sizeof(int));
    }
    parent->children[parent->child_count++] = child_idx;
    suffix_tree_nodes[child_idx].parent = parent_idx;
}

int newNode(int l, int r, int lcp) {
    if (node_count >= node_capacity) {
        node_capacity = (node_capacity == 0) ? 1024 : node_capacity * 2;
        suffix_tree_nodes = (SuffixTreeNode*)realloc(suffix_tree_nodes, 
                                                      node_capacity * sizeof(SuffixTreeNode));
    }
    
    int idx = node_count++;
    suffix_tree_nodes[idx].l = l;
    suffix_tree_nodes[idx].r = r;
    suffix_tree_nodes[idx].lcp = lcp;
    suffix_tree_nodes[idx].parent = -1;
    suffix_tree_nodes[idx].child_count = 0;
    suffix_tree_nodes[idx].capacity = 0;
    suffix_tree_nodes[idx].children = NULL;
    
    return idx;
}

void buildSuffixTreeFromLCP() {
    node_count = 0;
    node_capacity = 1024;
    suffix_tree_nodes = (SuffixTreeNode*)malloc(node_capacity * sizeof(SuffixTreeNode));
    
    sa_to_node = (int*)malloc(text_length * sizeof(int));
    for (int i = 0; i < text_length; i++) {
        sa_to_node[i] = -1;
    }
    
    int* stack = (int*)malloc(text_length * sizeof(int));
    int stack_top = -1;
    
    int root = newNode(0, text_length - 1, 0);
    stack[++stack_top] = root;
    
    int x = -1;
    
    for (int i = 1; i < text_length; i++) {
        int l = i - 1;
        
        while (stack_top >= 0 && LCP[i] < suffix_tree_nodes[stack[stack_top]].lcp) {
            x = stack[stack_top--];
            suffix_tree_nodes[x].r = i - 1;
            
            if (suffix_tree_nodes[x].l == suffix_tree_nodes[x].r) {
                sa_to_node[suffix_tree_nodes[x].l] = x;
            }
            
            l = suffix_tree_nodes[x].l;
            
            if (stack_top >= 0 && LCP[i] <= suffix_tree_nodes[stack[stack_top]].lcp) {
                addChild(stack[stack_top], x);
                x = -1;
            }
        }
        
        if (stack_top >= 0 && LCP[i] > suffix_tree_nodes[stack[stack_top]].lcp) {
            int new_node = newNode(l, i - 1, LCP[i]);
            stack[++stack_top] = new_node;
            
            if (x != -1) {
                addChild(new_node, x);
                x = -1;
            }
        }
    }
    
    while (stack_top > 0) {
        x = stack[stack_top--];
        suffix_tree_nodes[x].r = text_length - 1;
        
        if (suffix_tree_nodes[x].l == suffix_tree_nodes[x].r) {
            sa_to_node[suffix_tree_nodes[x].l] = x;
        }
        
        addChild(stack[stack_top], x);
    }
    
    if (stack_top == 0) {
        suffix_tree_nodes[stack[stack_top]].r = text_length - 1;
    }
    
    free(stack);
}

void buildLeafLists() {
    L = (LeafList*)malloc(Delta * sizeof(LeafList));
    for (int i = 0; i < Delta; i++) {
        L[i].size = 0;
        L[i].capacity = 16;
        L[i].leaves = (int*)malloc(L[i].capacity * sizeof(int));
        L[i].originals = (int*)malloc(L[i].capacity * sizeof(int));
    }
    
    for (int i = 0; i < text_length; i++) {
        int pos = SA[i];
        int str_id = position_to_string_id[pos];
        
        if (str_id < 0 || str_id >= Delta) continue;
        
        LeafList* list = &L[str_id];
        if (list->size >= list->capacity) {
            list->capacity <<= 1;
            list->leaves = (int*)realloc(list->leaves, list->capacity * sizeof(int));
            list->originals = (int*)realloc(list->originals, list->capacity * sizeof(int));
        }
        
        int newLeafId = list->size;
        list->leaves[newLeafId] = newLeafId;
        list->originals[newLeafId] = i;
        list->size++;
    }
}

void step1_phaseA() {
    buildLeafLists();
}

void dfs_euler_suffix_tree(int node_idx, int d) {
    first[node_idx] = euler_size;
    euler[euler_size] = node_idx;
    depth_euler[euler_size++] = d;
    
    SuffixTreeNode* node = &suffix_tree_nodes[node_idx];
    for (int i = 0; i < node->child_count; i++) {
        int child_idx = node->children[i];
        dfs_euler_suffix_tree(child_idx, d + 1);
        euler[euler_size] = node_idx;
        depth_euler[euler_size++] = d;
    }
}

static inline int compute_block_pattern_suffix(int block_id) {
    int start = block_id * block_size;
    int end = (block_id + 1) * block_size;
    if (end > euler_size) end = euler_size;
    
    int pattern = 0;
    for (int i = start + 1; i < end; i++) {
        if (depth_euler[i] > depth_euler[i-1]) {
            pattern |= (1 << (i - start - 1));
        }
    }
    return pattern;
}

static inline int queryLCA_suffix_tree(int u, int v) {
    if (u == v) return u;
    
    int pos_u = first[u];
    int pos_v = first[v];
    
    if (pos_u > pos_v) {
        int tmp = pos_u;
        pos_u = pos_v;
        pos_v = tmp;
    }
    
    int block_u = pos_u / block_size;
    int block_v = pos_v / block_size;
    
    if (block_u == block_v) {
        int pattern = block_pattern[block_u];
        int offset_u = pos_u % block_size;
        int offset_v = pos_v % block_size;
        return euler[block_u * block_size + pattern_table[pattern][offset_u][offset_v]];
    }
    
    int min_pos = pos_u;
    int min_depth = depth_euler[pos_u];
    
    if (pos_u % block_size != 0) {
        int pattern = block_pattern[block_u];
        int offset_u = pos_u % block_size;
        int offset_end = ((block_u + 1) * block_size <= euler_size) ? 
                         block_size - 1 : (euler_size - 1) % block_size;
        
        int cand_pos = block_u * block_size + pattern_table[pattern][offset_u][offset_end];
        if (depth_euler[cand_pos] < min_depth) {
            min_depth = depth_euler[cand_pos];
            min_pos = cand_pos;
        }
    }
    
    int left_block = block_u + 1;
    int right_block = block_v - 1;
    
    if (left_block <= right_block) {
        int k = integer_log2(right_block - left_block + 1);
        
        int cand1 = block_sparse[left_block][k];
        int cand2 = block_sparse[right_block - (1 << k) + 1][k];
        
        int cand_pos = (depth_euler[cand1] <= depth_euler[cand2]) ? cand1 : cand2;
        if (depth_euler[cand_pos] < min_depth) {
            min_depth = depth_euler[cand_pos];
            min_pos = cand_pos;
        }
    }
    
    if (pos_v >= block_v * block_size) {
        int pattern = block_pattern[block_v];
        int offset_v = pos_v % block_size;
        
        int cand_pos = block_v * block_size + pattern_table[pattern][0][offset_v];
        if (depth_euler[cand_pos] < min_depth) {
            min_pos = cand_pos;
        }
    }
    
    return euler[min_pos];
}

void precompute_pattern_table_suffix() {
    max_patterns = 1 << (block_size - 1);
    pattern_table = (int***)malloc(max_patterns * sizeof(int**));
    
    int* depths = (int*)malloc(block_size * sizeof(int));
    
    for (int pattern = 0; pattern < max_patterns; pattern++) {
        pattern_table[pattern] = (int**)malloc(block_size * sizeof(int*));
        for (int i = 0; i < block_size; i++) {
            pattern_table[pattern][i] = (int*)malloc(block_size * sizeof(int));
        }
        
        depths[0] = 0;
        for (int i = 1; i < block_size; i++) {
            depths[i] = depths[i-1] + (((pattern >> (i-1)) & 1) ? 1 : -1);
        }
        
        for (int i = 0; i < block_size; i++) {
            int min_idx = i;
            int min_depth = depths[i];
            pattern_table[pattern][i][i] = i;
            
            for (int j = i + 1; j < block_size; j++) {
                if (depths[j] < min_depth) {
                    min_depth = depths[j];
                    min_idx = j;
                }
                pattern_table[pattern][i][j] = min_idx;
            }
        }
    }
    
    free(depths);
}

void build_block_sparse_table_suffix() {
    int log_bc = integer_log2(block_count) + 1;
    
    block_sparse = (int**)malloc(block_count * sizeof(int*));
    for (int i = 0; i < block_count; i++) {
        block_sparse[i] = (int*)malloc(log_bc * sizeof(int));
        block_sparse[i][0] = block_min[i];
    }
    
    for (int j = 1; (1 << j) <= block_count; j++) {
        for (int i = 0; i + (1 << j) <= block_count; i++) {
            int left = block_sparse[i][j-1];
            int right = block_sparse[i + (1 << (j-1))][j-1];
            block_sparse[i][j] = (depth_euler[left] <= depth_euler[right]) ? left : right;
        }
    }
}

void buildLCA_suffix_tree() {
    euler = (int*)malloc((2 * node_count) * sizeof(int));
    depth_euler = (int*)malloc((2 * node_count) * sizeof(int));
    first = (int*)malloc(node_count * sizeof(int));
    
    euler_size = 0;
    dfs_euler_suffix_tree(0, 0);
    
    block_size = integer_log2(euler_size) >> 1;
    if (block_size < 1) block_size = 1;
    if (block_size > 31) block_size = 31;
    
    block_count = (euler_size + block_size - 1) / block_size;
    
    precompute_pattern_table_suffix();
    
    block_min = (int*)malloc(block_count * sizeof(int));
    block_pattern = (int*)malloc(block_count * sizeof(int));
    
    for (int b = 0; b < block_count; b++) {
        int start = b * block_size;
        int end = (b + 1) * block_size;
        if (end > euler_size) end = euler_size;
        
        int min_pos = start;
        int min_depth_val = depth_euler[start];
        for (int i = start + 1; i < end; i++) {
            if (depth_euler[i] < min_depth_val) {
                min_depth_val = depth_euler[i];
                min_pos = i;
            }
        }
        block_min[b] = min_pos;
        block_pattern[b] = compute_block_pattern_suffix(b);
    }
    
    build_block_sparse_table_suffix();
}

void buildSingleColorTree_Kasai(int c) {
    int leaf_count = L[c].size;
    
    SingleColorTree* tree = &T_color[c];
    tree->size = 0;
    tree->capacity = leaf_count * 2;
    tree->nodes = (TreeNode*)malloc(tree->capacity * sizeof(TreeNode));
    
    if (leaf_count == 0) return;
    
    int* st_nodes = (int*)malloc(leaf_count * sizeof(int));
    for (int i = 0; i < leaf_count; i++) {
        int sa_pos = L[c].originals[i];
        st_nodes[i] = sa_to_node[sa_pos];
    }
    
    if (leaf_count == 1) {
        if (st_nodes[0] >= 0) {
            tree->nodes[0] = (TreeNode){st_nodes[0], 0, suffix_tree_nodes[st_nodes[0]].lcp};
            tree->size = 1;
        }
        free(st_nodes);
        return;
    }
    
    TreeNode* stack = (TreeNode*)malloc(tree->capacity * sizeof(TreeNode));
    int stack_top = -1;
    
    stack[++stack_top] = (TreeNode){-1, -1, -1};
    
    for (int k = 0; k <= leaf_count; k++) {
        int lca_node, lca_depth;
        
        if (k < leaf_count) {
            if (k == 0) {
                if (st_nodes[0] >= 0) {
                    lca_node = st_nodes[0];
                    lca_depth = suffix_tree_nodes[lca_node].lcp;
                } else {
                    lca_node = -1;
                    lca_depth = -1;
                }
            } else {
                int u = st_nodes[k - 1];
                int v = st_nodes[k];
                if (u >= 0 && v >= 0) {
                    lca_node = queryLCA_suffix_tree(u, v);
                    lca_depth = suffix_tree_nodes[lca_node].lcp;
                } else {
                    lca_node = -1;
                    lca_depth = -1;
                }
            }
        } else {
            lca_node = -1;
            lca_depth = -1;
        }
        
        while (stack_top >= 0 && stack[stack_top].depth_val > lca_depth) {
            TreeNode node = stack[stack_top--];
            
            if (tree->size >= tree->capacity) {
                tree->capacity *= 2;
                tree->nodes = (TreeNode*)realloc(tree->nodes, tree->capacity * sizeof(TreeNode));
                stack = (TreeNode*)realloc(stack, tree->capacity * sizeof(TreeNode));
            }
            tree->nodes[tree->size++] = node;
        }
        
        if (k < leaf_count && lca_node >= 0) {
            if (stack_top < 0 || stack[stack_top].st_node != lca_node) {
                stack[++stack_top] = (TreeNode){lca_node, k, lca_depth};
            }
            
            if (st_nodes[k] >= 0) {
                stack[++stack_top] = (TreeNode){st_nodes[k], k, suffix_tree_nodes[st_nodes[k]].lcp};
            }
        }
    }
    
    free(stack);
    free(st_nodes);
}

void step1_phaseB() {
    buildLCA_suffix_tree();
    
    T_color = (SingleColorTree*)malloc(Delta * sizeof(SingleColorTree));
    for (int c = 0; c < Delta; c++) {
        T_color[c].size = 0;
        T_color[c].capacity = 0;
        T_color[c].nodes = NULL;
    }
    
    for (int c = 0; c < Delta; c++) {
        if (L[c].size == 0) continue;
        buildSingleColorTree_Kasai(c);
    }
}

void countLeafDescendantsForTree(int c) {
    SingleColorTree* tree = &T_color[c];
    
    ColorCounts* counts = &leaf_counts[c];
    counts->size = 0;
    counts->capacity = tree->size;
    counts->counts = (NodeCount*)malloc(counts->capacity * sizeof(NodeCount));
    
    if (tree->size == 0) return;
    
    int* node_counts = (int*)calloc(node_count, sizeof(int));
    bool* in_tree = (bool*)calloc(node_count, sizeof(bool));
    
    for (int i = 0; i < tree->size; i++) {
        in_tree[tree->nodes[i].st_node] = true;
    }
    
    for (int i = 0; i < tree->size; i++) {
        int st_node = tree->nodes[i].st_node;
        SuffixTreeNode* node = &suffix_tree_nodes[st_node];
        
        if (node->l == node->r) {
            int pos = SA[node->l];
            int str_id = position_to_string_id[pos];
            if (str_id == c) {
                node_counts[st_node] = 1;
            } else {
                node_counts[st_node] = 0;
            }
        } else {
            int count = 0;
            for (int j = 0; j < node->child_count; j++) {
                int child = node->children[j];
                if (in_tree[child]) {
                    count += node_counts[child];
                }
            }
            node_counts[st_node] = count;
        }
        
        counts->counts[counts->size++] = (NodeCount){st_node, node_counts[st_node]};
    }
    
    free(node_counts);
    free(in_tree);
}

void step2_countColors() {
    leaf_counts = (ColorCounts*)malloc(Delta * sizeof(ColorCounts));
    
    for (int c = 0; c < Delta; c++) {
        if (L[c].size == 0) continue;
        countLeafDescendantsForTree(c);
    }
}

void buildReverseMapping() {
    reverse_mapping = (MappedNodes*)malloc(node_count * sizeof(MappedNodes));
    
    for (int v = 0; v < node_count; v++) {
        reverse_mapping[v].count = 0;
        reverse_mapping[v].capacity = 4;
        reverse_mapping[v].pairs = (ColorNodePair*)malloc(
            reverse_mapping[v].capacity * sizeof(ColorNodePair));
    }
    
    for (int c = 0; c < Delta; c++) {
        if (L[c].size == 0) continue;
        
        ColorCounts* counts = &leaf_counts[c];
        
        for (int i = 0; i < counts->size; i++) {
            int v = counts->counts[i].st_node;
            int count = counts->counts[i].count;
            
            if (v < 0 || v >= node_count) continue;
            
            SuffixTreeNode* node = &suffix_tree_nodes[v];
            if (node->l == node->r) continue;
            
            MappedNodes* mapping = &reverse_mapping[v];
            
            if (mapping->count >= mapping->capacity) {
                mapping->capacity *= 2;
                mapping->pairs = (ColorNodePair*)realloc(
                    mapping->pairs, 
                    mapping->capacity * sizeof(ColorNodePair));
            }
            
            mapping->pairs[mapping->count].color_id = c;
            mapping->pairs[mapping->count].st_node = v;
            mapping->pairs[mapping->count].count = count;
            mapping->count++;
        }
    }
}

void mergeCounts_DFS(int v) {
    SuffixTreeNode* node = &suffix_tree_nodes[v];
    
    if (node->l == node->r) {
        int pos = SA[node->l];
        int str_id = position_to_string_id[pos];
        if (str_id >= 0 && str_id < Delta) {
            max_count[v].frequency = 1;
            max_count[v].string_id = str_id;
        } else {
            max_count[v].frequency = 0;
            max_count[v].string_id = -1;
        }
        return;
    }
    
    for (int i = 0; i < node->child_count; i++) {
        mergeCounts_DFS(node->children[i]);
    }
    
    int max_freq = (node->child_count > 0) ? max_count[node->children[0]].frequency : 0;
    int max_col = (node->child_count > 0) ? max_count[node->children[0]].string_id : -1;
    
    for (int i = 1; i < node->child_count; i++) {
        int child = node->children[i];
        if (max_count[child].frequency > max_freq) {
            max_freq = max_count[child].frequency;
            max_col = max_count[child].string_id;
        }
    }
    
    MappedNodes* mapping = &reverse_mapping[v];
    
    for (int i = 0; i < mapping->count; i++) {
        int count = mapping->pairs[i].count;
        int c = mapping->pairs[i].color_id;
        
        if (count > max_freq) {
            max_freq = count;
            max_col = c;
        }
    }
    
    max_count[v].frequency = max_freq;
    max_count[v].string_id = max_col;
}

void step3_mergeCounts() {
    buildReverseMapping();
    max_count = (MaxCount*)malloc(node_count * sizeof(MaxCount));
    mergeCounts_DFS(0);
}

void cleanup() {
    if (strings) {
        for (int i = 0; i < Delta; i++) {
            if (strings[i]) {
                free(strings[i]);
            }
        }
        free(strings);
        strings = NULL;
    }
    if (string_lengths) {
        free(string_lengths);
        string_lengths = NULL;
    }
    if (string_ids) {
        free(string_ids);
        string_ids = NULL;
    }
    if (concatenated_text) {
        free(concatenated_text);
        concatenated_text = NULL;
    }
    if (SA) {
        free(SA);
        SA = NULL;
    }
    if (ISA) {
        free(ISA);
        ISA = NULL;
    }
    if (LCP) {
        free(LCP);
        LCP = NULL;
    }
    if (position_to_string_id) {
        free(position_to_string_id);
        position_to_string_id = NULL;
    }
    if (position_to_string_pos) {
        free(position_to_string_pos);
        position_to_string_pos = NULL;
    }
    if (sa_to_node) {
        free(sa_to_node);
        sa_to_node = NULL;
    }
    
    if (suffix_tree_nodes) {
        for (int i = 0; i < node_count; i++) {
            if (suffix_tree_nodes[i].children) {
                free(suffix_tree_nodes[i].children);
            }
        }
        free(suffix_tree_nodes);
        suffix_tree_nodes = NULL;
    }
    if (max_count) {
        free(max_count);
        max_count = NULL;
    }
    if (reverse_mapping) {
        for (int i = 0; i < node_count; i++) {
            if (reverse_mapping[i].pairs) {
                free(reverse_mapping[i].pairs);
            }
        }
        free(reverse_mapping);
        reverse_mapping = NULL;
    }
    if (L) {
        for (int i = 0; i < Delta; i++) {
            if (L[i].leaves) free(L[i].leaves);
            if (L[i].originals) free(L[i].originals);
        }
        free(L);
        L = NULL;
    }
    if (T_color) {
        for (int i = 0; i < Delta; i++) {
            if (T_color[i].nodes) free(T_color[i].nodes);
        }
        free(T_color);
        T_color = NULL;
    }
    if (leaf_counts) {
        for (int i = 0; i < Delta; i++) {
            if (leaf_counts[i].counts) free(leaf_counts[i].counts);
        }
        free(leaf_counts);
        leaf_counts = NULL;
    }
    if (euler) {
        free(euler);
        euler = NULL;
    }
    if (depth_euler) {
        free(depth_euler);
        depth_euler = NULL;
    }
    if (first) {
        free(first);
        first = NULL;
    }
    if (block_min) {
        free(block_min);
        block_min = NULL;
    }
    if (block_pattern) {
        free(block_pattern);
        block_pattern = NULL;
    }
    if (block_sparse) {
        for (int i = 0; i < block_count; i++) {
            if (block_sparse[i]) free(block_sparse[i]);
        }
        free(block_sparse);
        block_sparse = NULL;
    }
    if (pattern_table) {
        for (int p = 0; p < max_patterns; p++) {
            if (pattern_table[p]) {
                for (int i = 0; i < block_size; i++) {
                    if (pattern_table[p][i]) free(pattern_table[p][i]);
                }
                free(pattern_table[p]);
            }
        }
        free(pattern_table);
        pattern_table = NULL;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }
    
    const char* input_file = argv[1];
    const char* output_file = argv[2];
    
    readInputStrings(input_file, -1);
    buildSuffixArray();
    buildLCPArray();
    buildSuffixTreeFromLCP();
    
    step1_phaseA();
    step1_phaseB();
    step2_countColors();
    step3_mergeCounts();
    
    int root_mode_str = max_count[0].string_id;
    int root_mode_freq = max_count[0].frequency;
    
    FILE* out = fopen(output_file, "w");
    if (out) {
        fprintf(out, "StringID: %d\nFrequency: %d\nNodes: %d\n", 
                root_mode_str, root_mode_freq, node_count);
        fclose(out);
    }
    
    cleanup();
    return 0;
}
