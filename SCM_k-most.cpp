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

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_DIVSUFSORT
    #include <divsufsort.h>
    #ifdef _USE_64
        #include <divsufsort64.h>
    #endif
#endif

#define MAXN 2000000
#define MAX_DELTA 2000
#define MAX_K 512
#define SEPARATOR_CHAR 1

int k_value = 1;

int Delta;
char** strings;
int* string_lengths;
int* string_ids;

unsigned char* concatenated_text;
int text_length;
int* SA;
int* ISA;
int* LCP;

int* position_to_string_id;
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
} ColorFreq;

typedef struct {
    ColorFreq* colors;
    int count;
} KMostFrequent;

KMostFrequent* k_most_frequent;

int* shared_freq_map;
ColorFreq* shared_temp;
int* shared_active_colors;
int* shared_count_array;
ColorFreq* shared_output_array;

static inline int integer_log2(int n) {
    int result = 0;
    while (n > 1) {
        n >>= 1;
        result++;
    }
    return result;
}

int compare_color_freq(const void* a, const void* b) {
    ColorFreq* ca = (ColorFreq*)a;
    ColorFreq* cb = (ColorFreq*)b;
    if (ca->frequency != cb->frequency) {
        return cb->frequency - ca->frequency;
    }
    return ca->string_id - cb->string_id;
}

int compare_color_freq_by_id(const void* a, const void* b) {
    ColorFreq* ca = (ColorFreq*)a;
    ColorFreq* cb = (ColorFreq*)b;
    return ca->string_id - cb->string_id;
}

int partition_desc(ColorFreq* arr, int left, int right, int pivot_idx) {

    ColorFreq pivot_val = arr[pivot_idx];
    
    ColorFreq tmp = arr[pivot_idx];
    arr[pivot_idx] = arr[right];
    arr[right] = tmp;
    
    int store_idx = left;
    for (int i = left; i < right; i++) {

        if (compare_color_freq(&arr[i], &pivot_val) <= 0) {
            tmp = arr[store_idx];
            arr[store_idx] = arr[i];
            arr[i] = tmp;
            store_idx++;
        }
    }
    
    tmp = arr[store_idx];
    arr[store_idx] = arr[right];
    arr[right] = tmp;
    
    return store_idx;
}

int quickselect(ColorFreq* arr, int left, int right, int k) {
    if (left == right) {
        return left;
    }
    
    int pivot_idx = left + (right - left) / 2;
    pivot_idx = partition_desc(arr, left, right, pivot_idx);
    
    if (k == pivot_idx) {
        return k;
    } else if (k < pivot_idx) {
        return quickselect(arr, left, pivot_idx - 1, k);
    } else {
        return quickselect(arr, pivot_idx + 1, right, k);
    }
}

void select_top_k_linear(ColorFreq* arr, int n, int k, ColorFreq* result) {
    if (n <= k) {

        for (int i = 0; i < n; i++) {
            result[i] = arr[i];
        }
        return;
    }
    
    quickselect(arr, 0, n - 1, k - 1);
    
    for (int i = 0; i < k; i++) {
        result[i] = arr[i];
    }
}

void global_radix_sort_by_frequency() {

    for (int v = 0; v < node_count; v++) {
        int n = k_most_frequent[v].count;
        if (n <= 1) continue;
        
        ColorFreq* arr = k_most_frequent[v].colors;
        
        int max_freq = 0;
        for (int i = 0; i < n; i++) {
            if (arr[i].frequency > max_freq) {
                max_freq = arr[i].frequency;
            }
        }
        
        if (max_freq == 0) continue;
        
        memset(shared_count_array, 0, (max_freq + 1) * sizeof(int));
        
        for (int i = 0; i < n; i++) {
            shared_count_array[arr[i].frequency]++;
        }
        
        int pos = n;
        for (int f = 0; f <= max_freq; f++) {
            int temp = shared_count_array[f];
            shared_count_array[f] = pos - temp;
            pos -= temp;
        }
        
        for (int i = 0; i < n; i++) {
            int freq = arr[i].frequency;
            shared_output_array[shared_count_array[freq]] = arr[i];
            shared_count_array[freq]++;
        }
        
        int start = 0;
        while (start < n) {
            int end = start;
            int current_freq = shared_output_array[start].frequency;
            while (end < n && shared_output_array[end].frequency == current_freq) {
                end++;
            }
            
            if (end - start > 1) {
                qsort(&shared_output_array[start], end - start, sizeof(ColorFreq), 
                      compare_color_freq_by_id);
            }
            
            start = end;
        }
        
        for (int i = 0; i < n; i++) {
            arr[i] = shared_output_array[i];
        }
    }
    
}

void mergeAndSelectTopK(KMostFrequent* result, KMostFrequent* sources[], int source_count) {

    ColorFreq* temp = (ColorFreq*)malloc(Delta * sizeof(ColorFreq));
    int temp_count = 0;
    
    int* freq_map = (int*)calloc(Delta, sizeof(int));
    bool* seen = (bool*)calloc(Delta, sizeof(bool));
    
    for (int s = 0; s < source_count; s++) {
        if (sources[s] == NULL) continue;
        for (int i = 0; i < sources[s]->count; i++) {
            int c = sources[s]->colors[i].string_id;
            int f = sources[s]->colors[i].frequency;
            if (c >= 0 && c < Delta) {
                freq_map[c] += f;
                if (!seen[c]) {
                    seen[c] = true;
                    temp_count++;
                }
            }
        }
    }
    
    temp_count = 0;
    for (int c = 0; c < Delta; c++) {
        if (freq_map[c] > 0) {
            temp[temp_count].string_id = c;
            temp[temp_count].frequency = freq_map[c];
            temp_count++;
        }
    }
    
    qsort(temp, temp_count, sizeof(ColorFreq), compare_color_freq);
    
    result->count = (temp_count < k_value) ? temp_count : k_value;
    for (int i = 0; i < result->count; i++) {
        result->colors[i] = temp[i];
    }
    
    free(temp);
    free(freq_map);
    free(seen);
}

void readInputStrings(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open file %s\n", filename);
        exit(1);
    }
    
    char line[100000];
    Delta = 0;
    
    while (fgets(line, sizeof(line), file) && Delta < MAX_DELTA) {
        int len = strlen(line);
        if (len > 0 && line[len-1] != '\n' && len == sizeof(line) - 1) {
            int ch;
            while ((ch = fgetc(file)) != EOF && ch != '\n') {

            }
        }
        if (len > 0 && line[len-1] == '\n') {
            len--;
            line[len] = '\0';
        }
        if (len > 0) {
            Delta++;
        }
    }
    
    if (Delta == 0) {
        printf("Error: No strings found in input file\n");
        exit(1);
    }
    
    if (Delta > MAX_DELTA) {
        printf("Error: Delta (%d) exceeds MAX_DELTA (%d)\n", Delta, MAX_DELTA);
        exit(1);
    }
    
    strings = (char**)malloc(Delta * sizeof(char*));
    string_lengths = (int*)malloc(Delta * sizeof(int));
    string_ids = (int*)malloc(Delta * sizeof(int));
    
    rewind(file);
    int idx = 0;
    while (fgets(line, sizeof(line), file) && idx < Delta) {
        int len = strlen(line);
        if (len > 0 && line[len-1] != '\n' && len == sizeof(line) - 1) {
            int ch;
            while ((ch = fgetc(file)) != EOF && ch != '\n') {

            }
        }
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
    
    if (idx != Delta) {
        printf("Error: Expected %d strings, but only read %d strings\n", Delta, idx);
        exit(1);
    }
    
    fclose(file);
    
    int actual_total_length = 0;
    for (int i = 0; i < Delta; i++) {
        actual_total_length += string_lengths[i] + 1;
    }
    
    text_length = actual_total_length;
    concatenated_text = (unsigned char*)malloc(text_length * sizeof(unsigned char));
    position_to_string_id = (int*)malloc(text_length * sizeof(int));
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
    
    if (pos != text_length) {
        printf("Error: Mismatch in text_length calculation: expected %d, got %d\n", text_length, pos);
        exit(1);
    }
}

int compare_suffixes_qsort(const void* a, const void* b) {
    int i = *(int*)a;
    int j = *(int*)b;
    
    while (i < text_length && j < text_length) {
        if (concatenated_text[i] != concatenated_text[j]) {
            return concatenated_text[i] - concatenated_text[j];
        }
        i++;
        j++;
    }
    return (i < text_length) - (j < text_length);
}

void buildSuffixArray() {
    SA = (int*)malloc(text_length * sizeof(int));
    ISA = (int*)malloc(text_length * sizeof(int));
    
    #ifdef USE_DIVSUFSORT
        #ifdef _USE_64
            if (divsufsort64((sauchar_t*)concatenated_text, (saidx64_t*)SA, (saidx64_t)text_length) != 0) {
                fprintf(stderr, "Error: SA computation failed (divsufsort64).\n");
                exit(1);
            }
        #else
            if (divsufsort((sauchar_t*)concatenated_text, (saidx_t*)SA, (saidx_t)text_length) != 0) {
                fprintf(stderr, "Error: SA computation failed (divsufsort).\n");
                exit(1);
            }
        #endif
    #else
        for (int i = 0; i < text_length; i++) {
            SA[i] = i;
        }
        qsort(SA, text_length, sizeof(int), compare_suffixes_qsort);
    #endif
    
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

void buildSAPositionToNodeMapping() {

    sa_to_node = (int*)malloc(text_length * sizeof(int));
    
    for (int i = 0; i < text_length; i++) {
        sa_to_node[i] = -1;
    }
    
    for (int i = 0; i < node_count; i++) {
        SuffixTreeNode* node = &suffix_tree_nodes[i];
        if (node->l == node->r) {

            sa_to_node[node->l] = i;
        }
    }
    
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

    buildSAPositionToNodeMapping();
    
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
            k_most_frequent[v].colors[0].frequency = 1;
            k_most_frequent[v].colors[0].string_id = str_id;
            k_most_frequent[v].count = 1;
        } else {
            k_most_frequent[v].count = 0;
        }
        return;
    }
    
    for (int i = 0; i < node->child_count; i++) {
        mergeCounts_DFS(node->children[i]);
    }
    
    memset(shared_freq_map, 0, Delta * sizeof(int));
    
    int active_colors_count = 0;
    
    for (int i = 0; i < node->child_count; i++) {
        int child = node->children[i];
        KMostFrequent* child_k = &k_most_frequent[child];
        for (int j = 0; j < child_k->count; j++) {
            int c = child_k->colors[j].string_id;
            int f = child_k->colors[j].frequency;
            if (c >= 0 && c < Delta) {
                if (shared_freq_map[c] == 0) {

                    shared_active_colors[active_colors_count++] = c;
                }
                shared_freq_map[c] += f;
            }
        }
    }
    
    MappedNodes* mapping = &reverse_mapping[v];
    for (int i = 0; i < mapping->count; i++) {
        int c = mapping->pairs[i].color_id;
        int f = mapping->pairs[i].count;
        if (c >= 0 && c < Delta) {
            if (shared_freq_map[c] == 0) {

                shared_active_colors[active_colors_count++] = c;
            }
            shared_freq_map[c] += f;
        }
    }
    
    int temp_count = 0;
    for (int i = 0; i < active_colors_count; i++) {
        int c = shared_active_colors[i];
        if (shared_freq_map[c] > 0) {
            shared_temp[temp_count].string_id = c;
            shared_temp[temp_count].frequency = shared_freq_map[c];
            temp_count++;
        }
    }
    
    if (temp_count == 0) {
        k_most_frequent[v].count = 0;
    } else {
        int k_to_select = (temp_count < k_value) ? temp_count : k_value;
        select_top_k_linear(shared_temp, temp_count, k_to_select, k_most_frequent[v].colors);
        k_most_frequent[v].count = k_to_select;
    }
}

void step3_mergeCounts() {

    buildReverseMapping();
    
    shared_freq_map = (int*)calloc(Delta, sizeof(int));
    shared_temp = (ColorFreq*)malloc(Delta * sizeof(ColorFreq));
    shared_active_colors = (int*)malloc(Delta * sizeof(int));
    
    int max_count_size = text_length + 1;
    shared_count_array = (int*)calloc(max_count_size, sizeof(int));
    shared_output_array = (ColorFreq*)malloc(k_value * sizeof(ColorFreq));
    
    k_most_frequent = (KMostFrequent*)malloc(node_count * sizeof(KMostFrequent));
    for (int i = 0; i < node_count; i++) {
        k_most_frequent[i].colors = (ColorFreq*)malloc(k_value * sizeof(ColorFreq));
        k_most_frequent[i].count = 0;
    }
    
    mergeCounts_DFS(0);
    
    global_radix_sort_by_frequency();
    
    free(shared_freq_map);
    shared_freq_map = NULL;
    free(shared_temp);
    shared_temp = NULL;
    free(shared_active_colors);
    shared_active_colors = NULL;
    free(shared_count_array);
    shared_count_array = NULL;
    free(shared_output_array);
    shared_output_array = NULL;
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
    if (k_most_frequent) {
        for (int i = 0; i < node_count; i++) {
            if (k_most_frequent[i].colors) {
                free(k_most_frequent[i].colors);
            }
        }
        free(k_most_frequent);
        k_most_frequent = NULL;
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
    if (argc != 4) {
        printf("Usage: %s <input_file> <output_file> <k>\n", argv[0]);
        printf("\nInput format: Each line is one string\n");
        printf("\nParameters:\n");
        printf("  k: Number of most frequent colors to compute\n");
        printf("     Must be between 1 and %d\n", MAX_K);
        printf("\nExample:\n");
        printf("  %s input.txt output.txt 3\n", argv[0]);
        return 1;
    }
    
    const char* input_file = argv[1];
    const char* output_file = argv[2];
    
    k_value = atoi(argv[3]);
    if (k_value < 1 || k_value > MAX_K) {
        printf("Error: k must be between 1 and %d\n", MAX_K);
        return 1;
    }
    
    readInputStrings(input_file);
    buildSuffixArray();
    buildLCPArray();
    buildSuffixTreeFromLCP();
    
    step1_phaseA();
    step1_phaseB();
    step2_countColors();
    step3_mergeCounts();
    
    KMostFrequent* root_k = &k_most_frequent[0];
    
    FILE* out = fopen(output_file, "w");
    if (out) {
        fprintf(out, "k = %d\n", k_value);
        fprintf(out, "Root Node - Top %d Most Frequent Colors:\n", root_k->count);
        for (int i = 0; i < root_k->count; i++) {
            fprintf(out, "%d %d\n", root_k->colors[i].string_id, root_k->colors[i].frequency);
        }
        fclose(out);
    }
    
    cleanup();
    return 0;
}
