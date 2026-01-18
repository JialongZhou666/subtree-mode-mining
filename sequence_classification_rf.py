#!/usr/bin/env python3
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import (confusion_matrix, accuracy_score, 
                            precision_recall_fscore_support, roc_auc_score)
from collections import Counter
from itertools import product

K = 7
N_ESTIMATORS = 100
CLASS_WEIGHT = 'balanced'
RANDOM_STATE = 42

def generate_kmers(sequence, k):
    sequence = sequence.upper().strip()
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if all(c in 'ATCGUN' for c in kmer):
            kmers.append(kmer)
    return kmers

def get_all_possible_kmers(k, alphabet='ATCG'):
    return [''.join(p) for p in product(alphabet, repeat=k)]

def sequence_to_kmer_vector(sequence, k, kmer_list):
    kmers = generate_kmers(sequence, k)
    kmer_counts = Counter(kmers)
    total = len(kmers) if len(kmers) > 0 else 1
    vector = np.array([kmer_counts.get(kmer, 0) / total for kmer in kmer_list])
    return vector

def extract_features_from_sequences(sequences, k):
    kmer_list = get_all_possible_kmers(k)
    features = []
    for seq in sequences:
        vector = sequence_to_kmer_vector(seq, k, kmer_list)
        features.append(vector)
    return np.array(features), kmer_list

def load_sequences_from_file(filepath):
    sequences = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('>'):
                line = line.replace(' ', 'N').replace('\t', 'N')
                sequences.append(line)
    return sequences

def load_data():
    sars_sequences = load_sequences_from_file("Sars-Cov-2.txt")
    infl_train_sequences = load_sequences_from_file("influenza_genomes_1000.txt")
    hum_train_sequences = load_sequences_from_file("extracted_human_sequences_1000.txt")
    infl_test_sequences = load_sequences_from_file("Influenza.txt")
    hum_test_sequences = load_sequences_from_file("Human.txt")
    
    train_sars = sars_sequences[:1950]
    train_infl = infl_train_sequences[:1000]
    train_hum = hum_train_sequences[:1000]
    
    train_sequences = train_sars + train_infl + train_hum
    train_labels = [1] * len(train_sars) + [0] * (len(train_infl) + len(train_hum))
    
    test_sars = sars_sequences[1950:2000]
    test_infl = infl_test_sequences[:50]
    test_hum = hum_test_sequences[:50]
    
    test_sequences = test_sars + test_infl + test_hum
    test_labels = [1] * len(test_sars) + [0] * (len(test_infl) + len(test_hum))
    
    return train_sequences, train_labels, test_sequences, test_labels

def train_random_forest(X_train, y_train):
    clf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        max_depth=None,
        criterion='gini',
        class_weight=CLASS_WEIGHT,
        min_samples_split=2,
        min_samples_leaf=1,
        random_state=RANDOM_STATE,
        n_jobs=1,
        verbose=0
    )
    clf.fit(X_train, y_train)
    return clf

def evaluate_model(clf, X_train, y_train, X_test, y_test):
    y_train_pred = clf.predict(X_train)
    y_test_pred = clf.predict(X_test)
    
    y_train_proba = clf.predict_proba(X_train)[:, 1]
    y_test_proba = clf.predict_proba(X_test)[:, 1]
    
    results = {}
    results['train_accuracy'] = accuracy_score(y_train, y_train_pred)
    results['train_auc'] = roc_auc_score(y_train, y_train_proba)
    results['test_accuracy'] = accuracy_score(y_test, y_test_pred)
    results['test_auc'] = roc_auc_score(y_test, y_test_proba)
    
    precision, recall, f1, _ = precision_recall_fscore_support(
        y_test, y_test_pred, average='binary'
    )
    results['test_precision'] = precision
    results['test_recall'] = recall
    results['test_f1'] = f1
    
    cv_scores = cross_val_score(clf, X_train, y_train, cv=10, scoring='accuracy')
    results['cv_accuracy'] = cv_scores.mean()
    results['cv_std'] = cv_scores.std()
    
    cm = confusion_matrix(y_test, y_test_pred)
    results['confusion_matrix'] = cm
    
    return results

def main():
    train_sequences, train_labels, test_sequences, test_labels = load_data()
    
    X_train, kmer_list = extract_features_from_sequences(train_sequences, K)
    y_train = np.array(train_labels)
    
    X_test_list = []
    for seq in test_sequences:
        vector = sequence_to_kmer_vector(seq, K, kmer_list)
        X_test_list.append(vector)
    X_test = np.array(X_test_list)
    y_test = np.array(test_labels)
    
    clf = train_random_forest(X_train, y_train)
    results = evaluate_model(clf, X_train, y_train, X_test, y_test)
    
    print("K =", K)
    print("Train Accuracy: {:.4f}".format(results['train_accuracy']))
    print("Train AUC: {:.4f}".format(results['train_auc']))
    print("Test Accuracy: {:.4f}".format(results['test_accuracy']))
    print("Test AUC: {:.4f}".format(results['test_auc']))
    print("Test Precision: {:.4f}".format(results['test_precision']))
    print("Test Recall: {:.4f}".format(results['test_recall']))
    print("Test F1: {:.4f}".format(results['test_f1']))
    print("CV Accuracy: {:.4f} (+/- {:.4f})".format(results['cv_accuracy'], results['cv_std']))
    cm = results['confusion_matrix']
    print("Confusion Matrix: TN={}, FP={}, FN={}, TP={}".format(cm[0,0], cm[0,1], cm[1,0], cm[1,1]))

if __name__ == "__main__":
    main()
