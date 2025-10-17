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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }
    
    string inputFile = argv[1];
    string outputFile = argv[2];
    
    ifstream inFile(inputFile);
    if (!inFile.is_open()) {
        cerr << "Error: Cannot open file " << inputFile << endl;
        return 1;
    }
    
    vector<string> sequences;
    string line;
    
    while (getline(inFile, line)) {
        if (line.length() > 0) {
            sequences.push_back(line);
        }
    }
    inFile.close();
    
    MultiColorESA mcesa(sequences);
    mcesa.runSCM();
    
    ofstream outFile(outputFile);
    if (outFile.is_open()) {
        outFile << "Algorithm completed successfully." << endl;
        outFile.close();
    }
    
    return 0;
}