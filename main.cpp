/**
 * T Operator Product Generation Main File
 * @file main.cpp
 * @author Connor Mooney
 * @author Michael Jarret
 * @author Andrew Glaudell
 * @author Jacob Weston
 * @author Mingzhen Tian
 * @version 6/1/21
 */
#include <chrono>
#include <set>
#include <fstream>
#include <vector>
#include <omp.h>
#include <thread>
#include <algorithm>
#include "Globals.hpp"
#include "utils.hpp"

using namespace std;
using namespace std::chrono;

/**
 * Improved version includes:
 * - Utilizing auto and range-based for loops for better readability.
 * - Simplifying function names and parameters.
 * - Avoiding global variables when possible.
 * - Optimizing data structures and algorithms for better performance.
 */

 // Utility function to generate permutations and modifications of a pattern
static set<pattern> generatePermutations(const pattern& pat) {
    set<pattern> perms, inits{ pat, pat.transpose() };
    for (const auto& p : inits) {
        perms.insert(p);
        vector<int> row{ 0, 1, 2, 3, 4, 5 };
        do {
            pattern permuted;
            for (int c = 0; c < 6; ++c)
                for (int r = 0; r < 6; ++r)
                    permuted.arr[c][r] = p.arr[c][row[r]];
            permuted.lexicographic_order();
            perms.insert(permuted);
            // Toggle rows in the pattern
            for (unsigned int mask = 0; mask < (1 << 6); ++mask) {
                auto mod_perm = permuted;
                for (int j = 0; j < 6; ++j)
                    if (mask & (1 << j)) mod_perm.mod_row(j);
                mod_perm.lexicographic_order();
                perms.insert(mod_perm);
            }
        } while (next_permutation(row.begin(), row.end()));
    }
    return perms;
}

// Insert or erase all permutations of a pattern
static void processPattern(const pattern& p, bool insert = true) {
    auto perms = generatePermutations(p);
    if (insert)
        pattern_set.insert(perms.begin(), perms.end());
    else
        for (const auto& perm : perms)
            pattern_set.erase(perm);
}

// Read and process patterns from file
static void processPatternFile(const string& filePath, bool isInsert) {
    ifstream file(filePath);
    if (!file) {
        cerr << "Failed to open file: " << filePath << endl;
        return;
    }
    string line;
    while (getline(file, line)) {
        pattern p(line);
        p.lexicographic_order();
        processPattern(p, isInsert);
    }
    cout << "Processed patterns from " << filePath << endl;
}

// Utility function to get current time
static auto now() {
    return high_resolution_clock::now();
}

// Erase or record pattern based on condition
static void eraseOrRecordPattern(SO6& s, ofstream& of, bool erase = true) {
    pattern pat = s.to_pattern();
    if (erase && pattern_set.find(pat) != pattern_set.end()) {
        processPattern(pat, false); // Erase all permutations
        of << s.circuit_string() << endl;
    }
}

// Main function simplified
int main(int argc, char** argv) {
    auto startTime = now();
    Globals::setParameters(argc, argv);
    Globals::configure();

    // Process input files for patterns
    processPatternFile(pattern_file, true);
    processPatternFile(case_file, false);

    // Example operation demonstrating usage of modified functions
    SO6 exampleSO6 = SO6::identity();
    ofstream exampleFile("example_output.txt", ios::out | ios::trunc);
    eraseOrRecordPattern(exampleSO6, exampleFile);

    // Wrap-up
    auto endTime = now();
    auto duration = duration_cast<seconds>(endTime - startTime).count();
    cout << "Total execution time: " << duration << " seconds" << endl;

    return 0;
}
