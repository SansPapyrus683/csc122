#include <algorithm>
#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include "hashing.h"

using namespace std;

const string PROJ = "c_sample";

constexpr int KMER_LEN = 30;
constexpr int HASH_FUNCS = 10;
constexpr long long MOD = 2875519;

int main() {
    hashing_init(KMER_LEN, HASH_FUNCS, MOD);

    vector<bitset<MOD>> filters;

    string line;
    for (int g = 0;; g++) {
        ostringstream raw_path;
        raw_path << "1" << PROJ << "/project1" << PROJ << "_genome_" << g << ".fasta";
        ifstream fin(raw_path.str());
        if (!fin) { break; }

        string genome;
        getline(fin, line);
        while (getline(fin, line)) { genome += line; }
        fin.close();

        filters.push_back(bitset<MOD>());
        bitset<MOD>& f = filters.back();

        for (const vector<long long>& hash_set : kmer_hashes(genome)) {
            for (long long hsh : hash_set) { f.set(hsh); }
        }
        cerr << g << endl;
    }

    vector<pair<int, int>> freq(filters.size());
    for (int i = 0; i < filters.size(); i++) { freq[i].second = i; }

    ostringstream read_path;
    read_path << "1" << PROJ << "/project1" << PROJ << "_reads.fasta";
    ifstream fin(read_path.str());
    int v = 0;
    while (getline(fin, line)) {
        v++;
        if (v % 2 == 1 || line.size() < KMER_LEN) { continue; }

        vector<vector<long long>> hashes(kmer_hashes(line));

        pair<int, vector<int>> best{-1, {}};
        for (int i = 0; i < filters.size(); i++) {
            int amt = 0;
            for (const vector<long long>& mer : hashes) {
                bool good = true;
                for (long long h : mer) {
                    if (!filters[i][h]) {
                        good = false;
                        break;
                    }
                }
                amt += good;
            }

            if (amt > best.first) {
                best = {amt, {i}};
            } else if (amt == best.first) {
                best.second.push_back(i);
            }
        }

        for (int g : best.second) {
            freq[g].first--;  // too lazy to sort in reverse
        }
    }
    fin.close();

    sort(freq.begin(), freq.end());
    for (const auto& [amt, g] : freq) {
        if (amt == 0) { break; }
        cout << -amt << ' ' << g << endl;
    }
}
