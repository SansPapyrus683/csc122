#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "hashing.h"

using namespace std;

// const string PROJ = "c";
// const vector<int> USE{24, 49, 64, 94, 1, 37, 67, 42, 78, 14};
const string PROJ = "d";
const vector<int> USE{965,  3767, 2664, 1972, 2307, 1137, 2052, 388,  1172, 4258,
                      2463, 2413, 1945, 4119, 4363, 4325, 2103, 2911, 593,  2373,
                      4689, 1372, 2367, 2881, 2251, 2016, 2932, 3144, 4353, 971,
                      4556, 2900, 1624, 1792, 4626, 209,  777,  2083, 341,  969,
                      2975, 4442, 4080, 4708, 4892, 4480, 683,  119,  1141};

constexpr int KMER_LEN = 15;
constexpr int HASH_FUNCS = 20;
constexpr long long MOD = 1e9 + 7;

int main() {
    hashing_init(KMER_LEN, HASH_FUNCS, MOD);

    map<int, bitset<MOD>> filters;

    string line;
    for (int g : USE) {
        ostringstream raw_path;
        raw_path << "1" << PROJ << "/project1" << PROJ << "_genome_" << g << ".fasta";
        ifstream fin(raw_path.str());
        if (!fin) { break; }

        string genome;
        getline(fin, line);
        while (getline(fin, line)) { genome += line; }
        fin.close();

        filters[g] = bitset<MOD>();
        bitset<MOD>& f = filters[g];

        for (const vector<long long>& hash_set : kmer_hashes(genome)) {
            for (long long hsh : hash_set) { f.set(hsh); }
        }
    }

    ostringstream read_path;
    read_path << "1" << PROJ << "/project1" << PROJ << "_reads.fasta";
    ifstream fin(read_path.str());
    int v = 0;
    while (getline(fin, line)) {
        v++;
        if (v % 2 == 1 || line.size() < KMER_LEN) { continue; }

        vector<vector<long long>> hashes(kmer_hashes(line));

        pair<int, int> best{0, 0};
        for (const auto& [i, f] : filters) {
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
            best = max(best, {amt, i});
        }
        printf(">read_%i\tGenome_Number_%i\n", v / 2 - 1, best.second);
    }
    fin.close();
}
