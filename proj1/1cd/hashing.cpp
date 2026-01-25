#include "hashing.h"

#include <vector>
#include <string>

using namespace std;

const std::map<char, int> CHARS{{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};

int kmer_len;
long long mod;

vector<long long> bases;
vector<long long> pows;

// is this best practice? hell if i know.
void hashing_init(int kmer_len_, int hash_funcs, long long mod_) {
    kmer_len = kmer_len_;
    mod = mod_;

    bases.clear();
    for (int i = 0; i < hash_funcs; i++) {
        bases.push_back(4 + i);
    }

    pows = vector<long long>(bases.size(), 1);
    for (int i = 0; i < kmer_len; i++) {
        for (int j = 0; j < pows.size(); j++) { pows[j] = pows[j] * bases[j] % mod; }
    }
}

vector<vector<long long>> kmer_hashes(const string& genome) {
    vector<vector<long long>> prefix(genome.size() + 1,
                                     vector<long long>(bases.size()));
    vector<vector<long long>> ret(genome.size() - kmer_len + 1,
                                  vector<long long>(bases.size()));
    for (int i = 0; i < genome.size(); i++) {
        int val = CHARS.at(genome[i]);
        vector<long long>& curr = prefix[i + 1];
        for (int j = 0; j < bases.size(); j++) {
            curr[j] = (prefix[i][j] * bases[j] % mod + val) % mod;
        }
        if (i >= kmer_len - 1) {
            const vector<long long>& prev = prefix[i - kmer_len + 1];
            for (int j = 0; j < bases.size(); j++) {
                long long hash = curr[j] - prev[j] * pows[j] % mod;
                ret[i - kmer_len + 1][j] = (hash + mod) % mod;
            }
        }
    }
    return ret;
}
