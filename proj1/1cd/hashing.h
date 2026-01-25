#ifndef HASHING_UTILS
#define HASHING_UTILS 1

#include <map>
#include <string>
#include <vector>

void hashing_init(int kmer_len, int hash_funcs, long long mod);

std::vector<std::vector<long long>> kmer_hashes(const std::string& genome);

#endif
