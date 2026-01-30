/**
 * this is a version i tried to make more performant but its like
 * REALLY awful in terms of actual correctness so unlucky ig
 */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

const string REF = "1a/project1a_reference_genome.fasta";
const string READS = "1a/project1a_with_error_paired_reads.fasta";

// const string REF = "1b/project1b-u_reference_genome.fasta";
// const string READS = "1b/project1b-u_with_error_paired_reads.fasta";

vector<pair<string, int>> kmers(const string& s, int k = 20) {
    vector<pair<string, int>> ret;
    for (int i = 0; i + k <= (int)s.size(); i++) {
        ret.emplace_back(s.substr(i, k), i);
    }
    return ret;
}

// you can delete one char from a to make its hamming distance closer to b
pair<int, int> best_del(const string& a, const string& b) {
    assert(b.size() + 1 == a.size());

    vector<int> pref(b.size() + 1);
    vector<int> suff(b.size() + 1);
    for (int i = 0; i < b.size(); i++) {
        pref[i + 1] = pref[i] + (b[i] != a[i]);
        suff[i + 1] = suff[i] + (b[b.size() - i - 1] != a[a.size() - i - 1]);
    }

    pair<int, int> best{a.size(), -1};
    for (int i = 0; i < a.size(); i++) {
        best = min(best, {pref[i] + suff[a.size() - i - 1], i});
    }
    return best;
}

int main() {
    ifstream gf(REF);
    string genome, line;
    getline(gf, line);
    while (getline(gf, line)) { genome += line; }
    gf.close();

    unordered_map<string, vector<int>> mins;
    for (auto& [m, ind] : kmers(genome)) { mins[m].push_back(ind); }

    vector<string> queries;
    ifstream rf(READS);
    int v = 0;
    long long coverage = 0;
    while (getline(rf, line)) {
        if (v % 2 == 1) {
            queries.push_back(line);
            coverage += line.size();
        }
        v++;
    }
    rf.close();

    coverage /= genome.size();

    auto val_cmp = [](const pair<int, int>& a, const pair<int, int>& b) {
        return a.second < b.second;
    };
    map<pair<int, char>, int> subs;
    map<int, int> dels;
    map<pair<int, char>, int> ins;
    for (const string& q : queries) {
        map<int, int> offsets;
        for (auto& [m, ind] : kmers(q)) {
            for (int match : mins[m]) {
                int off = match - ind;
                if (0 <= off && off < (int)genome.size()) { offsets[off]++; }
            }
        }
        if (offsets.empty()) { continue; }

        int best_ind = max_element(offsets.begin(), offsets.end(), val_cmp)->first;

        for (int i = 0; i < q.size(); i++) {
            if (q[i] != genome[best_ind + i]) { subs[{best_ind + i, q[i]}]++; }
        }

        string less = genome.substr(best_ind, q.size() - 1);
        auto [err, ind] = best_del(q, less);
        if (err < q.size() / 3) {}
        ins[{best_ind + ind - 1, q[ind]}]++;

        if (best_ind + q.size() >= genome.size()) { continue; }
        string more = genome.substr(best_ind, q.size() + 1);
        tie(err, ind) = best_del(more, q);
        if (err < q.size() / 3) {}
        dels[best_ind + ind]++;
    }

    vector<pair<int, string>> sorted;
    for (const auto& [s, amt] : subs) {
        ostringstream out;
        out << 'S' << s.first << ' ' << genome[s.first] << ' ' << s.second;
        sorted.push_back({amt, out.str()});
    }
    for (const auto& [s, amt] : ins) {
        ostringstream out;
        out << 'I' << s.first << ' ' << s.second;
        sorted.push_back({amt, out.str()});
    }
    for (const auto& [ind, amt] : dels) {
        ostringstream out;
        out << 'D' << ind << ' ' << genome[ind];
        sorted.push_back({amt, out.str()});
    }
    sort(sorted.rbegin(), sorted.rend());

    for (const auto& [amt, s] : sorted) {
        if (amt < coverage * .65) { break; }
        cout << '>' << s << '\n';
    }
}
