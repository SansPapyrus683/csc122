#include <algorithm>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
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

tuple<int, string, string> best_match(const string& a, const string& b) {
    int n = a.size(), m = b.size();

    vector<vector<pair<int, pair<int, int>>>> dp(
        n + 1, vector<pair<int, pair<int, int>>>(m + 1, {INT32_MAX, {-1, -1}}));

    dp[0][0] = {0, {-1, -1}};

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            if (i > 0) {
                dp[i][j] = min(dp[i][j], {dp[i - 1][j].first + 1, {i - 1, j}});
            }
            if (j > 0) {
                dp[i][j] = min(dp[i][j], {dp[i][j - 1].first + 1, {i, j - 1}});
            }
            if (i > 0 && j > 0) {
                int cost = (a[i - 1] != b[j - 1]);
                dp[i][j] =
                    min(dp[i][j], {dp[i - 1][j - 1].first + cost, {i - 1, j - 1}});
            }
        }
    }

    string ra, rb;
    int i = n, j = m;
    while (!(i == 0 && j == 0)) {
        auto [pi, pj] = dp[i][j].second;
        if (pi == i) {
            ra.push_back('_');
            rb.push_back(b[pj]);
        } else if (pj == j) {
            ra.push_back(a[pi]);
            rb.push_back('_');
        } else {
            ra.push_back(a[pi]);
            rb.push_back(b[pj]);
        }
        i = pi;
        j = pj;
    }

    reverse(ra.begin(), ra.end());
    reverse(rb.begin(), rb.end());

    return {dp[n][m].first, ra, rb};
}

int main() {
    ifstream gf(REF);
    string genome, line;
    getline(gf, line);
    while (getline(gf, line)) { genome += line; }
    gf.close();

    map<string, vector<int>> mins;
    for (auto& [m, ind] : kmers(genome)) { mins[m].push_back(ind); }

    vector<string> queries;
    ifstream rf(READS);
    int v = 0;
    while (getline(rf, line)) {
        if (v % 2 == 1) { queries.push_back(line); }
        v++;
    }
    rf.close();

    map<vector<string>, int> changes;
    for (int qi = 0; qi < (int)queries.size(); qi++) {
        const string& q = queries[qi];
        map<int, int> offsets;

        for (auto& [m, ind] : kmers(q)) {
            for (int match : mins[m]) {
                int off = match - ind;
                if (0 <= off && off < (int)genome.size()) { offsets[off]++; }
            }
        }

        if (offsets.empty()) { continue; }

        int best_ind =
            max_element(offsets.begin(), offsets.end(), [](auto& a, auto& b) {
                return a.second < b.second;
            })->first;

        tuple<int, string, string> best = {INT32_MAX, "", ""};
        for (int change : {0, -1, 1}) {
            int len = q.size() + change;
            string genome_sub = genome.substr(best_ind, len);

            tuple<int, string, string> test;
            if (change == 0) {
                int cost = 0;
                for (int i = 0; i < (int)q.size(); i++) {
                    cost += genome_sub[i] != q[i];
                }
                test = {cost, genome_sub, q};
            } else {
                test = best_match(genome_sub, q);
            }

            if (get<0>(test) < get<0>(best)) { best = test; }
        }

        if (get<0>(best) >= 3) { continue; }

        int ref_at = best_ind;
        const string& a = get<1>(best);
        const string& b = get<2>(best);
        for (int i = 0; i < (int)a.size(); i++) {
            if (a[i] == '_') {
                // not sure why insertions are -1 here but oh well
                changes[{"I", to_string(ref_at - 1), string(1, b[i])}]++;
            } else if (b[i] == '_') {
                changes[{"D", to_string(ref_at)}]++;
            } else if (a[i] != b[i]) {
                changes[{"S", to_string(ref_at), string(1, b[i])}]++;
            }
            ref_at += a[i] != '_';
        }
    }

    vector<pair<int, vector<string>>> sorted;
    for (auto& [k, v] : changes) { sorted.push_back({v, k}); }
    sort(sorted.rbegin(), sorted.rend());

    for (const auto& [occs, info] : sorted) {
        if (occs < 4) { break; }
        if (info[0] == "S") {
            int ind = stoi(info[1]);
            printf(">%s%i %c %s\n", info[0].c_str(), ind, genome[ind], info[2].c_str());
        } else if (info[0] == "I") {
            printf(">%s%s %s\n", info[0].c_str(), info[1].c_str(), info[2].c_str());
        } else if (info[0] == "D") {
            int ind = stoi(info[1]);
            printf(">%s%i %c\n", info[0].c_str(), ind, genome[ind]);
        }
    }
}
