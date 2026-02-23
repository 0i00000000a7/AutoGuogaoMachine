#include <fstream>
#include "PPSNotation.h"
#include <omp.h>
using namespace std;
const int ENUM_LAYERS = 78;
const int MAX_NONLOOP_LEN = ENUM_LAYERS / 4;
const int CHECK_COPY_NUM=5;
const int PRINT_THRESHOLD = ENUM_LAYERS - 13;

PPSNotation smallest_nonterminate_seq;
const vector<INT> known_smallest_seq =
{
0,1,0,2,2,0,5,2,0,8,8,6,0,12,8,6,0,16,16,15,0,16,16,15,0,15,0,26,26,15,0,26,26,13,6,2,0,15,0,38,38,15,0,38,15,0,45,45,13,8,6,0,0,38,15,0,55,55,13,0,59,57,13,0,0,59,53,15,0,0,13,0,71,55,0,74,74,74,72,55,0,80,80,0,83,80,0,86,86,84,0,90,86,84,0,94,94,93,0,94,94,93,0,93,0,104,104,93,0,104,104,91,84,80,0,93,0,116,116,93,0,116,93,0,123,123,91,86,84,0,0,116,93,0,133,133,91,0,137,135,91,0,0,137,131,93,0,0,91,0,149,133,0,152,152,152,150,133,0,158,158,0,161,158,0,164,164,162,0,168,164,162,0,172,172,171,0,172,172,171,0,171,0,182,182,171,0,182,182,169,162,158,0,171,0,194,194,171,0,194,171,0,201,201,169,164,162,0,0,194,171,0,211,211,169,0,215,213,169,0,0,215,209,171,0,0,169,0,227,211,0,230,230,230,228,211,0,236,236,0,239,236,0,242,242,240,0,246,242,240,0,250,250,249,0,250,250,249,0,249,0,260,260,249,0,260,260,247,240,236,0,249,0,272,272,249,0,272,249,0,279,279,247,242,240,0,0,272,249,0,289,289,247,0,293,291,247,0,0,293,287,249,0,0,247,0,305,289,0,308,308,308,306,289,0,314,314,0,317,314,0,320,320,318,0,324,320,318,0,328,328,327,0,328,328,327,0,327,0,338,338,327,0,338,338,325,318,314,0,327,0,350,350,327,0,350,327,0,357,357,325,320,318,0,0,350,327,0,367,367,325,0,371,369,325,0,0,371,365,327,0,0,325,0,383,367,0,386,386,386,384,367,0,392,392,0,395,392,0,398,398,396,0,402,398,396,0,406,406,405,0,406,406,405,0,405,0,416,416,405,0,416,416,403,396,392,0,405,0,428,428,405,0,428,405,0,435,435,403,398,396,0,0,428,405,0,445,445,403,0,449,447,403,0,0,449,443,405,0,0,403,0,461,445,0,464,464,464,462,445
};
const int PARALLEL_DEPTH = 70;
const int BASE_PREFIX_LEN = 60;

void check_nonterminate_seq(const vector<INT>& seq, const vector<INT>& loop_dif, PPSNotation& local_best);
void check_0nn(const vector<INT>& seq, PPSNotation& local_best);
void enum_pps_recursive(vector<INT>& prev, PPSNotation& local_best);

void check_nonterminate_seq(const vector<INT>& seq, const vector<INT>& loop_dif, PPSNotation& local_best)
{
    INT loop_len = loop_dif.size();
    PPSNotation copied_seq;
    copied_seq.seq = seq;
    INT loop_start = seq.size() - loop_len;
    for (INT i = 0; i < CHECK_COPY_NUM; i++) {
        for (INT j = 0; j < loop_len; j++) {
            copied_seq.seq.push_back(seq[loop_start + j] + (i + 1) * loop_dif[j]);
        }
    }
    assert(copied_seq.seq.size() == seq.size() + CHECK_COPY_NUM * loop_len);
    vector<bool> res;
    bool is_standard = copied_seq._checkStandardAndNonMaximum(res);
    if (!is_standard) return;

    INT last_loop_start = copied_seq.seq.size() - 2 * loop_len;
    bool is_same = true;
    bool is_empty = true;
    for (INT i = 0; i < loop_len; i++) {
        if (res[last_loop_start + i] != res[last_loop_start + loop_len + i]) {
            is_same = false;
            break;
        }
        if (res[last_loop_start + i]) {
            is_empty = false;
        }
    }
    if (is_same && !is_empty) {
        #pragma omp critical
        {
            if (smallest_nonterminate_seq.seq.empty() || smallest_nonterminate_seq > copied_seq) {
                smallest_nonterminate_seq = copied_seq;
                ofstream ofs("pps_infinite.txt", ios::app);
                cout << "Found smaller nonterminate seq: ";
                copied_seq.print(cout);
                copied_seq.print(ofs);
                cout << endl;
                ofs << endl;
            }
        }

        if (local_best.seq.empty() || local_best > copied_seq) {
            local_best = copied_seq;
        }
    }
}

void check_loop_dif(const vector<INT>& seq, PPSNotation& local_best)
{
    INT min_loop_len = max((INT)1, (INT)(seq.size() - MAX_NONLOOP_LEN));
    for (int loop_len = min_loop_len; loop_len <= (int)seq.size(); loop_len++) {
        vector<INT> loop_dif;
        for (int j = 0; j < loop_len; j++) {
            loop_dif.push_back(seq[seq.size() - loop_len + j] > 0 ? loop_len : 0);
        }
        check_nonterminate_seq(seq, loop_dif, local_best);
    }
}

void enum_pps_recursive(vector<INT>& prev, PPSNotation& local_best)
{
    PPSNotation pps;
    pps.seq = prev;

    if (!local_best.seq.empty() && local_best <= pps)
        return;

    vector<bool> res;
    bool is_standard = pps._checkStandardAndNonMaximum(res);
    if (!is_standard) return;

    if (!pps._isSuccessor()) {
        check_loop_dif(pps.seq, local_best);
    }

    if (prev.size() <= PRINT_THRESHOLD) {
        #pragma omp critical
        {
            pps.print(cout);
            cout << endl;
        }
    }

    if (prev.size() < ENUM_LAYERS) {
        INT max_i = prev.size();
        if (prev.size() < PARALLEL_DEPTH) {
            for (int i = max_i; i >= 0; i--) {
                PPSNotation local_best_copy = local_best;
                #pragma omp task firstprivate(local_best_copy)
                {
                    vector<INT> new_prev = prev;
                    new_prev.push_back(i);
                    enum_pps_recursive(new_prev, local_best_copy);
                }
            }
        } else {
            for (int i = max_i; i >= 0; i--) {
                prev.push_back(i);
                enum_pps_recursive(prev, local_best);
                prev.pop_back();
            }
        }
    }
}

int main_enumpps1()
{
    smallest_nonterminate_seq.seq = known_smallest_seq;

    vector<INT> current_prefix(known_smallest_seq.begin(),
                               known_smallest_seq.begin() + BASE_PREFIX_LEN);

    while (true) {
        #pragma omp parallel
        {
            #pragma omp single
            {
                vector<INT> prev = current_prefix;
                PPSNotation local_best;
                #pragma omp critical
                {
                    local_best = smallest_nonterminate_seq;
                }
                enum_pps_recursive(prev, local_best);
            }
        }
        
        if (current_prefix.empty()) return 0;

        PPSNotation tmp;
        tmp.seq = current_prefix;
        tmp._reduce();
        current_prefix = tmp.seq;
    }
    return 0;
}
