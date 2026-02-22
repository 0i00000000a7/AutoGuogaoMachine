#include <fstream>
#include "PPSNotation.h"
#include <omp.h>
using namespace std;
const int ENUM_LAYERS = 65;
const int MAX_NONLOOP_LEN = ENUM_LAYERS / 4;
const int CHECK_COPY_NUM=5;
const int PRINT_THRESHOLD = ENUM_LAYERS - 13;

PPSNotation smallest_nonterminate_seq;
const vector<INT> known_smallest_seq =
{
    0,1,0,2,2,0,5,2,0,8,8,6,0,12,8,6,0,16,16,15,0,16,16,15,0,15,0,26,26,15,0,26,26,13,6,2,0,26,13,8,6,0,0,15,0,44,44,15,0,44,15,0,51,51,0,54,54,49,15,0,59,59,0,62,59,0,65,65,63,0,69,65,63,0,73,73,72,0,73,73,72,0,72,0,83,83,72,0,83,83,70,63,59,0,83,70,65,63,0,0,72,0,101,101,72,0,101,72,0,108,108,0,111,111,106,72,0,116,116,0,119,116,0,122,122,120,0,126,122,120,0,130,130,129,0,130,130,129,0,129,0,140,140,129,0,140,140,127,120,116,0,140,127,122,120,0,0,129,0,158,158,129,0,158,129,0,165,165,0,168,168,163,129,0,173,173,0,176,173,0,179,179,177,0,183,179,177,0,187,187,186,0,187,187,186,0,186,0,197,197,186,0,197,197,184,177,173,0,197,184,179,177,0,0,186,0,215,215,186,0,215,186,0,222,222,0,225,225,220,186,0,230,230,0,233,230,0,236,236,234,0,240,236,234,0,244,244,243,0,244,244,243,0,243,0,254,254,243,0,254,254,241,234,230,0,254,241,236,234,0,0,243,0,272,272,243,0,272,243,0,279,279,0,282,282,277,243,0,287,287,0,290,287,0,293,293,291,0,297,293,291,0,301,301,300,0,301,301,300,0,300,0,311,311,300,0,311,311,298,291,287,0,311,298,293,291,0,0,300,0,329,329,300,0,329,300,0,336,336,0,339,339,334,300
};
const int PARALLEL_DEPTH = ENUM_LAYERS - 10;
const int BASE_PREFIX_LEN = ENUM_LAYERS - 20;

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
