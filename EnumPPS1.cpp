#include <fstream>
#include "PPSNotation.h"
#include <omp.h>
using namespace std;
const int ENUM_LAYERS = 75;
const int MAX_NONLOOP_LEN = ENUM_LAYERS / 4;
const int CHECK_COPY_NUM=5;
const int PRINT_THRESHOLD = ENUM_LAYERS - 13;

PPSNotation smallest_nonterminate_seq;
const vector<INT> known_smallest_seq =
{
0,1,0,2,2,0,5,2,0,8,8,6,0,12,8,6,0,16,16,15,0,16,16,15,0,15,0,26,26,15,0,26,26,13,6,2,0,15,0,38,38,15,0,38,15,0,45,45,13,8,6,0,0,38,15,0,55,55,13,0,59,59,53,15,0,59,59,53,0,68,68,0,71,68,0,74,74,72,0,78,74,72,0,82,82,81,0,82,82,81,0,81,0,92,92,81,0,92,92,79,72,68,0,81,0,104,104,81,0,104,81,0,111,111,79,74,72,0,0,104,81,0,121,121,79,0,125,125,119,81,0,125,125,119,0,134,134,0,137,134,0,140,140,138,0,144,140,138,0,148,148,147,0,148,148,147,0,147,0,158,158,147,0,158,158,145,138,134,0,147,0,170,170,147,0,170,147,0,177,177,145,140,138,0,0,170,147,0,187,187,145,0,191,191,185,147,0,191,191,185,0,200,200,0,203,200,0,206,206,204,0,210,206,204,0,214,214,213,0,214,214,213,0,213,0,224,224,213,0,224,224,211,204,200,0,213,0,236,236,213,0,236,213,0,243,243,211,206,204,0,0,236,213,0,253,253,211,0,257,257,251,213,0,257,257,251,0,266,266,0,269,266,0,272,272,270,0,276,272,270,0,280,280,279,0,280,280,279,0,279,0,290,290,279,0,290,290,277,270,266,0,279,0,302,302,279,0,302,279,0,309,309,277,272,270,0,0,302,279,0,319,319,277,0,323,323,317,279,0,323,323,317,0,332,332,0,335,332,0,338,338,336,0,342,338,336,0,346,346,345,0,346,346,345,0,345,0,356,356,345,0,356,356,343,336,332,0,345,0,368,368,345,0,368,345,0,375,375,343,338,336,0,0,368,345,0,385,385,343,0,389,389,383,345,0,389,389,383
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
