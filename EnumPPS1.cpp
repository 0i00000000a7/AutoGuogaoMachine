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
    0,1,0,2,2,0,5,2,0,8,8,6,0,12,8,6,0,16,16,15,0,16,16,15,0,15,0,26,26,15,0,26,26,13,6,2,0,15,0,38,38,15,0,38,15,0,45,45,13,8,6,0,51,43,15,0,0,45,45,13,0,60,60,0,63,60,0,66,66,64,0,70,66,64,0,74,74,73,0,74,74,73,0,73,0,84,84,73,0,84,84,71,64,60,0,73,0,96,96,73,0,96,73,0,103,103,71,66,64,0,109,101,73,0,0,103,103,71,0,118,118,0,121,118,0,124,124,122,0,128,124,122,0,132,132,131,0,132,132,131,0,131,0,142,142,131,0,142,142,129,122,118,0,131,0,154,154,131,0,154,131,0,161,161,129,124,122,0,167,159,131,0,0,161,161,129,0,176,176,0,179,176,0,182,182,180,0,186,182,180,0,190,190,189,0,190,190,189,0,189,0,200,200,189,0,200,200,187,180,176,0,189,0,212,212,189,0,212,189,0,219,219,187,182,180,0,225,217,189,0,0,219,219,187,0,234,234,0,237,234,0,240,240,238,0,244,240,238,0,248,248,247,0,248,248,247,0,247,0,258,258,247,0,258,258,245,238,234,0,247,0,270,270,247,0,270,247,0,277,277,245,240,238,0,283,275,247,0,0,277,277,245,0,292,292,0,295,292,0,298,298,296,0,302,298,296,0,306,306,305,0,306,306,305,0,305,0,316,316,305,0,316,316,303,296,292,0,305,0,328,328,305,0,328,305,0,335,335,303,298,296,0,341,333,305,0,0,335,335,303
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
