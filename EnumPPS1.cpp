#include <fstream>
#include "PPSNotation.h"
#include <omp.h>
using namespace std;
const int ENUM_LAYERS = 95;
const int MAX_NONLOOP_LEN = ENUM_LAYERS / 4;
const int CHECK_COPY_NUM=5;
const int PRINT_THRESHOLD = ENUM_LAYERS - 11;

PPSNotation smallest_nonterminate_seq;
const vector<INT> known_smallest_seq =
{
0,1,0,2,2,0,5,2,0,8,8,6,0,12,8,6,0,16,16,15,0,16,16,15,0,15,0,26,26,15,0,26,26,13,6,2,0,15,0,38,38,15,0,38,15,0,45,45,13,8,6,0,0,38,15,0,55,55,13,0,59,57,13,0,0,59,53,15,0,0,13,0,71,55,0,74,74,72,55,0,74,74,72,0,83,83,0,86,83,0,89,89,87,0,93,89,87,0,97,97,96,0,97,97,96,0,96,0,107,107,96,0,107,107,94,87,83,0,96,0,119,119,96,0,119,96,0,126,126,94,89,87,0,0,119,96,0,136,136,94,0,140,138,94,0,0,140,134,96,0,0,94,0,152,136,0,155,155,153,136,0,155,155,153,0,164,164,0,167,164,0,170,170,168,0,174,170,168,0,178,178,177,0,178,178,177,0,177,0,188,188,177,0,188,188,175,168,164,0,177,0,200,200,177,0,200,177,0,207,207,175,170,168,0,0,200,177,0,217,217,175,0,221,219,175,0,0,221,215,177,0,0,175,0,233,217,0,236,236,234,217,0,236,236,234,0,245,245,0,248,245,0,251,251,249,0,255,251,249,0,259,259,258,0,259,259,258,0,258,0,269,269,258,0,269,269,256,249,245,0,258,0,281,281,258,0,281,258,0,288,288,256,251,249,0,0,281,258,0,298,298,256,0,302,300,256,0,0,302,296,258,0,0,256,0,314,298,0,317,317,315,298,0,317,317,315,0,326,326,0,329,326,0,332,332,330,0,336,332,330,0,340,340,339,0,340,340,339,0,339,0,350,350,339,0,350,350,337,330,326,0,339,0,362,362,339,0,362,339,0,369,369,337,332,330,0,0,362,339,0,379,379,337,0,383,381,337,0,0,383,377,339,0,0,337,0,395,379,0,398,398,396,379,0,398,398,396,0,407,407,0,410,407,0,413,413,411,0,417,413,411,0,421,421,420,0,421,421,420,0,420,0,431,431,420,0,431,431,418,411,407,0,420,0,443,443,420,0,443,420,0,450,450,418,413,411,0,0,443,420,0,460,460,418,0,464,462,418,0,0,464,458,420,0,0,418,0,476,460,0,479,479,477,460,0,479,479,477
};
const int PARALLEL_DEPTH = 85;
const int BASE_PREFIX_LEN = 75;

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
    int candidates[] = { (int)seq.size() - 3, (int)seq.size() - 6 };
    for (int loop_len : candidates) {
        if (loop_len < 1 || loop_len > (int)seq.size()) continue; // 跳过不合法的
        vector<INT> loop_dif(loop_len);
        for (int j = 0; j < loop_len; j++) {
            loop_dif[j] = seq[seq.size() - loop_len + j] > 0 ? loop_len : 0;
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
