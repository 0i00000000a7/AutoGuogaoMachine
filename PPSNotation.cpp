#include "PPSNotation.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>

// Check if successor ordinal (last term is 0)
bool PPSNotation::_isSuccessor() {
    if (seq.empty()) return false;
    return seq.back() == 0;
}

// Cut the last number
void PPSNotation::_cut() {
    if (!seq.empty()) {
        seq.pop_back();
    }
}

// Cut the last n numbers
void PPSNotation::_cut(INT n) {
    if (n < 0) return;
    if ((size_t)n >= seq.size()) {
        seq.clear();
    } else {
        seq.resize(seq.size() - n);
    }
}

// Expand non-successor n times
// Fundamental sequence [n] expands to the y + nL - 1 term
void PPSNotation::_expand(INT n) {
    if (seq.empty()) return;
    
    // If successor, undefined by user for _expand(n), but n=0 is _cut.
    // Assuming standard behavior for successor is just to cut (predecessor).
    if (_isSuccessor()) {
        _cut();
        return;
    }

    if (n == 0) {
        _cut();
        return;
    }

    // 1. Identify parameters
    // User definitions:
    // Last term value: x
    // Last term index (1-based): y
    // Bad root: The x-th term
    // L = y - x
    
    INT x = seq.back();
    INT y = (INT)seq.size();
    
    // 0-based index for bad root is x - 1
    INT bad_idx = x - 1;

    // Validation
    if (bad_idx < 0 || bad_idx >= y - 1) {
        std::cout << "Invalid PPS structure or logic error: bad_idx out of range." << std::endl;
        throw std::invalid_argument("bad_idx out of range");
        return;
    }

    INT b = seq[bad_idx]; // Bad root value
    INT L = y - x;        // Period length

    // 2. Last term expansion
    // "If last term and bad root between (exclusive) exists a term == b, set last = b"
    // Range (0-based): (bad_idx + 1) to (y - 2)
    
    bool found_b = false;
    for (INT i = bad_idx + 1; i <= y - 2; ++i) {
        if (seq[i] == b) {
            found_b = true;
            break;
        }
    }

    if (found_b) {
        seq.back() = b;
    } else {
        seq.back() = x - 1;
    }

    // 3. Other terms expansion
    // Target length: y + nL - 1
    // We need to generate terms starting from index y up to target_len - 1
    
    INT target_len = y + n * L - 1;
    
    if (target_len > y) {
        seq.reserve(target_len);
        for (INT k = y; k < target_len; ++k) {
            // Recurrence: "For any i > y - L + 1... i+L term..."
            // Maps to: seq[k] is derived from seq[k-L]
            // If seq[k-L] >= x (original x), then + L
            
            INT source_idx = k - L;
            // source_idx should be >= bad_idx + 1 usually
            
            INT val = seq[source_idx];
            if (val >= x) {
                seq.push_back(val + L);
            } else {
                seq.push_back(val);
            }
        }
    } else if (target_len < y) {
        // Possible if n=1 and L=1 -> y + 1 - 1 = y (no append)
        // If target_len < y (unlikely for n>=1, L>=1), resize down
        seq.resize(target_len);
    }
}

// Expand until length is at least x + n, then cut to x + n
// Here 'x' refers to the original length (as per user context: "If original length is x...")
void PPSNotation::_expandLen(INT n) {
    if (seq.empty()) return;

    if (_isSuccessor()) {
        _cut();
        return;
    }

    INT orig_len = (INT)seq.size();
    INT target_len = orig_len + n;

    // Calculate parameters
    INT x = seq.back();
    INT y = orig_len;
    INT bad_idx = x - 1;

    if (bad_idx < 0 || bad_idx >= y - 1) return;

    INT b = seq[bad_idx];
    INT L = y - x;

    // Modification Step (Always happens for non-successor expansion)
    bool found_b = false;
    for (INT i = bad_idx + 1; i <= y - 2; ++i) {
        if (seq[i] == b) {
            found_b = true;
            break;
        }
    }

    if (found_b) {
        seq.back() = b;
    } else {
        seq.back() = x - 1;
    }

    // Append until length >= target_len
    // If n=0, target_len = orig_len. We are done.
    
    if (target_len > orig_len) {
        seq.reserve(target_len);
        for (INT k = orig_len; k < target_len; ++k) {
            INT source_idx = k - L;
            INT val = seq[source_idx];
            if (val >= x) {
                seq.push_back(val + L);
            } else {
                seq.push_back(val);
            }
        }
    } else {
        // Just in case modification resulted in correct length (always true for n=0)
        // But if target_len < orig_len (impossible since n>=0), we would cut.
    }
}

// Same as _expandLen(0)
void PPSNotation::_reduce() {
    _expandLen(0);
}

// Dictionary order comparison
int PPSNotation::_compare(const SequenceNotation& other) const {
    size_t len1 = seq.size();
    size_t len2 = other.seq.size();
    size_t min_len = std::min(len1, len2);

    for (size_t i = 0; i < min_len; ++i) {
        if (seq[i] > other.seq[i]) return 1;
        if (seq[i] < other.seq[i]) return -1;
    }

    if (len1 > len2) return 1;
    if (len1 < len2) return -1;
    return 0;
}

// Initialize to [0, 1, ..., n-1]
void PPSNotation::_init(INT n) {
    seq.clear();
    seq.reserve(n);
    for (INT i = 0; i < n; ++i) {
        seq.push_back(i);
    }
}

// Find standard form and check non-maximum
bool PPSNotation::_checkStandardAndNonMaximum(std::vector<bool>& res) {
    // res must be same size as current seq and initially false
    res.resize(seq.size());
    std::fill(res.begin(), res.end(), false);

    // 1. Find the first fundamental sequence > current sequence
    // The user description: "Start from the initial fundamental sequence, find the first one > current sequence"
    // However, the example "[0,1,2,2,1] find [0,1,2,3] > [0,1,2,2,1]" implies finding a common prefix + next element > current.
    
    // Actually, the user says "From initial fundamental sequence...".
    // Initial fundamental sequence for length N is usually [0, 1, ..., N-1]? Or maybe the one constructed by _init(N)?
    // The user example: For [0,1,2,2,1], find [0,1,2,3].
    // [0,1,2,3] has length 4. [0,1,2,2,1] has length 5.
    // [0,1,2,3] is > [0,1,2,2,1] because at index 3: 3 > 2.
    
    // It seems we should iterate through indices.
    // The user logic implies we want to find a "parent" or "upper bound" sequence that can expand to the current one.
    
    // "Find the first difference, cut to this position and start _expand, mark this position true."
    
    // Let's interpret "Start from initial fundamental sequence" as constructing a standard sequence that is "larger" or "simpler" 
    // and then refining it down.
    
    // Actually, looking at the example: [0,1,2,2,1]. First diff is at index 3 (value 2). 
    // The standard sequence [0,1,2,3] differs at index 3 (value 3).
    // So we found an index `i` such that `standard[i] > current[i]`.
    
    // Strategy:
    // We maintain a "working" sequence `work`. Initially it should be a "standard" sequence.
    // What is the "standard" sequence? Maybe [0, 1, 2, 3, ...]
    
    PPSNotation work;
    // We don't know the exact length to start with.
    // But based on "find the first difference", we can compare [0, 1, 2, ...] with `seq`.
    
    work._init(seq.size());
    
    INT idx=0;

    while(idx < seq.size()) {
        if(work.seq[idx] < seq[idx]) //non-standard
            return false;
        if(work.seq[idx] == seq[idx]) 
        {
            idx++;
            continue;
        }
        assert(work.seq[idx] > seq[idx]);
        res[idx] = true;

        //cut
        work.seq.resize(idx + 1);
        
        INT expanded_terms=work._expandUntilLarger(*this, false);
        if(expanded_terms == -1)
        {
            return false;
        }
        idx += expanded_terms;
    }
    return true;
}

INT PPSNotation::_expandUntilLarger(const SequenceNotation& target, bool selfcheck) {
    if (selfcheck) {
        if (seq.size() > target.seq.size()) {
            throw std::runtime_error("Sequence length larger than target length");
        }
        for (size_t i = 0; i < seq.size() - 1; ++i) {
            if (seq[i] != target.seq[i]) {
                throw std::runtime_error("Sequence prefix does not match target");
            }
        }
        if(seq.empty())
        {
            throw std::runtime_error("Sequence is empty");
        }
        if (_isSuccessor()) {
            throw std::runtime_error("Sequence is a successor");
        }

        
    }


    INT orig_len = (INT)seq.size();
    
    // Capture parameters from original sequence
    INT x = seq.back();
    INT y = orig_len;
    INT bad_idx = x - 1;

    if (bad_idx < 0 || bad_idx >= y - 1) return false;

    INT b = seq[bad_idx];
    INT L = y - x;

    // Modification Step
    bool found_b = false;
    for (INT i = bad_idx + 1; i <= y - 2; ++i) {
        if (seq[i] == b) {
            found_b = true;
            break;
        }
    }

    if (found_b) {
        seq.back() = b;
    } else {
        seq.back() = x - 1;
    }


    INT val = seq.back();
    INT target_val = target.seq[orig_len - 1];

    if (val > target_val) return 0;
    if (val < target_val) return -1;

    // If equal, we must append terms until we differ or reach target length
    assert(val == target_val);
    // If lengths are equal and we match exactly
    if (orig_len == target.seq.size()) return 0;

    // Append loop
    for (INT k = orig_len; k < (INT)target.seq.size(); ++k) {
        INT source_idx = k - L;
        INT next_val = seq[source_idx];
        if (next_val >= x) {
            next_val += L;
        }
        seq.push_back(next_val);
        
        if (next_val > target.seq[k]) return seq.size()-orig_len; //larger
        if (next_val < target.seq[k]) return -1; //smaller, non-standard
        //same, continue
    }
    
    // If we finished loop, we are exactly equal to target
    return seq.size()-orig_len;


}
void PPSNotation::print(std::ostream& os) const {
    for (INT i = 0; i < (INT)seq.size(); ++i) {
        os << seq[i] << ",";
    }
}
