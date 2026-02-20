#pragma once
#include "SequenceNotation.h"

class PPSNotation : public SequenceNotation
{
public:
    // Check if successor ordinal (last term is 0)
    bool _isSuccessor() override;

    // Cut the last number
    void _cut() override;

    // Cut the last n numbers
    void _cut(INT n) override;

    // Expand non-successor n times
    void _expand(INT n) override;

    // Expand until length is at least x + n, then cut to x + n
    void _expandLen(INT n) override;

    // Same as _expandLen(0)
    void _reduce() override;

    // Dictionary order comparison
    int _compare(const SequenceNotation& other) const override;

    // Initialize to [0, 1, ..., n-1]
    void _init(INT n) override;

    // Find standard form and check non-maximum
    bool _checkStandardAndNonMaximum(std::vector<bool>& res) override;
    INT _expandUntilLarger(const SequenceNotation& target, bool selfcheck = true) override;
    void print(std::ostream& os) const override;
};
