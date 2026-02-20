#pragma once
#include <vector>
#include <cstdint>
#include <cassert>
#include <iostream>

typedef int16_t INT;

class SequenceNotation
{
public:
    std::vector<INT> seq;

    virtual ~SequenceNotation() = default;

    // Interfaces
    virtual bool _isSuccessor() = 0;
    virtual void _cut() = 0;
    virtual void _cut(INT n) = 0;
    virtual void _expand(INT n) = 0;
    virtual void _expandLen(INT n) = 0;
    virtual void _reduce() = 0;
    
    // New Interfaces
    virtual int _compare(const SequenceNotation& other) const = 0;
    virtual void _init(INT n) = 0;
    virtual bool _checkStandardAndNonMaximum(std::vector<bool>& res) = 0;
    virtual INT _expandUntilLarger(const SequenceNotation& target, bool selfcheck = true) = 0;

    virtual void print(std::ostream& os) const = 0;

    // Operators
    bool operator>(const SequenceNotation& other) const {
        return _compare(other) > 0;
    }
    bool operator<(const SequenceNotation& other) const {
        return _compare(other) < 0;
    }
    bool operator>=(const SequenceNotation& other) const {
        return _compare(other) >= 0;
    }
    bool operator<=(const SequenceNotation& other) const {
        return _compare(other) <= 0;
    }
    bool operator==(const SequenceNotation& other) const {
        return _compare(other) == 0;
    }
    bool operator!=(const SequenceNotation& other) const {
        return _compare(other) != 0;
    }
};
