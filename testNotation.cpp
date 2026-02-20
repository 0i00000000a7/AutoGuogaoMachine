#include "PPSNotation.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

void print_seq(const std::vector<INT>& seq) {
    std::cout << "Current Sequence: (";
    for (size_t i = 0; i < seq.size(); ++i) {
        std::cout << seq[i] << (i == seq.size() - 1 ? "" : ", ");
    }
    std::cout << ")" << std::endl;
}

int main_testpps() {
    PPSNotation pps;
    pps.seq = {0, 1, 2, 3, 4, 5};

    std::cout << "PPS Notation Tester" << std::endl;
    std::cout << "Commands:" << std::endl;
    std::cout << "  c       : cut last term" << std::endl;
    std::cout << "  c<n>    : cut last n terms (e.g. c100)" << std::endl;
    std::cout << "  e<n>    : expand n times (e.g. e5)" << std::endl;
    std::cout << "  r       : reduce" << std::endl;
    std::cout << "  check   : run _checkStandardAndNonMaximum" << std::endl;
    std::cout << "  expand_until <target_seq>: run _expandUntilLarger (e.g. expand_until 0,1,2,3,4,5)" << std::endl;
    std::cout << "  read <seq>: read sequence (e.g. read 0,1,2,3,4,5)" << std::endl;
    std::cout << "  init <n>: run _init(n)" << std::endl;
    std::cout << "  reset   : reset to [0, 1, 2, 3, 4, 5]" << std::endl;
    std::cout << "  exit    : exit program" << std::endl;

    std::string line;
    while (true) {
        print_seq(pps.seq);
        std::cout << "> ";
        if (!std::getline(std::cin, line)) break;
        if (line.empty()) continue;

        if (line == "exit") break;

        if (line == "reset") {
            pps.seq = {0, 1, 2, 3, 4, 5};
            std::cout << "Sequence reset." << std::endl;
            continue;
        }

        if (line.substr(0, 4) == "read") {
            std::string nums_str = line.length() > 5 ? line.substr(5) : "";
            std::vector<INT> new_seq;
            std::stringstream ss(nums_str);
            std::string segment;
            bool error = false;
            while(std::getline(ss, segment, ',')) {
                try {
                    new_seq.push_back(std::stoll(segment));
                } catch (...) {
                    error = true;
                    break;
                }
            }
            if (error || new_seq.empty()) {
                std::cout << "Invalid sequence format. Usage: read 0,1,2,3" << std::endl;
            } else {
                pps.seq = new_seq;
                std::cout << "Sequence updated." << std::endl;
            }
            continue;
        }

        if (line.substr(0, 12) == "expand_until") {
            std::string nums_str = line.substr(12);
            std::vector<INT> target_seq;
            std::stringstream ss(nums_str);
            std::string segment;
            bool error = false;
            while(std::getline(ss, segment, ',')) {
                try {
                    target_seq.push_back(std::stoll(segment));
                } catch (...) {
                    error = true;
                    break;
                }
            }
            if (error || target_seq.empty()) {
                std::cout << "Invalid target sequence. Usage: expand_until 0,1,2,3" << std::endl;
                continue;
            }

            PPSNotation target;
            target.seq = target_seq;
            try {
                bool res = pps._expandUntilLarger(target);
                std::cout << "Result: " << (res ? "True" : "False") << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
            continue;
        }

        if (line.substr(0, 4) == "init") {
             try {
                INT n = std::stoll(line.substr(5));
                pps._init(n);
                std::cout << "Initialized to fundamental sequence of length " << n << std::endl;
            } catch (...) {
                std::cout << "Invalid init command. Usage: init <n>" << std::endl;
            }
            continue;
        }

        if (line == "check") {
            std::vector<bool> res(pps.seq.size(), false);
            bool result = pps._checkStandardAndNonMaximum(res);
            std::cout << "Check result: " << (result ? "True" : "False") << std::endl;
            if (result) {
                std::cout << "Res: (";
                for (size_t i = 0; i < res.size(); ++i) {
                    std::cout << (res[i] ? "T" : "F") << (i == res.size() - 1 ? "" : ", ");
                }
                std::cout << ")" << std::endl;
            }
            continue;
        }

        if (line == "r") {
            if (pps._isSuccessor()) {
                std::cout << "Successor detected. Cutting last term." << std::endl;
            }
            pps._reduce();
            continue;
        }

        if (line[0] == 'c' && line != "check") {
            if (line.length() == 1) {
                pps._cut();
            } else {
                try {
                    INT n = std::stoll(line.substr(1));
                    pps._cut(n);
                } catch (...) {
                    std::cout << "Invalid cut command." << std::endl;
                }
            }
            continue;
        }

        if (line[0] == 'e') {
            try {
                INT n = std::stoll(line.substr(1));
                if (pps._isSuccessor()) {
                    std::cout << "Successor detected. Cutting last term." << std::endl;
                }
                pps._expand(n);
            } catch (...) {
                std::cout << "Invalid expand command. Usage: e<number>" << std::endl;
            }
            continue;
        }

        std::cout << "Unknown command." << std::endl;
    }

    return 0;
}
