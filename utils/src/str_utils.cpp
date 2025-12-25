#include "utils/str_utils.h"

#include <algorithm>

#include <iostream>
namespace str_utils {
std::vector<std::string> split(std::string s, const std::string& delimiter) {
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    tokens.push_back(s);

    return tokens;
}

std::vector<std::vector<std::string>> parse_matrix(std::string_view str, char delim, char br_op, char br_cl) {
    std::vector<std::vector<std::string>> res;

    auto start = str.find(br_op);
    if (start == std::string::npos) {
        return {};
    }
    ++start;
    int64_t opened = 1;
    size_t pos = start;
    while (pos < str.size()) {
        if (str[pos] == br_op) {
            opened++;
            auto row_begin = pos + 1;
            if (row_begin >= str.size()) {
                return {};
            }

            while (str[pos] != br_cl && pos < str.size()) {
                pos++;
            }

            if (pos >= str.size()) {
                return {};
            }

            std::string row_str(pos - row_begin, '\0');
            str.copy(row_str.data(), pos - row_begin, row_begin);
            row_str.erase(std::remove_if(row_str.begin(),
                                         row_str.end(),
                                         [](auto&& e) {
                                             return std::isspace(e);
                                         }),
                          row_str.end());
            res.push_back(split(row_str, std::string{ delim }));
            opened--;
        } else if (str[pos] == br_cl) {
            opened--;
        }

        pos++;
    }

    if (opened != 0) {
        return {};
    }

    return res;
}
}  // namespace str_utils
