#pragma once

#include <string>
#include <vector>

namespace str_utils {

std::vector<std::string> split(std::string s, const std::string& delimiter);

std::vector<std::vector<std::string>> parse_matrix(std::string_view str, char delim, char br_op, char br_cl);

}  // namespace str_utils