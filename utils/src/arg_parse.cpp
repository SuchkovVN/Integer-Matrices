#include "utils/arg_parse.h"

#include "utils/str_utils.h"

#include <string>
#include <vector>
#include <iostream>

void Params::printParams() {
    for (auto p : _params) {
        auto* param = p.second;
        std::cout << "\"" << param->name << "\" - \"" << param->description << "\"\n";
    }
}

ParamBase* Params::get(std::string_view name) {
    auto it = _params.find(std::string(name));
    if (it == _params.end()) {
        return nullptr;
    }

    return it->second;
}

template <>
void Param<std::string>::parse(std::string_view str) {
    *value = std::string(str);
}

template <>
void Param<int64_t>::parse(std::string_view str) {
    *value = std::stoll(std::string(str));
}

template <>
void Param<bool>::parse(std::string_view str) {
    if (str == "0" || str == "false") {
        *value = false;
    } else if (str == "1" || str == "true") {
        *value = true;
    } else {
       throw std::runtime_error("Invalid bool param value");
    }
}

template <>
void Param<int>::parse(std::string_view str) {
    *value = std::stoi(std::string(str));
}

template <>
void Param<StrMatrix>::parse(std::string_view str) {
    *value = str_utils::parse_matrix(str, ',', '[', ']');
}

bool parseParams(char** argv, int argc, Params& params) {
    for (int i = 1; i < argc; ++i) {
        std::string_view s(argv[i]);

        if (s == "--help") {
            params.printParams();
            return false;
        }

        auto pos = s.find('=');
        if (pos == std::string_view::npos) {
            std::cout << "Invalid option: " << s << '\n';
            return false;
        }

        auto name = s.substr(0, pos);
        auto value = s.substr(pos + 1);

        auto p = params.get(name);
        if (p == nullptr) {
            std::cout << "Unrecognized option ignored: " << name << '\n';
            continue;
        }

        try {
            p->parse(value);
        } catch (const std::exception& ex) {
            std::cout << "Failed to parse param: " << ex.what();
            return false;
        } catch (...) {
            std::cout << "Failed to parse param due to unknown error";
            return false;
        }
    }

    return true;
}
