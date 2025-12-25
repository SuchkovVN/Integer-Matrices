#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>

typedef std::vector<std::vector<std::string>> StrMatrix;

struct ParamBase {
    std::string name;
    std::string description;

    ParamBase(std::string_view name, std::string_view description) : name(name), description(description) {}

    virtual ~ParamBase() = default;
    virtual void parse(std::string_view) = 0;
};

template <typename Ty>
struct Param : public ParamBase {
    Ty* value;

    Param(std::string_view name, std::string_view description, Ty* value, Ty defaultValue) : ParamBase(name, description), value(value) {
        if (value == nullptr) {
            throw std::runtime_error("Invalid param value ptr");
        }
        *value = defaultValue;
    }
    void parse(std::string_view str) override;
};

class Params {
public:
    Params() = default;
    ~Params() {
        for (auto& p : _params) {
            delete p.second;
        }
        _params.clear();
    }

    template <typename Ty>
    void createParam(std::string_view name, Ty defaultValue, Ty* value, std::string_view description = "") {
        _params.try_emplace(std::string(name), new Param<Ty>(std::string(name), std::string(description), value, defaultValue));
    }

    ParamBase* get(std::string_view name);

    void printParams();

private:
    std::unordered_map<std::string, ParamBase*> _params;
};

bool parseParams(char** argv, int argc, Params& params);
