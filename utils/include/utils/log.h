#pragma once

#include <fstream>
#include <iostream>

enum class LogMode { NONE, STDOUT, FILE };

class Logger {
public:
    Logger(LogMode mode = LogMode::NONE, std::string filename = "intmat.log") : mode_(mode), filename_(filename) {
        if (mode == LogMode::FILE) {
            fs_ = std::fstream{filename, fs_.out};
        }
    };
    Logger(Logger&& log) : mode_(log.mode_), fs_(std::move(log.fs_)), filename_(log.filename_) {};
    Logger(Logger& log) : mode_(log.mode_), fs_(log.filename_, log.fs_.out), filename_(log.filename_) {};

    Logger& operator=(const Logger& log) {
        mode_ = log.mode_;
        fs_ = std::fstream{log.filename_, fs_.out};
        filename_ = log.filename_;

        return *this;
    }

    template <typename T>
    Logger& operator<<(const T& t) {
        switch (mode_) {
        case LogMode::NONE:
        break;
        case LogMode::STDOUT: {
            std::cout << t <<'\n';
            break;
        }
        case LogMode::FILE: {
            fs_ << t << '\n';
            fs_.flush();
            break;
        }
        }

        return *this;
    }

private:
    LogMode mode_;
    std::fstream fs_;
    std::string filename_;
};
