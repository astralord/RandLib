#include "data.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace randlib_test
{

double parseNumber(const std::string &text)
{
    if (text == "inf")
        return std::numeric_limits<double>::infinity();
    if (text == "-inf")
        return -std::numeric_limits<double>::infinity();
    if (text == "nan" || text == "-nan")
        return std::numeric_limits<double>::quiet_NaN();
    char *end = nullptr;
    double value = std::strtod(text.c_str(), &end);
    if (end == text.c_str())
        throw std::invalid_argument("not a number: " + text);
    return value;
}

double Record::number() const
{
    return parseNumber(value);
}

std::string dataDirectory()
{
    return RANDLIB_TEST_DATA_DIR;
}

std::vector<std::string> splitCsvLine(const std::string &line)
{
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (quoted) {
            if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') {
                field += '"';
                ++i;
            } else if (c == '"') {
                quoted = false;
            } else {
                field += c;
            }
        } else if (c == '"') {
            quoted = true;
        } else if (c == ',') {
            fields.push_back(field);
            field.clear();
        } else if (c != '\r') {
            field += c;
        }
    }
    fields.push_back(field);
    return fields;
}

std::vector<Record> loadRecords(const std::string &subdirectory)
{
    namespace fs = std::filesystem;
    std::vector<fs::path> files;
    for (const auto &entry : fs::directory_iterator(fs::path(dataDirectory()) / subdirectory)) {
        if (entry.path().extension() == ".csv" && entry.path().filename() != "known_deviations.csv")
            files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());

    std::vector<Record> records;
    for (const auto &path : files) {
        std::ifstream in(path);
        std::string line;
        std::getline(in, line);
        auto header = splitCsvLine(line);
        auto column = [&header, &path](const std::string &name) {
            auto it = std::find(header.begin(), header.end(), name);
            if (it == header.end())
                throw std::runtime_error(path.string() + ": missing column " + name);
            return static_cast<size_t>(it - header.begin());
        };
        size_t caseCol = column("case"), functionCol = column("function"), argCol = column("arg"),
               valueCol = column("value");
        auto statusIt = std::find(header.begin(), header.end(), "status");
        while (std::getline(in, line)) {
            if (line.empty())
                continue;
            auto fields = splitCsvLine(line);
            Record r;
            r.caseId = fields.at(caseCol);
            r.function = fields.at(functionCol);
            r.arg = fields.at(argCol);
            r.value = fields.at(valueCol);
            r.status = (statusIt == header.end()) ? "ok" : fields.at(static_cast<size_t>(statusIt - header.begin()));
            records.push_back(std::move(r));
        }
    }
    return records;
}

} // namespace randlib_test
