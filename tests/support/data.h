#ifndef RANDLIB_TESTS_DATA_H
#define RANDLIB_TESTS_DATA_H

#include <string>
#include <vector>

namespace randlib_test
{

/// One row of tests/data/golden_v1/*.csv or tests/data/reference/*.csv.
struct Record
{
    std::string caseId;
    std::string function;
    std::string arg;
    std::string value; ///< textual value as stored
    std::string status; ///< "ok", "throw:<type>" or "timeout" (golden only; "ok" for reference)

    double number() const;
};

/// Parses a number written by the data generators, including "inf", "nan" and
/// subnormal values (which std::stod rejects as out of range).
double parseNumber(const std::string &text);

/// Directory with the committed test data (tests/data).
std::string dataDirectory();

/// All rows of every CSV file in tests/data/<subdirectory>.
std::vector<Record> loadRecords(const std::string &subdirectory);

/// Splits one CSV line, honouring double-quoted fields.
std::vector<std::string> splitCsvLine(const std::string &line);

} // namespace randlib_test

#endif // RANDLIB_TESTS_DATA_H
