#pragma once

#include <stdexcept>
#include <string>

namespace randlib
{

/**
 * @brief The not_implemented class
 * Thrown when a quantity is not yet available for the given parameters, although the
 * distribution itself is valid (for example, the stable density for exponents close to 1).
 */
class not_implemented : public std::runtime_error
{
public:
    explicit not_implemented(const std::string &what) : std::runtime_error(what) {}
};

/// Reason why a distribution could not be fitted to a sample
enum class FitErrorCode {
    WrongSample, ///< the sample could not have been drawn from this distribution
    NotApplicable, ///< the method cannot be applied to this distribution or sample
    WrongReturn, ///< the method returned invalid parameters
    TooFewElements, ///< the sample is too small
    WrongLevel, ///< the significance level is not in (0, 1)
    Undefined ///< any other failure
};

/**
 * @brief The fit_error class
 * Thrown when a distribution cannot be fitted to a sample; code() tells why
 */
class fit_error : public std::runtime_error
{
public:
    fit_error(FitErrorCode code, const std::string &what) : std::runtime_error(what), errorCode(code) {}

    [[nodiscard]] FitErrorCode code() const noexcept { return errorCode; }

private:
    FitErrorCode errorCode;
};

} // namespace randlib
