#ifndef RANDLIB_ERRORS_H
#define RANDLIB_ERRORS_H

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

}

#endif // RANDLIB_ERRORS_H
