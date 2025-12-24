#ifndef UTILS_HH
#define UTILS_HH

#include <iostream>
#include <cmath>
#include <functional>
#include <optional>
#include <sstream>
#include <string>

// in the .h file you declare the functions, classes, ..., while in the .cpp you define them

using type_fun = std::function<double(double const &)>;

struct Root_finder_params
{
    type_fun function_root;
    int max_iter = 1000;
    double rtol = 1e-6;
    double stol = 1e-6;
};

class Root_finder
{
public:
    Root_finder(type_fun, int, double, double);
    Root_finder(Root_finder_params);

    virtual ~Root_finder() = default; // Always make virtual the destructor when working with class inheritance

    virtual double solve() = 0;                  // Pure virtual class
    constexpr virtual const char *getName() = 0; // Make also this virtual to resolve at compile time the name
    virtual const std::string toString() = 0;
    std::vector<double> getXs();
    int getIter();
    double getRes();

    double current_root_estimate;
    Root_finder_params params;
    int iter = 0;
    double residual;
    double delta_x;
    std::vector<double> xs;
};

class Newton : public Root_finder
{
private:
    double initial_guess;
    type_fun der_fun;
    double step_size;
    bool is_Newton;

public:
    Newton(const double, type_fun, int, double, double, std::optional<type_fun> der_func_in = std::nullopt, double step_size = 1e-5);

    virtual ~Newton() override = default;
    double solve() override
    {
        // Implement here the method!
        current_root_estimate = initial_guess;
        do
        {
            delta_x = params.function_root(current_root_estimate) / der_fun(current_root_estimate);
            current_root_estimate -= delta_x;
            xs.push_back(current_root_estimate);
            residual = std::abs(params.function_root(current_root_estimate));
        } while (++iter < params.max_iter && std::abs(delta_x) > params.stol && residual > params.rtol);
        return (current_root_estimate);
    };
    constexpr virtual const char *getName() override { return (is_Newton ? "Newton" : "Secant"); };

    virtual const std::string toString() override
    {
        std::ostringstream oss;
        oss << getName() << " method\n"
            << " - approximate root : " << current_root_estimate << "\n"
            << " - # iterations     : " << iter << "\n"
            << " - residual         : " << residual << "\n";
        return oss.str();
    }
};

class Bisection : public Root_finder
{
private:
    double a, b;

public:
    Bisection(const double, const double, type_fun, int, double, double);
    virtual ~Bisection() override = default;
    double solve() override
    {
        // Check if root is bracketed
        if (params.function_root(a) * params.function_root(b) > 0)
        {
            std::cout << "Bisection method cannot find the zero of this function" << std::endl;
            exit(1);
        }

        do
        {
            current_root_estimate = .5 * (a + b);

            // Calculate function values
            double fa = params.function_root(a);
            double fc = params.function_root(current_root_estimate);

            // --- FIXED LOGIC ---
            // Check where the sign change occurs
            if (fa * fc < 0)
            {
                // Root is in the left half [a, c]
                b = current_root_estimate;
            }
            else
            {
                // Root is in the right half [c, b]
                a = current_root_estimate;
            }
            // -------------------

            delta_x = std::abs(b - a);
            xs.push_back(current_root_estimate);
            residual = std::abs(fc);

        } while (++iter < params.max_iter && std::abs(delta_x) > params.stol && residual > params.rtol);

        return (current_root_estimate);
    };
    constexpr virtual const char *getName() override { return "Bisection"; };
    virtual const std::string toString() override
    {
        std::ostringstream oss;
        oss << getName() << " method\n"
            << " - approximate root : " << current_root_estimate << "\n"
            << " - # iterations     : " << iter << "\n"
            << " - residual         : " << residual << "\n";
        return oss.str();
    }
};

#endif
