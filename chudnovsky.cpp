/* Comput π using the Chudnovsky Algorithm:
 *
 * 1/π = 12 Σ ((-1)^k * (6k)! * (545140134k + 13591409)) / ((3k)!(k!)^3 * (640320)^(3k+3/2))
 *
 * Where k = 0 => ∞
 *
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <iostream>
#include <execution>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <chrono>
#include <gmp.h>

void calculatePI(int, mpf_t&);

using namespace std::chrono;

class ExecutionTimer
{
public:
    // Use the best steady clock available
    using Clock = std::conditional_t<high_resolution_clock::is_steady,
        high_resolution_clock,
        steady_clock>;

    ~ExecutionTimer()
    {
        std::string units = " microSeconds";
        // Determine whether to print uSecs or mSecs or Secs
        double count = duration_cast<microseconds>(Clock::now() - mStart).count();
        if (count > 1000)
        {
            // Convert to milliseconds
            units = " milliSeconds";
            count /= 1000.0f;
            if (count > 1000)
            {
                // Convert to seconds
                units = " Seconds";
                count /= 1000.0f;
            }
        }

        std::cout
            << "Elapsed: " << count << units.data() << std::endl;
    }
private:
    Clock::time_point mStart = Clock::now();
};

using boost::multiprecision::cpp_int;

static constexpr int log2n = 4;

static cpp_int factorial(const cpp_int& num) 
{
    cpp_int fact = 1;
    for(cpp_int i = 1; i <= num; ++i)
        fact *= i;
    return fact;
}

static cpp_int numerator(const cpp_int& k) 
{
	const auto six_k_fact = factorial(6*k);
    
    return (k & 1 ? -1 : 1) * six_k_fact * (545140134*k + 13591409);
}

static cpp_int denominator_a(const cpp_int& k) 
{
	const auto factorial_k = factorial(k);
    return factorial(3*k) * factorial_k * factorial_k * factorial_k;
}

static cpp_int pow_3k(int64_t k) 
{
    cpp_int base = 640320;
    cpp_int ret = 1;
    int64_t exponent = 3 * k;

    for (int64_t i = 0; i < exponent; ++i)
        ret *= base;

    return ret;
}

static cpp_int newpow_3k(int64_t k)
{
    cpp_int base = 640320;
    cpp_int ret = 1;
    int64_t exponent = 3 * k;

    if (exponent < 0)
    {
        base = 1 / base;
        exponent = -exponent;
    }

    while (exponent > 0)
    {
        if (exponent % 2 == 1)
        {
            // If exponent is odd
            ret *= base;
        }
        base *= base; // Square the base
        exponent /= 2; // Divide exponent by 2
    }
    return ret;

}

int calcPrecision(int num_terms) 
{
    constexpr int PLACES_PER_TERM = 20;

    return num_terms * PLACES_PER_TERM;
}

// Function to calculate multiple of sigma series with required precision
boost::multiprecision::mpfr_float calcConstant(int precision)
{
    using boost::multiprecision::mpfr_float;
    mpfr_float::default_precision(precision * log2n);
    const mpfr_float numerator = 1;
    const mpfr_float denominator_a = 426880;
    const mpfr_float denominator_b_squared = 10005;
    mpfr_float result = numerator / (denominator_a * sqrt(denominator_b_squared));
    return result;
}


int main(const int argc, char* argv[]) 
{
    if(argc != 2) 
    {
        std::cerr << "Usage: " << argv[0] << " number_of_terms\n";
        return 1;
    }
    using boost::multiprecision::mpfr_float;

    try {
	    const int64_t num_terms = std::stoi(argv[1]);
        
        if (num_terms <= 0) 
        {
            throw std::invalid_argument("Number of terms must be greater than zero.");
        }


        mpfr_float pi;
        const int precision = calcPrecision(num_terms);

        mpfr_float::default_precision(precision * log2n);

        {
	        mpfr_float constant = calcConstant(precision);
	        mpfr_float pi_inverse = 0;
	        ExecutionTimer timer;
#if 1
            unsigned num_threads = std::thread::hardware_concurrency();

            std::vector<mpfr_float> partial_sums(num_threads);

            std::for_each(std::execution::par_unseq, partial_sums.begin(), partial_sums.end(), [&](mpfr_float& local_sum)
                {
	                const auto thread_id = &local_sum - partial_sums.data();
                    for (int64_t k = thread_id; k < num_terms; k += num_threads)
                    {
	                    const mpfr_float temp = mpfr_float(numerator(k)) / mpfr_float(denominator_a(k) * pow_3k(k));
                        local_sum += temp;
                    }
                });

            // Combine the partial sums
            for (const auto& sum : partial_sums)
            {
                pi_inverse += sum;
            }

#else
            unsigned num_threads = std::thread::hardware_concurrency();
            std::vector<std::thread> threads(num_threads);
			std::mutex mtx;

            for (unsigned i = 0; i < num_threads; ++i)
            {
                threads[i] = std::thread([&, i]
                    {
                        for (int64_t k = i; k < num_terms; k += num_threads)
                        {
                            const mpfr_float temp = mpfr_float(numerator(k)) / mpfr_float(denominator_a(k) * pow_3k(k));
                            {
                                std::lock_guard lock(mtx);
                                pi_inverse += temp;
                            }
                        }
                    });
            }
            for (auto& th : threads)
                th.join();


#endif
            pi = mpfr_float(1) / (pi_inverse * constant);

        }

        std::cout << std::setprecision(precision) << pi << "\n";

    } 
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}




void calculatePI(int pi_len, mpf_t &_pi)
{

    /*
       mpf_t is the type defined for GMP integers.
       It is a pointer to the internals of the GMP integer data structure
     */


    mpf_t _prev; // previous eq. ans.
    mpf_t temp; // Temp var.
    mpf_t temp1; // Temp var.
    mpf_t one; // ONE
    mpf_t two; //TWO

    mpf_set_default_prec(pi_len);

    /* 1. Initialize the number */
    mpf_init(_pi);
    mpf_init(_prev);
    mpf_init(temp);
    mpf_init(temp1);
    mpf_init(one);
    mpf_init(two);

    mpf_set_d(_pi, 0);
    mpf_set_d(_prev, 0);
    mpf_set_d(temp, 0);
    mpf_set_d(temp1, 0);
    mpf_set_d(one, 1);
    mpf_set_d(two, 2);

    for (int i = pi_len; i > 0; i--) {
        mpf_set_d(temp, i);

        mpf_set(temp1, two);

        mpf_mul_ui(temp1, temp1, i);
        mpf_add_ui(temp1, temp1, 1);

        mpf_div(temp, temp, temp1);

        mpf_mul(_prev, _prev, temp);
        mpf_add(_pi, _prev, one);


        mpf_set(_prev, _pi);
    }

    mpf_mul(_pi, _pi, two);
    gmp_printf("PI = %.*Ff", pi_len / 14 + 2, _pi);

    printf("\n");
    /* 6. Clean up the mpf_t handles or else we will leak memory */
    mpf_clear(_pi);
    mpf_clear(_prev);
    mpf_clear(temp);
    mpf_clear(temp1);
    mpf_clear(one);
    mpf_clear(two);

}