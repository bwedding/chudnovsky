/* Comput π using the Chudnovsky Algorithm:
 *
 * 1/π = 12 Σ ((-1)^k * (6k)! * (545140134k + 13591409)) / ((3k)!(k!)^3 * (640320)^(3k+3/2))
 *
 * Where k = 0 => ∞
 *
 */

#include <algorithm>
#include <chrono>
#include <execution>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <thread>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <gmp.h>

using namespace std::chrono;
using boost::multiprecision::cpp_int;
using boost::multiprecision::mpfr_float;

const std::string filename = "pi-cache.bin";
const std::string keyTrackingFile = "last_key.txt";

static constexpr int log2n = 4;
std::mutex cache_mutex;

static cpp_int numerator(const cpp_int& k);
static cpp_int denominator_a(const cpp_int& k);

int64_t initializeLastWrittenKey(const std::string& keyTrackingFile)
{
    std::ifstream file(keyTrackingFile);
    int64_t lastKey = 0;
    if (file.is_open()) 
    {
        file >> lastKey;
    }
    file.close();
    return lastKey;
}

// Optional: Function to save the last written key for future use
void saveLastWrittenKey(const std::string& keyTrackingFile, int64_t lastWrittenKey)
{
    std::ofstream file(keyTrackingFile);
    if (file.is_open()) 
    {
        file << lastWrittenKey;
    }
    file.close();
}

void appendNewEntriesToFile(const std::vector<mpfr_float>& cache, const std::string& filename, int64_t& lastWrittenKey)
{
    std::ofstream file(filename, std::ios::app); // Open file in append mode
    if (!file.is_open()) 
    {
        std::cerr << "Failed to open file for appending.\n";
        return;
    }

    // Iterate over new elements based on lastWrittenKey
    for (size_t i = lastWrittenKey; i < cache.size(); ++i) 
    {
        file << (i + 1) << "," << cache[i].str() << std::endl; // i + 1 because keys start from 1
        lastWrittenKey = i + 1; // Update lastWrittenKey to match the key just written
    }

    file.close();
}

void writeObjectsToFile(const std::vector<mpfr_float>& cache, const std::string& filename, int64_t& lastWrittenKey)
{
    std::ofstream file(filename, std::ios::binary | std::ios::app);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file for appending.\n";
        return;
    }

    if(cache.size() == lastWrittenKey)
    {
    	return;
	}

    for (size_t i = lastWrittenKey; i < cache.size(); ++i)
    {
        file << cache[i] << "\n"; 
        lastWrittenKey = i + 1; 
    }
}

void readVectorFromFile(std::vector<mpfr_float>& cache, const std::string& filename, int64_t &lastWrittenKey)
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Failed to open file for reading. Starting with an empty cache.\n";
        return;
    }

    std::string line;
    int64_t i = 0;
    while (std::getline(file, line) && i <=  lastWrittenKey)
    {
        std::istringstream iss(line);
        std::string part;
        std::getline(iss, part, ','); // Extract key
        int64_t key = std::stoll(part);

        std::getline(iss, part); // Extract value
        mpfr_float value(part);

        if (key > cache.size()) 
        {
            cache.resize(key); // Ensure vector is large enough
        }
        cache[key - 1] = value; // Subtract 1 because vector indices are 0-based
    }

    file.close();
}

void readObjectsFromFile(std::vector<mpfr_float>& vec, const std::string& filename, int64_t& lastWrittenKey)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) 
    {
        std::cerr << "Failed to open file for reading. Starting with an empty cache.\n";
        return;

    }

    mpfr_float obj;
    while (file >> obj)
    {
        vec.push_back(obj);
    }
}
boost::multiprecision::cpp_int newpow_3k(int64_t k);

boost::multiprecision::cpp_int cached_newpow_3k(int64_t k, std::vector<boost::multiprecision::cpp_int> &cache)
{
    // Check if k is within the bounds of the cache vector
    if (k >= 0 && static_cast<size_t>(k) <= cache.size()) 
    {
        // Adjust k to be zero-based for vector access
        if(cache[k] == 0)
        {
            cache[k] = newpow_3k(k);
            return cache[k];
        }
			
        return cache[k]; // Return cached result
    }

    // Calculate result for new k
    boost::multiprecision::cpp_int result = newpow_3k(k);

    // Resize the vector if needed, to accommodate the new k value
    if (static_cast<size_t>(k) > cache.size()) {
        cache.resize(k);
    }

    // Cache the new result (adjust k to be zero-based)
    cache[k] = result;

    return result;
}


mpfr_float cached_result(const int64_t k, std::vector<mpfr_float>& cache)
{
    std::lock_guard lock(cache_mutex); // Ensure thread safety

    // Check if k is within the bounds of the cache vector
    if (k >= 0 && static_cast<size_t>(k) <= cache.size())
    {
        
        if (k >= cache.size())
        {
            cache.push_back( mpfr_float(numerator(k)) / mpfr_float(denominator_a(k) * newpow_3k(k)));
            return cache[k];
        }

        return cache[k]; // Return cached result
    }

    // Cache the new result (adjust k to be zero-based)
    cache.push_back(mpfr_float(numerator(k)) / mpfr_float(denominator_a(k) * newpow_3k(k)));
    return cache[k];
}


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

static cpp_int altpow_3k(int64_t k)
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
    std::vector<mpfr_float> cache;

    if(argc != 2) 
    {
        std::cerr << "Usage: " << argv[0] << " number_of_terms\n";
        return 1;
    }

    mpfr_t result;
    unsigned long int base = 640320;
    unsigned long int exponent = 3;

    // Initialize MPFR variable with a default precision of 53 bits (double precision)
    mpfr_init2(result, 5300);

    // Calculate base^exponent and store the result in 'result'
    mpfr_ui_pow_ui(result, base, exponent, MPFR_RNDN); // RNDN is the rounding mode (to nearest)

    // Print the result
    printf("640320^3 = ");
    mpfr_out_str(stdout, 10, 100000, result, MPFR_RNDN); // Output the result in base 10
    printf("\n");

    // Clear the MPFR variable
    mpfr_clear(result);


    int64_t lastWrittenKey = initializeLastWrittenKey(keyTrackingFile);

    try {
	    const int64_t num_terms = std::stoi(argv[1]);
        
        if (num_terms <= 0) 
        {
            throw std::invalid_argument("Number of terms must be greater than zero.");
        }

        const int precision = calcPrecision(num_terms);
        mpfr_float pi;

        {
            //ExecutionTimer timer;
	        cache.reserve(num_terms);
            mpfr_float::default_precision(precision * log2n);

            const mpfr_float constant = calcConstant(precision);

            readObjectsFromFile(cache, filename, lastWrittenKey);
	        std::cout << "Cache size: " << cache.size() << std::endl;
	                
	        mpfr_float pi_inverse = 0;

            for (int64_t k = 0; k < num_terms; k++)
            {
                const mpfr_float temp = cached_result(k, cache);
                pi_inverse += temp;
            }
			pi = mpfr_float(1) / (pi_inverse * constant);

        }

        std::cout << std::setprecision(precision) << pi << "\n";
        writeObjectsToFile(cache, filename, lastWrittenKey);
        saveLastWrittenKey(keyTrackingFile, lastWrittenKey);

    }

    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}




