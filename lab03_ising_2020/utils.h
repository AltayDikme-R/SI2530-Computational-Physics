// MC simulation of a 2D Ising model
// Andreas Andersson, 2010
// Petter Johansson, Tymofiy Lukinov 2015

#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

using namespace std;

namespace utils {

class UniformRNG {
public:
    // Set-up the seeded generator
    UniformRNG (const int seed, const int L)
    :engine(seed),
     uniform_real(0.0, 1.0),
     L{ L } {}

    // Generate an integer in the range [0, L), ie. a 1D lattice site
    int gen_lattice_site() { return static_cast<int>(uniform_real(engine) * L); }

    // Generate a real in the range [0, 1)
    float gen_real() { return uniform_real(engine); }

private:
    default_random_engine engine; // Generation algorithm
    uniform_int_distribution<int> uniform_int;
    uniform_real_distribution<float> uniform_real;
    int L;
};

// Clear line (or at least 50 chars ... )
static const string clearChars(const int num) noexcept
{
    ostringstream oss;
    oss << '\r' << setw(num) << ' ' << '\r';
    return oss.str();
}
static const string clearLine {clearChars(50)};

static string getTemp(const double T) noexcept
{
    ostringstream oss;
    oss << "temp = "
        << left << fixed << setw(5) << setprecision(1) << T
        << setw(4) << ' ';
    return oss.str();
}

void printEquil(const double T) noexcept
{
    clog << clearLine << getTemp(T) << "equilibrating";
}

void printStep(const int step, const int Nsample, const double T) noexcept
{
    if (step % 5000 == 0)
    {
        clog << clearLine << getTemp(T) << "step " << step << " of " << Nsample;
    }
}

void printParams(const int seed,
                 const int L,
                 const int Nwarmup,
                 const int Nsample,
                 const double Tmin,
                 const double Tmax,
                 const double dT) noexcept
{
    clog << "Parameters:\n"
         << "seed    = " << seed << '\n'
         << "L       = " << L << " (total size is L*L)\n"
         << "Nwarmup = " << Nwarmup << " (equilibration steps)\n"
         << "Nsample = " << Nsample << " (sampling steps)\n"
         << "Tmin    = " << Tmin << '\n'
         << "Tmax    = " << Tmax << '\n'
         << "dT      = " << dT << '\n'
         << '\n';
}

void printUsage(const char* progname) noexcept
{
    cerr << "Usage: " << progname << " SEED L NWARMUP NSAMPLE TMIN TMAX DT FILE\n"
         << '\n'
         << "   SEED is a seed for the random number generator\n"
         << "   L is the system lattice size (total size is L*L)\n"
         << "   NWARMUP is the number of equilibration rounds per temperature\n"
         << "   NSAMPLE is the number of sampling rounds per temperature\n"
         << "   TMIN is theconst string initial simulation temperature\n"
         << "   TMAX is the final simulation temperature\n"
         << "   DT is the temperature difference per step\n"
         << "   FILE is the output file\n";
}

template<typename T>
T parseInput(const char* argv[],
             const int num)
{
    T value;
    istringstream iss(argv[num]);
    iss >> value;

    if (iss.fail())
    {
        cerr << "Could not read argument number " << num
             << " ('" << argv[num] << "'): verify that it is correct.\n"
             << '\n';
        printUsage(argv[0]);
        exit(1);
    }

     return value;
}

ofstream openOutputFile(const string filename) noexcept
{
    clog << "Writing data to file '" << filename << "'.\n";

    const ifstream file(filename);
    if (file.good()) {
        cerr << "Warning: file already exists! Content will be overwritten.\n";
    }

    ofstream outfile(filename, ios::trunc);
    if (outfile.fail()) {
        cerr << "Error: could not open file for writing.\n";
        exit(1);
    }

    clog << '\n';

    return outfile;
}

void outputData(ofstream& out,
                const double T,
                const double m,
                const double mabs,
                const double chi,
                const double e,
                const double c,
                const double binder,
                const int L,
                const int seed)
{
    out << setprecision(3)
        << setw(4) << T << ' '
        << setw(9) << m << ' '
        << setw(9) << mabs << ' '
        << setw(9) << chi << ' '
        << setw(9) << e << ' '
        << setw(9) << c << ' '
        << setw(9) << binder << ' '
        << setw(5) << L << ' '
        << setw(6) << seed << ' '
        << '\n' << flush;
}

} // namespace io
