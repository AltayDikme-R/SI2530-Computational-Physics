// MC simulation of a 2D Ising model
// Andreas Andersson, 2010
// Petter Johansson, Tymofiy Lukinov 2015
// Marcus Pernow, 2017

#include <cmath>
#include <vector>

#include "utils.h"

using namespace std;

constexpr int Lmax = 100;

// Task for the reader: Figure out what this variable does! ->
// This variable is the number of neighbours each point in the lattice has. 2D = 8 3D = 18 (Next to nearest neighbour)
// ****Change the following line for 3D model!
constexpr int dHmax = 18;

// A 2D L-by-L lattice is represented by a 1D vector of length L
// where each element is another 1D vector of length L.
// **** Change the following line for a 3D model.
using Lattice = vector<vector<vector<int>>>;

static void MCupdate(Lattice& lattice,
                     double& m,
                     double& e,
                     const vector<double>& wtable,
                     const vector<int>& plus,
                     const vector<int>& minus,
                     const int L,
                     utils::UniformRNG& rng)
{
    // ****Change the following line for 3D model!
    const auto num_rounds = L * L * L;

    // ****Change the following line for 3D model!
    //Havnt changed this. What do I need to change? It looks ok??
    for (int round = 0; round < num_rounds; ++round)
    {
        const auto i = rng.gen_lattice_site();
        const auto j = rng.gen_lattice_site();
        const auto k = rng.gen_lattice_site();
        // ****Add a line here for 3D model!

        const auto i0 = plus.at(i);
        const auto i1 = minus.at(i);
        const auto j0 = plus.at(j);
        const auto j1 = minus.at(j);
        // 3D line added
        const auto k0 = plus.at(k);
        const auto k1 = minus.at(k);
    	// ****Add two lines here for 3D model!

        // Compute energy difference dH for a flip
    	// ****Change the following line for 3D model!
        const auto dH = 2 * lattice[i][j][k] * (
                  lattice[i0][j][k] + lattice[i1][j][k]
                + lattice[i][j0][k] + lattice[i][j1][k]
                + lattice[i][j][k0] + lattice[i][j][k1]
        );

        // Get transition probability w
        const auto w = wtable.at(dH + dHmax);

        // Test for acceptance
        if (w > rng.gen_real())
        {
    	    // ****Change the following 2 lines for 3D model!
            lattice[i][j][k] = -lattice[i][j][k];
            m += 2 * lattice[i][j][k];
            e += dH;
        }
    }
}

static Lattice init_lattice(const int L)
{
    Lattice lattice;

    // ****Add another loop for a 3D model!
    for (auto i = 0; i < L; ++i)
    {
        // Create a (row) vector of L ones and add it to the lattice.
        const auto row = vector<vector<int>>(L, vector<int>(L, 1));
        lattice.push_back(row);
    }



    return lattice;
}

static void reset_lattice(Lattice& lattice)
{
	// ****Add a nested loop for 3D model!
    for (auto& row : lattice)
    {
        // Assign the value 1 to the entire row to reset it
        row.assign(row.size(), vector<int>(row[0].size(), 1));
    }
}

int main(const int argc, const char* argv[])
{
    if (argc < 9)
    {
        cerr << "Error: not all required parameters were entered.\n\n";
        utils::printUsage(argv[0]);
        exit(1);
    }

    // Set input parameters
    const auto seed        = utils::parseInput<int>(argv, 1);
    const auto L           = utils::parseInput<int>(argv, 2);
    const auto num_warmup  = utils::parseInput<int>(argv, 3);
    const auto num_sample  = utils::parseInput<int>(argv, 4);
    const auto Tmin     = utils::parseInput<double>(argv, 5);
    const auto Tmax     = utils::parseInput<double>(argv, 6);
    const auto dT       = utils::parseInput<double>(argv, 7);
    const auto filename = utils::parseInput<string>(argv, 8);

    // Total area/volume of lattice
    // ****Change the following line for 3D model!
    const auto total_size = L * L * L;

    if (L > Lmax)
    {
        cerr << "Input box size L (" << L << ") "
             << "is larger than maximum allowed (" << Lmax << ").\n";
        exit(1);
    }

    utils::printParams(seed, L, num_warmup, num_sample, Tmin, Tmax, dT);
    auto outfile = utils::openOutputFile(filename);

    // Initiate random number generator
    auto rng = utils::UniformRNG(seed, L);

    // Initialize lattice
    Lattice lattice = init_lattice(L);

    // Declare neighbour vectors
    vector<int> plus;
    vector<int> minus;

    // Initialize nearest neighbour indices
    for (int i = 0; i < L; i++)
    {
        plus.push_back(i + 1);
        minus.push_back(i - 1);
    }

    // Lattice PBC
    minus.front() = L - 1;
    plus.back() = 0;

    // Declare table of acceptance probabilities
    vector<double> wtable;
    wtable.reserve(2 * dHmax + 1);

    for (auto T = Tmin; T <= Tmax + 1e-14; T += dT)
    {
        // Contruct table of acceptance probabilities
        wtable.clear();
        for (auto dH = -dHmax; dH <= dHmax; dH++)
        {
            wtable.push_back(exp(-dH / T));
        }

        // Initialize system with spin up
        reset_lattice(lattice);

        // Initial magnetization and energy
        double dm = total_size;
        double de = -2 * total_size;

        // Equilibrate system
        utils::printEquil(T);
        for (int step = 0; step < num_warmup; step++)
        {
            MCupdate(lattice, dm, de, wtable, plus, minus, L, rng);
        }

        // Initialize expectation values
        auto m     = 0.0;
        auto mabs  = 0.0;
        auto mabs2 = 0.0;
        auto mabs4 = 0.0;
        auto e     = 0.0;
        auto e2    = 0.0;

        // Sample equilibrium data
        for (int step = 0; step < num_sample; step++)
        {
            utils::printStep(step, num_sample, T);
            MCupdate(lattice, dm, de, wtable, plus, minus, L, rng);

            m     += dm;
            mabs  += abs(dm);
            mabs2 += dm * dm;
            mabs4 += dm * dm * dm * dm;
            e     += de;
            e2    += de * de;
        }

        // Calculate expectation values
        m     /= num_sample;
        mabs  /= num_sample;
        mabs2 /= num_sample;
        mabs4 /= num_sample;
        e     /= num_sample;
        e2    /= num_sample;

        // Calculate remaining variables
        auto chi = (mabs2 - mabs * mabs) / (T * total_size);
        auto c = (e2 - e * e) / (T * T * total_size);
        auto binder = 1.0 - mabs4 / (3 * mabs2 * mabs2);

        // Normalize to system size
        m    /= total_size;
        mabs /= total_size;
        e    /= total_size;

        // Output
        utils::outputData(outfile, T, m, mabs, chi, e, c, binder, L, seed);
    }

    clog << '\n';

    return 0;
}
