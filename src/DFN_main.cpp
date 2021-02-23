#include <iostream>
#include "Loop_DFN_WL.h"
#include <sys/time.h>
#include <chrono>

int main()
{
    auto start = std::chrono::steady_clock::now();

    const gsl_rng_type *T;
    gsl_rng *random_seed;
    struct timeval tv;
    gettimeofday(&tv, 0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_default;
    random_seed = gsl_rng_alloc(T);
    gsl_rng_set(random_seed, mySeed);
    //gsl_rng_env_setup();

    //----
    string orientation_distribution;
    string percolation_direction;
    string fracture_size_distriution;
    double ratio_S;      // ratio of tri-ele length to minimum fracture size
    double min_frac_len; // minimum fracture radius
    DFN::Loop_DFN loop_1;

    std::ifstream oii("DFN_WL.inp", std::ios::in);
    if (!oii)
    {
        std::cout << "Please define an input.inp!\n";
        exit(0);
    }

    //---loop parameters
    oii >> ratio_S;
    oii.ignore(300, '\n');
    oii >> min_frac_len;
    oii.ignore(300, '\n');

    oii >> loop_1.Nproc;
    oii.ignore(300, '\n');
    oii >> loop_1.times;
    oii.ignore(300, '\n');
    oii >> loop_1.nt;
    oii.ignore(300, '\n');
    oii >> loop_1.nk;
    oii.ignore(300, '\n');
    oii >> loop_1.nv_MC_TIMES;
    oii.ignore(300, '\n');
    oii >> loop_1.nx;
    oii.ignore(300, '\n');
    //---model size
    oii >> loop_1.L;
    oii.ignore(300, '\n');
    //---num of fracture sets
    oii >> loop_1.NumofFsets;
    oii.ignore(300, '\n');
    //---orientation distribution pattern
    oii >> orientation_distribution;
    oii.ignore(300, '\n');
    oii >> fracture_size_distriution;
    oii.ignore(300, '\n');
    oii >> percolation_direction;
    oii.ignore(300, '\n');
    //---fracture size distribution parameters
    if (orientation_distribution == "uniform")
    {
        if (fracture_size_distriution == "powerlaw")
        {
            loop_1.array12.resize(1);
            oii >> loop_1.array12[0][0];
            oii.ignore(300, '\n');
            oii >> loop_1.array12[0][1];
            oii.ignore(300, '\n');
            oii >> loop_1.array12[0][2];
            oii.ignore(300, '\n');
        }
        else if (fracture_size_distriution == "lognormal")
        {
            loop_1.array12.resize(1);
            oii >> loop_1.array12[0][0];
            oii.ignore(300, '\n');
            oii >> loop_1.array12[0][1];
            oii.ignore(300, '\n');
            oii >> loop_1.array12[0][2];
            oii.ignore(300, '\n');
            oii >> loop_1.array12[0][3];
            oii.ignore(300, '\n');
        }
        else if (fracture_size_distriution == "uniform")
        {
            loop_1.array12.resize(1);
            oii >> loop_1.array12[0][0];
            oii.ignore(300, '\n');
            oii >> loop_1.array12[0][1];
            oii.ignore(300, '\n');
        }
        else if (fracture_size_distriution == "single")
        {
            loop_1.array12.resize(1);
            oii >> loop_1.array12[0][0];
            oii.ignore(300, '\n');
        }
    }
    else if (orientation_distribution == "fisher")
    {
        loop_1.array12.resize(loop_1.NumofFsets);
        if (fracture_size_distriution == "powerlaw")
        {
            for (size_t i = 0; i < loop_1.NumofFsets; ++i)
            {
                oii >> loop_1.array12[i][0];
                oii.ignore(300, '\n');
                oii >> loop_1.array12[i][1];
                oii.ignore(300, '\n');
                oii >> loop_1.array12[i][2];
                oii.ignore(300, '\n');
            }
        }
        else if (fracture_size_distriution == "lognormal")
        {
            for (size_t i = 0; i < loop_1.NumofFsets; ++i)
            {
                oii >> loop_1.array12[i][0];
                oii.ignore(300, '\n');
                oii >> loop_1.array12[i][1];
                oii.ignore(300, '\n');
                oii >> loop_1.array12[i][2];
                oii.ignore(300, '\n');
                oii >> loop_1.array12[i][3];
                oii.ignore(300, '\n');
            }
        }
        else if (fracture_size_distriution == "uniform")
        {
            for (size_t i = 0; i < loop_1.NumofFsets; ++i)
            {
                oii >> loop_1.array12[i][0];
                oii.ignore(300, '\n');
                oii >> loop_1.array12[i][1];
                oii.ignore(300, '\n');
            }
        }
        else if (fracture_size_distriution == "single")
        {
            for (size_t i = 0; i < loop_1.NumofFsets; ++i)
            {
                oii >> loop_1.array12[i][0];
                oii.ignore(300, '\n');
            }
        }
    }
    //---orientation distribution parameters
    if (orientation_distribution == "fisher")
    {
        loop_1.array13.resize(loop_1.NumofFsets);
        loop_1.DenWeight.resize(loop_1.NumofFsets);
        for (size_t i = 0; i < loop_1.NumofFsets; ++i)
        {
            oii >> loop_1.DenWeight[i];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][0];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][1];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][2];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][3];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][4];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][5];
            oii.ignore(300, '\n');
            oii >> loop_1.array13[i][6];
            oii.ignore(300, '\n');
        }
    }
    oii.close();
    double avg_ele_len = ratio_S * pow(2 * pow(min_frac_len, 2), 0.5);
    loop_1.Loop_create_DFNs(random_seed, orientation_distribution, fracture_size_distriution, percolation_direction, avg_ele_len);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::micro> elapsed = end - start; // std::micro time (us)
    std::cout << "Running time: " << (((double)(elapsed.count() * 1.0) * (0.000001)) / 60.00) / 60.00 << "h" << std::endl;
    gsl_rng_free(random_seed);
    return 0;
}