/* Copyright (c) 2018 M. Grady Saunders
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer in the documentation and/or other materials
 *      provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// for std::exp, std::log, ...
#include <cmath>

// for std::exit
#include <cstdlib>

// for std::strcmp, std::sscanf
#include <cstring>

// for std::swap
#include <algorithm>

// for std::cout, std::cerr, ...
#include <iostream>

// for std::pair
#include <utility>

// for std::generate_canonical
#include <random>

// for pcg32_k64
#include <pcg_random.hpp>

// minimum ground sample distance in inches
static double GSDmin = 3.93701e-1;

// maximum ground sample distance in inches
static double GSDmax = 3.93701e+2/5;

// minimum relative edge response
static double RERmin = 0.15;

// maximum relative edge response
static double RERmax = 0.95;

// minimum noise gain due to sharpening
static double Gmin = 1.0;

// maximum noise gain due to sharpening
static double Gmax = 50.0;

// minimum edge overshoot due to sharpening
static double Hmin = 1.0;

// maximum edge overshoot due to sharpening
static double Hmax = 2.0;

// minimum signal to noise ratio
static double SNRmin = 1.0;

// maximum signal to noise ratio
static double SNRmax = 100.0;

// giqev4
double g4(const double* X)
{
    double GSD = (1 - X[0]) * GSDmin + X[0] * GSDmax;
    double RER = (1 - X[1]) * RERmin + X[1] * RERmax;
    double G = (1 - X[2]) * Gmin + X[2] * Gmax;
    double H = (1 - X[3]) * Hmin + X[3] * Hmax;
    double SNR = (1 - X[4]) * SNRmin + X[4] * SNRmax;
    return 
        10.251 - 
        (RER >= 0.9 ? 3.320 : 3.160) * std::log10(GSD) + 
        (RER >= 0.9 ? 1.559 : 2.817) * std::log10(RER) - 
        0.656 * H - 
        0.344 * G / SNR;
}

double step(double x)
{
    if (x < 0) {
        return 0;
    }
    else {
        return 1;
    }
}

// integral of GSD partial
double int_dg4_dGSD(double GSD, double RER)
{
    return -std::log10(GSD) * 
        ((0.16 * RER - 0.144) * step(RER - 0.9) + 3.16 * RER);
}

// integral of RER partial
double int_dg4_dRER(double RER, double GSD)
{
    return GSD / std::log(10) * 
                (step(RER - 0.9) * (0.16 * 
                (1 - std::log(GSD)) - 1.258 * std::log(RER)) + 
                                      2.817 * std::log(RER));
}

// integral of G partial
double int_dg4_dG(double G, double SNR)
{
    return -0.344 * G * std::log(SNR);
}

// integral of H partial
double int_dg4_dH(double H)
{
    return -0.656 * H;
}

// integral of SNR partial
double int_dg4_dSNR(double SNR, double G)
{
    return -0.172 * G * G / SNR;
}

// giqev5
double g5(const double* X)
{
    double GSD = (1 - X[0]) * GSDmin + X[0] * GSDmax;
    double RER = (1 - X[1]) * RERmin + X[1] * RERmax;
    double SNR = (1 - X[2]) * SNRmin + X[2] * SNRmax;
    return 
        9.57 - 
        3.32 * std::log10(GSD) - 
        3.32 * std::expm1(-1.9 / SNR) * std::log10(RER) -
        2.0 * std::pow(std::log10(RER), 4.0) -
        1.8 / SNR;
}

// integral of GSD partial
double int_dg5_dGSD(double GSD)
{
    return -3.32 * std::log10(GSD);
}

// integral of RER partial
double int_dg5_dRER(double RER, double SNR)
{
    return -3.32 * std::log10(RER) * 
        (-1.9 * std::expint(-1.9 / SNR) + 
                SNR * (1 - std::exp(-1.9 / SNR))) - 
            2 * SNR * std::pow(std::log10(RER), 4);
}

// integral of SNR partial
double int_dg5_dSNR(double SNR, double RER)
{
    return RER * (-3.32 / std::log(10) * std::exp(-1.9 / SNR) * 
        (std::log(RER) - 1) - 1 / SNR);
}

void genSamples(pcg32_k64& gen, int N, int d, double* X)
{
    // stratify
    for (int k = 0; k < N; k++)
    for (int j = 0; j < d; j++) {
        X[k * d + j] = (k + std::generate_canonical<double, 53>(gen)) / N;
    }

    // shuffle per column
    for (int k = N - 1; k > 0; k--)
    for (int j = 0; j < d; j++) {
        std::swap(X[k * d + j], X[gen(k + 1) * d + j]);
    }
}

void run(
        int seed,
        int stream,
        int N, 
        int d, 
        double (*f)(const double*),
        double* m1, 
        double* m2, 
        double* V,
        double* E)
{
    // zero outputs
    *m1 = 0;
    *m2 = 0;
    for (int i = 0; i < d; i++) {
        V[i] = 0;
        E[i] = 0;
    }

    // random generator
    pcg32_k64 gen(seed, stream);

    // allocate matrices
    double* A = new double[N * d];
    double* B = new double[N * d];
    double* fA = new double[N];
    double* fB = new double[N];
    double** C = new double*[d];
    double** fC = new double*[d];
    for (int i = 0; i < d; i++) {
        C[i] = new double[N * d];
        fC[i] = new double[N];
    }

    // sample and evaluate
    genSamples(gen, N, d, A);
    genSamples(gen, N, d, B);
    for (int k = 0; k < N; k++) {
        fA[k] = f(&A[k * d]);
        fB[k] = f(&B[k * d]);
    }
    for (int i = 0; i < d; i++) {
        for (int k = 0; k < N; k++)
        for (int j = 0; j < d; j++) {
            C[i][k * d + j] = (j == i) ? B[k * d + j] : A[k * d + j];
        }
        for (int k = 0; k < N; k++) {
            fC[i][k] = f(&C[i][k * d]);
        }
    }

    // estimate variances
    for (int i = 0; i < d; i++)
    for (int k = 0; k < N; k++) {
        V[i] += fB[k] * (fC[i][k] - fA[k]) / N;
        E[i] += (fA[k] - fC[i][k]) * (fA[k] - fC[i][k]) / (2 * N);
    }

    // estimate 1st moment
    for (int k = 0; k < N; k++) {
        *m1 += fA[k] / (2 * N);
        *m1 += fB[k] / (2 * N);
    }

    // estimate 2nd moment
    for (int k = 0; k < N; k++) {
        *m2 += (fA[k] - *m1) * (fA[k] - *m1) / (2 * N - 1);
        *m2 += (fB[k] - *m1) * (fB[k] - *m1) / (2 * N - 1);
    }

    // clean up
    delete[] A;
    delete[] B;
    delete[] fA;
    delete[] fB;
    for (int i = 0; i < d; i++) {
        delete[] C[i];
        delete[] fC[i];
    }
    delete[] C;
    delete[] fC;
}

int main(int argc, char** argv)
{
    int seed = 0xABCDE;
    int Niters = 128;
    int Nsamps = 4096;

    // options
    static std::pair<const char*, int*> opt_int[] = {
        {"--seed", &seed},
        {"--Nsamps", &Nsamps},
        {"--Niters", &Niters}
    };
    static std::pair<const char*, double*> opt_double[] = {
        {"--GSDmin", &GSDmin},
        {"--GSDmax", &GSDmax},
        {"--RERmin", &RERmin},
        {"--RERmax", &RERmax},
        {"--Gmin", &Gmin},
        {"--Gmax", &Gmax},
        {"--Hmin", &Hmin},
        {"--Hmax", &Hmax},
        {"--SNRmin", &SNRmin},
        {"--SNRmax", &SNRmax}
    };
    const char* prog = *argv;
    --argc;
    ++argv;
    while (argc-- > 0) {
        for (auto& opt : opt_int) {
            if (std::strcmp(argv[0], opt.first) == 0) {
                if (argc == 0 || 
                    std::sscanf(argv[1], "%d", opt.second) != 1) {
                    // error
                    std::cerr << "Error parsing " << opt.first << ", ";
                    std::cerr << "expected int.\n";
                    std::exit(1);
                }
                // consume
                --argc;
                ++argv;
                continue;
            }
        }
        for (auto& opt : opt_double) {
            if (std::strcmp(argv[0], opt.first) == 0) {
                if (argc == 0 || 
                    std::sscanf(argv[1], "%lf", opt.second) != 1) {
                    // error
                    std::cerr << "Error parsing " << opt.first << ", ";
                    std::cerr << "expected double.\n";
                    std::exit(1);
                }
                // consume
                --argc;
                ++argv;
                continue;
            }
        }

        if (std::strcmp(argv[0], "-h") == 0 ||
            std::strcmp(argv[0], "--help") == 0) {

            // print help
            std::cout << "Usage: " << prog << " [OPTIONS]\n\n";
            std::cout << "--seed INT\n";
            std::cout << "\tSet random seed.\n\n";
            std::cout << "--Nsamps INT\n";
            std::cout << "\tSet number of samples per iteration.\n";
            std::cout << "\tBy default, 4096.\n\n";
            std::cout << "--Niters INT\n";
            std::cout << "\tSet number of iterations.\n";
            std::cout << "\tBy default, 128.\n\n";
            std::cout << "--GSDmin DOUBLE\n";
            std::cout << "\tSet minimum ground sample distance in inches.\n";
            std::cout << "\tBy default, 3.93701e-1 (1 centimeter).\n\n";
            std::cout << "--GSDmax DOUBLE\n";
            std::cout << "\tSet maximum ground sample distance in inches.\n";
            std::cout << "\tBy default, 3.93701e+2/2 (5 meters).\n\n";
            std::cout << "--RERmin DOUBLE\n";
            std::cout << "\tSet minimum relative edge response.\n";
            std::cout << "\tBy default, 0.15.\n\n";
            std::cout << "--RERmax DOUBLE\n";
            std::cout << "\tSet maximum relative edge response.\n";
            std::cout << "\tBy default, 0.95.\n\n";
            std::cout << "--Gmin DOUBLE\n";
            std::cout << "\tSet minimum noise gain due to sharpening.\n";
            std::cout << "\tBy default, 1.\n\n";
            std::cout << "--Gmax DOUBLE\n";
            std::cout << "\tSet maximum noise gain due to sharpening.\n";
            std::cout << "\tBy default, 50.\n\n";
            std::cout << "--Hmin DOUBLE\n";
            std::cout << "\tSet minimum edge overshoot due to sharpening.\n";
            std::cout << "\tBy default, 1.\n\n";
            std::cout << "--Hmax DOUBLE\n";
            std::cout << "\tSet maximum edge overshoot due to sharpening.\n";
            std::cout << "\tBy default, 2.\n\n";
            std::cout << "--SNRmin DOUBLE\n";
            std::cout << "\tSet minimum signal to noise ratio.\n";
            std::cout << "\tBy default, 1.\n\n";
            std::cout << "--SNRmax DOUBLE\n";
            std::cout << "\tSet maximum signal to noise ratio.\n";
            std::cout << "\tBy default, 100.\n\n";

            // exit
            std::exit(0);
        }

        // error
        std::cerr << "Error parsing " << *argv << ", unknown option.\n";
        std::exit(1);
    }

    // print
    std::cout << "input:\n\n";
    std::cout << "\tseed = " << seed << "\n";
    std::cout << "\tNiters = " << Niters << "\n";
    std::cout << "\tNsamps = " << Nsamps << "\n";
    std::cout << "\tGSDmin = " << GSDmin << "\n";
    std::cout << "\tGSDmax = " << GSDmax << "\n";
    std::cout << "\tRERmin = " << RERmin << "\n";
    std::cout << "\tRERmax = " << RERmax << "\n";
    std::cout << "\tGmin = " << Gmin << "\n";
    std::cout << "\tGmax = " << Gmax << "\n";
    std::cout << "\tHmin = " << Hmin << "\n";
    std::cout << "\tHmax = " << Hmax << "\n";
    std::cout << "\tSNRmin = " << SNRmin << "\n";
    std::cout << "\tSNRmax = " << SNRmax << "\n\n";

    {
        // run
        double m1 = 0;
        double m2 = 0;
        double V[5] = {};
        double E[5] = {};
        for (int iter = 0; iter < Niters; iter++) {
            double tmpm1;
            double tmpm2;
            double tmpV[5];
            double tmpE[5];
            run(
                seed, 
                iter, 
                Nsamps, 5, 
                g4, 
                &tmpm1, 
                &tmpm2, 
                &tmpV[0], 
                &tmpE[0]);
            m1 += tmpm1 / Niters;
            m2 += tmpm2 / Niters;
            for (int i = 0; i < 5; i++) {
                V[i] += tmpV[i] / Niters;
                E[i] += tmpE[i] / Niters;
            }
        }

        // deriv expectations
        double E_dg4_dGSD = 
            (int_dg4_dGSD(GSDmax, RERmax) - 
             int_dg4_dGSD(GSDmax, RERmin) - 
             int_dg4_dGSD(GSDmin, RERmax) +
             int_dg4_dGSD(GSDmin, RERmin)) /
            (GSDmax - GSDmin) / (RERmax - RERmin);
        double E_dg4_dRER = 
            (int_dg4_dRER(RERmax, GSDmax) - 
             int_dg4_dRER(RERmax, GSDmin) - 
             int_dg4_dRER(RERmin, GSDmax) +
             int_dg4_dRER(RERmin, GSDmin)) /
            (RERmax - RERmin) / (GSDmax - GSDmin);
        double E_dg4_dG = 
            (int_dg4_dG(Gmax, SNRmax) - 
             int_dg4_dG(Gmax, SNRmin) - 
             int_dg4_dG(Gmin, SNRmax) +
             int_dg4_dG(Gmin, SNRmin)) /
            (Gmax - Gmin) / (SNRmax - SNRmin);
        double E_dg4_dH = 
            (int_dg4_dH(Hmax) - 
             int_dg4_dH(Hmin)) /
            (Hmax - Hmin);
        double E_dg4_dSNR = 
            (int_dg4_dSNR(SNRmax, Gmax) - 
             int_dg4_dSNR(SNRmax, Gmin) - 
             int_dg4_dSNR(SNRmin, Gmax) +
             int_dg4_dSNR(SNRmin, Gmin)) /
            (SNRmax - SNRmin) / (Gmax - Gmin);

        // print
        std::cout << "giqev4 output:\n\n";
        std::cout << "\tm1 = " << m1 << "\n";
        std::cout << "\tm2 = " << m2 << "\n\n";
        std::cout << "\tmain effects:\n";
        std::cout << "\t\tS(GSD) = " << V[0] / m2 << "\n";
        std::cout << "\t\tS(RER) = " << V[1] / m2 << "\n";
        std::cout << "\t\tS(G) = " << V[2] / m2 << "\n";
        std::cout << "\t\tS(H) = " << V[3] / m2 << "\n";
        std::cout << "\t\tS(SNR) = " << V[4] / m2 << "\n\n";
        std::cout << "\ttotal effects:\n";
        std::cout << "\t\tST(GSD) = " << E[0] / m2 << "\n";
        std::cout << "\t\tST(RER) = " << E[1] / m2 << "\n";
        std::cout << "\t\tST(G) = " << E[2] / m2 << "\n";
        std::cout << "\t\tST(H) = " << E[3] / m2 << "\n";
        std::cout << "\t\tST(SNR) = " << E[4] / m2 << "\n\n";
        std::cout << "\tderivatives:\n";
        std::cout << "\t\tE(dg4/dGSD) = " << E_dg4_dGSD << "\n";
        std::cout << "\t\tE(dg4/dRER) = " << E_dg4_dRER << "\n";
        std::cout << "\t\tE(dg4/dG) = " << E_dg4_dG << "\n";
        std::cout << "\t\tE(dg4/dH) = " << E_dg4_dH << "\n";
        std::cout << "\t\tE(dg4/dSNR) = " << E_dg4_dSNR << "\n\n";
    }
    {
        // run
        double m1 = 0;
        double m2 = 0;
        double V[3] = {};
        double E[3] = {};
        for (int iter = 0; iter < Niters; iter++) {
            double tmpm1;
            double tmpm2;
            double tmpV[3];
            double tmpE[3];
            run(
                seed, 
                iter, 
                Nsamps, 3, 
                g5, 
                &tmpm1, 
                &tmpm2, 
                &tmpV[0], 
                &tmpE[0]);
            m1 += tmpm1 / Niters;
            m2 += tmpm2 / Niters;
            for (int i = 0; i < 3; i++) {
                V[i] += tmpV[i] / Niters;
                E[i] += tmpE[i] / Niters;
            }
        }

        // deriv expectations
        double E_dg5_dGSD = 
            (int_dg5_dGSD(GSDmax) - 
             int_dg5_dGSD(GSDmin)) /
            (GSDmax - GSDmin);
        double E_dg5_dRER = 
            (int_dg5_dRER(RERmax, SNRmax) - 
             int_dg5_dRER(RERmax, SNRmin) - 
             int_dg5_dRER(RERmin, SNRmax) +
             int_dg5_dRER(RERmin, SNRmin)) /
            (RERmax - RERmin) / (SNRmax - SNRmin);
        double E_dg5_dSNR = 
            (int_dg5_dSNR(SNRmax, RERmax) - 
             int_dg5_dSNR(SNRmax, RERmin) - 
             int_dg5_dSNR(SNRmin, RERmax) +
             int_dg5_dSNR(SNRmin, RERmin)) /
            (SNRmax - SNRmin) / (RERmax - RERmin);

        // print
        std::cout << "giqev5 output:\n\n";
        std::cout << "\tm1 = " << m1 << "\n";
        std::cout << "\tm2 = " << m2 << "\n\n";
        std::cout << "\tmain effects:\n";
        std::cout << "\t\tS(GSD) = " << V[0] / m2 << "\n";
        std::cout << "\t\tS(RER) = " << V[1] / m2 << "\n";
        std::cout << "\t\tS(SNR) = " << V[2] / m2 << "\n\n";
        std::cout << "\ttotal effects:\n";
        std::cout << "\t\tST(GSD) = " << E[0] / m2 << "\n";
        std::cout << "\t\tST(RER) = " << E[1] / m2 << "\n";
        std::cout << "\t\tST(SNR) = " << E[2] / m2 << "\n\n";
        std::cout << "\tderivatives:\n";
        std::cout << "\t\tE(dg5/dGSD) = " << E_dg5_dGSD << "\n";
        std::cout << "\t\tE(dg5/dRER) = " << E_dg5_dRER << "\n";
        std::cout << "\t\tE(dg5/dSNR) = " << E_dg5_dSNR << "\n\n";
    }
    return 0;
}
