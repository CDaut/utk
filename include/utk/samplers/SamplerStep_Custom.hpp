/*
 * Coded by Hélène Perrier helene.perrier@liris.cnrs.fr
 * and Bastien Doignies bastien.doignies@liris.cnrs.fr 
 * and David Coeurjolly David.coeurjolly@liris.cnrs.fr
 *
 * Copyright (c) 2023 CNRS Université de Lyon
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the UTK project.
 */
#pragma once

#include "../../../externals/Step/common.hh"
#include "utk/metrics/RadialSpectrum.hpp"
#include "utk/utils/PointsetIO.hpp"

#include <cmath>
#include <iomanip>
#include <random>
#include <utk/utils/FastPRNG.hpp>
#include <utk/utils/Pointset.hpp>


#ifndef STEP_COMMON_IMPLEMENTATION
#define STEP_COMMON_IMPLEMENTATION

#ifdef _OPENMP
#include <omp.h>
#include <optional>
#endif

namespace heckOptimized {

    struct PointWithDistance {
        float_t contrib;
        float_t distance;

        explicit PointWithDistance(float_t contrib = 0.f, float_t dist = -1.f) : contrib(contrib), distance(dist) {}
    };


    int nbins = -1;// default: number of points
    bool randomize_force = false;
    int maxattempts = 5;
    float Tmax = 1, Tmin = 1e-3, Tfactor = 0.9;

    // A flag set on SIGKILL to stop optimization
    bool abort_flag = false;


    template<typename T>
    void writeRadspecToFile(const std::string &filename,
                            std::pair<std::vector<T>, std::vector<T>> radspec) {
        std::ofstream file;
        file.open(filename);

        auto xs = radspec.first;
        auto ys = radspec.second;

        if (xs.size() != ys.size()) {
            std::cerr << "Dimensions of radial spactrum are unequal: xDim: "
                      << xs.size() << " yDim: " << ys.size() << std::endl;
            std::terminate();
        }

        for (int i = 0; i < xs.size(); ++i) {
            file << std::setprecision(std::numeric_limits<long double>::digits10 + 2)
                 << std::fixed;
            file << xs[i] << ", " << ys[i] << std::endl;
        }
    }

    float CalcGradients(const std::vector<heck_Point> &pts,
                        const Curve &rdf, const Curve &target,
                        std::vector<heck_Point> &gradient,
                        int index,
                        float r_tailcorrect,
                        bool tailcorrect,
                        const Curve &tailcurve) {
        Curve force(nbins, 0, 0.5f);

        //initialize fore curve
        //Equation 9??
        for (int j = 0; j < force.size(); j++) {
            force[j] = force.dx * (rdf[j] - target[j]);
        }

        //sum over all entries in force (basiclly integral)
        for (int j = 1; j < force.size(); j++) {
            force[j] += force[j - 1];
        }

        for (int j = 1; j < force.size(); j++) {
            float x = force.ToX(j);
            force[j] /= (pts.size() * x * x);
            // if (force[j] < 0) force[j] = 0;
        }
        //force.Write(fmt::format("force{}.dat", index));

        float maxforce = 0.0f;
        gradient.resize(pts.size());

        //remember how much each point contributed towards the overall force
        auto pointContributions = std::vector<std::vector<PointWithDistance>>{};
        pointContributions.resize(pts.size());

        //#ifdef _OPENMP
        //#pragma omp parallel
        //#endif
        {
            //#ifdef _OPENMP
            //#pragma omp for schedule(static)
            //#endif
            for (uint32_t i = 0; i < pts.size(); i++) {
                heck_Point grad;

                //initialize histogramm
                //pointContributions[i] = std::vector<PointWithDistance>{};

                for (uint32_t j = 0; j < pts.size(); j++) {
                    if (i == j)
                        continue;

                    float dx = pts[j].x - pts[i].x;
                    float dy = pts[j].y - pts[i].y;
                    //toroidal distance if you interpret bool as 0 or 1
                    dx += (dx < -0.5f) - (dx > 0.5f);
                    dy += (dy < -0.5f) - (dy > 0.5f);
                    float dist2 = dx * dx + dy * dy;

                    float f;
                    // RDF is only reliable up to 0.5, if we go higher, the periodic
                    // boundary conditions cause anisotropies in the final spectrum
                    if (dist2 > 0.49f * 0.49f || dist2 == 0) continue;
                    float dist = sqrtf(dist2);

                    //tail correction or cutoff
                    if (tailcorrect && dist > r_tailcorrect) {
                        //use tail correction curve
                        f = tailcurve.At(dist);
                    } else {
                        //use force
                        f = force.At(dist);// / dist;
                    }
                    grad.x -= dx * f;
                    grad.y -= dy * f;

                    // put the strength (2-Norm) contributing point into the histogram
                    //pointContributions[i].emplace_back(std::sqrt((dx * f * dx * f) + (dy * f * dy * f)), dist);
                }
                float ff = grad.x * grad.x + grad.y * grad.y;
                if (ff > maxforce)
                    maxforce = ff;
                gradient[i] = grad;

                //sort the point contributions by distance
//                std::sort(pointContributions[i].begin(), pointContributions[i].end(), [](auto a, auto b) {
//                    return a.distance < b.distance;
//                });
            }
        }

        ///////////////////////////
//        //dump the whole contribution matrix
//        std::ofstream contributions;
//        contributions.open(fmt::format("contribution_histogram{}.hst", index));
//
//        int i = 0;
//        for (const auto &contrib: pointContributions) {
//            contributions << i++ << std::endl;
//            for (const auto &pointContrib: contrib) {
//                contributions << pointContrib.distance << " " << pointContrib.contrib << std::endl;
//            }
//        }
//        contributions.close();
        ///////////////////////////

        return sqrtf(maxforce);
    }

    std::vector<heck_Point> MovePoints(const std::vector<heck_Point> &pts,
                                       const std::vector<heck_Point> &gradient, float stepsize) {
        std::vector<heck_Point> newpts(pts.size());
        for (unsigned i = 0; i < pts.size(); i++) {
            heck_Point force = gradient[i];
            if (randomize_force) {
                force.x *= drand48();
                force.y *= drand48();
            }
            float x = pts[i].x + stepsize * force.x;
            float y = pts[i].y + stepsize * force.y;

            while (x >= 1) x -= 1;
            while (y >= 1) y -= 1;
            while (x < 0) x += 1;
            while (y < 0) y += 1;

            newpts[i].x = x;
            newpts[i].y = y;
        }
        return newpts;
    }


    float CalcEnergy(const Curve &a, const Curve &b) {
        Curve c(a);
        for (int i = 0; i < c.size(); i++)
            c[i] = powf(c[i] - b[i], 2.0f);
        return sqrtf(Integrate(c));
        //    fprintf(stderr, "e = %g\n", x);
        //    return x;
    }

    float MainOptimization(const std::vector<heck_Point> &pts,
                           Curve &target,
                           float smoothing,
                           std::vector<heck_Point> &output,
                           float maxdist,
                           const Curve &tailcurve,
                           const bool tailcorrect) {

        std::vector<heck_Point> current = pts, best = pts;
        std::vector<heck_Point> gradient;
        Curve rdf = CalcRDF(nbins, current.size(), &current[0].x, smoothing);
        float bestenergy = CalcEnergy(rdf, target);
        float T = Tmax;// temperature
        int attempts = 0;
        int i = 0;
        while (!abort_flag && T >= Tmin) {
            // Calculate gradients and move points
            float maxgrad = CalcGradients(current,
                                          rdf,
                                          target,
                                          gradient,
                                          i,
                                          maxdist,
                                          tailcorrect,
                                          tailcurve);

            /////////////////////////////////
            //            std::ofstream forces;
            //            forces.open(fmt::format("forces_{}.txt", i));
            //            //dump gradients
            //            for (size_t index = 0; index < gradient.size(); index++) {
            //                heck_Point grad_i = gradient[index];
            //                heck_Point point_i = current[index];
            //                forces << point_i.x << " " << point_i.y << " " << grad_i.x << " " << grad_i.y << std::endl;
            //            }
            //
            //            forces.close();
            //
            //            std::ofstream pointset;
            //            pointset.open(fmt::format("pointset_{}.txt", i));
            //            for (const auto &item: current) {
            //                pointset << item.x << " " << item.y << std::endl;
            //            }
            //            pointset.close();
            /////////////////////////////////////////////

            float stepsize = T / (sqrt(current.size()) * maxgrad);
            current = MovePoints(current, gradient, stepsize);

            // Perform iteration control
            rdf = CalcRDF(nbins, current.size(), &current[0].x, smoothing);
            float energy = CalcEnergy(rdf, target);
            attempts++;
            if (energy < bestenergy) {
                attempts = 0;
                best = current;
                bestenergy = energy;
            } else if (energy > bestenergy * 1.2) {
                attempts = maxattempts;
                current = best;
            }

//            std::cout << "T: " << T << " stepsize: " << stepsize << " energy: " << energy << " bestenergy: "
//                      << bestenergy << std::endl;

            if (attempts >= maxattempts) {
                attempts = 0;
                T *= Tfactor;
            }
            ++i;
        }
        output = best;
        return bestenergy;
    }


    void FunctionStep(float critx, Curve &c, int npts) {
        float x0 = critx * sqrt(2 / (sqrt(3) * npts));
        for (int i = 0; i < c.size(); i++) {
            c[i] = (c.ToX(i) >= x0) ? 1 : 0;
        }
    }

    void FunctionJinc(float critx, Curve &c, int npts) {
        float x0 = critx / sqrt(2 / (sqrt(3) * npts));
        for (int i = 0; i < c.size(); i++) {
            float xx0 = c.ToX(i);

            float a = 0.0;
            const int N = 20;
            for (int j = 0; j < N; j++) {
                float x = xx0 + j * c.dx / N;
                if (x < 1e-5) a += 1 - M_PI * x0 * x0 / npts;
                else {
                    float xx = 2 * M_PI * x0 * x;
                    a += 1 - 2 * M_PI * x0 * x0 / npts * jn(1, xx) / xx;
                }
            }
            c[i] = a / N;
        }
    }

    void FunctionPeak(float critx, float peaky, Curve &c, int npts) {
        for (int i = 0; i < c.size(); i++) {
            c[i] = (c.ToX(i) > critx) ? 1 : 0;
            c[i] = ((int) c.ToX(i) == (int) critx) ? peaky : c[i];
        }
    }

}// namespace heckOptimized

#endif//  STEP_COMMON_IMPLEMENTATION

namespace utk {

    class CustomHeckSampler {
    protected:
    public:
        explicit CustomHeckSampler(float critFreq = 0.606f,
                                   float smooth = 8.f,
                                   float maxdist = 1.f,
                                   bool tailcorrect = false,
                                   const Curve &mediancurve = Curve{},
                                   bool hardCutoff = false) : critFrequency(critFreq),
                                                              smoothing(smooth),
                                                              maxdist(maxdist),
                                                              do_tailcorrection(tailcorrect),
                                                              mediancurve(mediancurve),
                                                              hardCutoff(hardCutoff) { setRandomSeed(); }

        [[nodiscard]] uint32_t GetDimension() const { return 2; }

        void setRandomSeed(uint64_t arg_seed) { gen.seed(arg_seed); }

        void setRandomSeed() { gen.seed(std::random_device{}()); }

        void setCriticalFrequency(float cf) { critFrequency = cf; }

        void setSmoothing(float sm) { smoothing = sm; }

        template<typename T>
        bool generateSamples(Pointset<T> &arg_pts, uint32_t N) {
            // Note : heck_Points are of type float
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            std::vector<heck_Point> pts, result;
            pts.resize(N);

            for (uint32_t i = 0; i < N; i++) {
                pts[i].x = dist(gen);
                pts[i].y = dist(gen);
            }

            heckOptimized::nbins = (int) pts.size();

            Curve target(heckOptimized::nbins, 0, 0.5f);
            heckOptimized::FunctionJinc(critFrequency, target, (int) pts.size());

            //initialize the tailcurve with all zeroes (this is a hard cutoff)
            const int nbinsTailcurve = 512;
            Curve tailcurve(nbinsTailcurve, maxdist, 0.5);

            if (!hardCutoff && do_tailcorrection) {
                const float zeroPoint = 0.25;
                const float yAtTailBegin = mediancurve.At(maxdist);
                const float linear_correction_m = (-yAtTailBegin) / (zeroPoint - maxdist);
                const int nBinsTailNonzero = std::floor((zeroPoint - maxdist) / (0.5 - maxdist) * nbinsTailcurve);
                const float tailcurveDxUntilZero = (zeroPoint - maxdist) / nBinsTailNonzero;

                //fill tailcurve linearly
                for (int i = 0; i < tailcurve.ToIndex(zeroPoint); ++i) {
                    tailcurve[i] = linear_correction_m * i * tailcurveDxUntilZero + yAtTailBegin;
                }
            }

            //tailcurve.Write(fmt::format("tailcurve_{}.crv", maxdist));

            result = pts;
            heckOptimized::MainOptimization(pts, target, smoothing, result, maxdist, tailcurve, do_tailcorrection);

            arg_pts.Resize(N, 2);
            for (uint32_t i = 0; i < N; i++) {
                arg_pts[i][0] = static_cast<T>(result[i].x);
                arg_pts[i][1] = static_cast<T>(result[i].y);
            }

            return true;
        };

    protected:
        utk::PCG32 gen;
        float critFrequency;
        float smoothing;
        float maxdist;
        bool do_tailcorrection;
        const Curve &mediancurve;
        const bool hardCutoff;
    };

}// namespace utk
