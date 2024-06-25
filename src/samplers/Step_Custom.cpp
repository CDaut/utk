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
#include <utk/utils/SamplerArgumentParser.hpp>
#include <utk/utils/PointsetIO.hpp>
#include <utk/utils/Pointset.hpp>
#include <utk/samplers/SamplerStep_Custom.hpp>
#include <format>
#include "utk/metrics/RDF_Heck.hpp"

template<typename T>
void writeRDFtoFile(const std::string &filename,
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

int main(int argc, char **argv) {
    CLI::App app{"Step sampler"};
    auto *args = utk::add_arguments(app, 2);

    float criticalFrequency = 0.606f;
    float smoothing = 8.f;
    app.add_option("--criticalFreq", criticalFrequency, "Critical frequency")->default_val(criticalFrequency);
    app.add_option("--smoothing", smoothing, "Smoothing")->default_val(smoothing);

    CLI11_PARSE(app, argc, argv);

    for (float maxdist = 0.166016; maxdist <= 0.5; maxdist += 0.0009765625) {
        utk::Pointset<double> pointset;
        utk::CustomHeckSampler sampler(criticalFrequency, smoothing, maxdist);
        sampler.setRandomSeed(args->seed);

        std::cout << "generating pointset with " << args->N << " points and maxdist " << maxdist << "…" << std::endl;
        std::cout << "progress: " << (maxdist / 0.5) * 100 << "%" << std::endl;

        if (!sampler.generateSamples(pointset, args->N)) {
            std::cerr << "Sampler returned non-zero output" << std::endl; // No log here, must be visible whatsoever
            return 1;
        }

        write_text_pointset(fmt::format("pointset_{}.txt", maxdist), pointset);

        std::cout << "generating rdf…" << std::endl;

        auto rdf = utk::RDF{true, 4096, 0.25}.compute(pointset);
        writeRDFtoFile(fmt::format("rdf_{}.txt", maxdist), rdf);
    }

    delete args;
    return 0;
}