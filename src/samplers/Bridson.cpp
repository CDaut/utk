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
#include <utk/samplers/SamplerBridson.hpp>

int main(int argc, char** argv)
{
    CLI::App app { "Bridson algorithm for Poisson disk sampling" };
    auto* args = utk::add_arguments(app);
    
    bool euclidean = true;
    uint32_t trials = 30;
    double distance = -1.0;

    app.add_flag("--euclidean", euclidean, "Use euclidean distance (default is wrap around)")->default_val(false);
    app.add_option("--trials", trials, "Number of trials candidates to try")->default_val(30);
    app.add_option("--distance", distance, "Minimal distance to use (will be scaled by N^-D !). If negative use precomputed sphere packing values (meaningfull only when D < 9)")->default_val(-1.0);

    CLI11_PARSE(app, argc, argv);
    
    std::vector<utk::Pointset<double>> pts = args->GetPointsets();
    
    utk::SamplerBridson brd(
        args->D, !euclidean, trials, distance
    );
    brd.setRandomSeed(args->seed);

    for (uint32_t i = 0; i < pts.size(); i++)
    {
        if(!brd.generateSamples(pts[i], args->N))
        {
            std::cerr << "Sampler returned non-zero output" << std::endl; // No log here, must be visible whatsoever
            return 1;
        }
    }

    args->WritePointsets(pts);

    delete args;
    return 0;
}