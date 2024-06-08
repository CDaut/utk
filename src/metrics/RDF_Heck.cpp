#include "CLI/App.hpp"
#include "utk/utils/MetricsArgumentParser.hpp"
#include "utk/metrics/RDF_Heck.hpp"

int main(int argc, char** argv)
{
    CLI::App app { "RDF calculator" };
    auto* margs = utk::add_arguments(app);
    bool toroidal = true;
    uint32_t bins = 1024;
    //TODO: This value is totally arbitrary
    float_t maxdist = 0.25;

    app.add_flag("--toroidal", toroidal, "When set, use toroidal distance. Euclidian otherwise")->default_val(toroidal);
    app.add_option("-b,--bins", bins, "Number of bins (0 means automatic)")->default_val(bins);
    app.add_option("-m,--maxdist", maxdist, "Maximum distance to consider")->default_val(maxdist);

    CLI11_PARSE(app, argc, argv);

    auto ptss = margs->GetAllPointsets();
    if (!utk::CheckPointsets(ptss))
        return 1;

    auto rslts = utk::RDF(bins, toroidal, maxdist).compute(ptss[0]);

    auto N = rslts.first.size();

    auto& ostream = margs->GetOutputStream();
    for (decltype(N) i = 0; i < N; i++)
        ostream << rslts.first[i] << " " << rslts.second[i] << '\n';

    delete margs;
}