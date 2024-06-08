//
// Created by clemens on 08.06.24.
//

#ifndef UNICORNTK_RDF_HECK_HPP
#define UNICORNTK_RDF_HECK_HPP
namespace utk {
    class RDF {
    public:

        explicit RDF(
                bool toroidal = true,
                int32_t nbins = 1024,
                float_t maxdist = 0.25
        ) : nbins(nbins), toroidal(toroidal), maxdist(maxdist) {
            if (nbins == 0) {
                this->nbins = 1024;
            }
        }

        template<typename T>
        std::pair<std::vector<T>, std::vector<int32_t>> compute(const Pointset<T> &points) {
            const int npoints = points.Npts();

            std::vector<T> values(this->nbins, 0);
            std::vector<int32_t> frequencies(this->nbins, 0);

            auto distMetric = toroidal ? &RDF::toricDistance<T> : &RDF::distance<T>;

            for (int i = 0; i < npoints; ++i) {
                for (int j = i + 1; j < npoints; ++j) {
                    T dist = distMetric(points.Ndim(), points[i], points[j]);
                    //calculate which bin this distance goes into
                    int idx = std::floor((dist / maxdist) * static_cast<T>(nbins));
                    if (idx > 0 && idx < nbins)
                        frequencies[idx]++;
                }
            }

            //TODO: This probably deals with renormalization. Skip it for now because it is undocumented
            /*
            const float scale = npoints * (npoints - 1) / 2 * PI * rdf->dx * rdf->dx;
            for (int i = 0; i < rdf->size(); ++i)
                (*rdf)[i] = bins[i] / (scale * (2 * i + 1));
            */

            //save all the distances to the bins
            T binwidth = maxdist / static_cast<T>(nbins);

            for (int i = 0; i < nbins; ++i) {
                values[i] = static_cast<T>(i) * binwidth;
            }

            return std::make_pair(values, frequencies);
        }

    private:

        int32_t nbins;
        bool toroidal;
        float_t maxdist;

        template<typename T>
        static T toricDistance(uint32_t D, const T *a, const T *b) {
            T dist = 0;
            for (uint32_t d = 0; d < D; d++) {
                T tmp = std::abs(a[d] - b[d]);
                T diff = std::min(tmp, 1 - tmp);
                dist += diff * diff;
            }
            return std::sqrt(dist);
        }

        template<typename T>
        static T distance(uint32_t D, const T *a, const T *b) {
            T dist = 0;
            for (uint32_t d = 0; d < D; d++) {
                T diff = a[d] - b[d];
                dist += diff * diff;
            }
            return std::sqrt(dist);
        }
    };
}
#endif //UNICORNTK_RDF_HECK_HPP
