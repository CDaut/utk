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
        std::pair<std::vector<T>, std::vector<T>> compute(const Pointset<T> &points) {
            constexpr T PI = 3.141592653589793238;
            const int npoints = points.Npts();

            std::vector<T> values(this->nbins, 0);
            std::vector<int32_t> bins(this->nbins, 0);
            std::vector<T> frequencies(this->nbins, 0);

            auto distMetric = toroidal ? &RDF::toricDistance<T> : &RDF::distance<T>;

            for (int i = 0; i < npoints; ++i) {
                for (int j = i + 1; j < npoints; ++j) {
                    T dist = distMetric(points.Ndim(), points[i], points[j]);
                    //calculate which bin this distance goes into
                    int idx = std::floor((dist / maxdist) * static_cast<T>(nbins));
                    if (idx > 0 && idx < nbins)
                        bins[idx]++;
                }
            }

            //renormalization
            const T scale = npoints * (npoints - 1) / 2 * PI * maxdist * maxdist;

            for (int i = 0; i < nbins; ++i) {
                frequencies[i] = static_cast<T>(bins[i]) / (scale * (2 * i + 1));
            }


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
