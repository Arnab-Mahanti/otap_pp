#pragma once
#include "OTypes.h"
#include <memory>
#include <valarray>
// Assumption that P,T are the independent variable
namespace OTAP
{
    struct Upstream
    {
        double P, T, M;
    };

    struct Downstream
    {
        double P, T, P0, T0, M;
    };

    using EdgeProp = Downstream;

    class FluidBase
    {
    private:
        /* data */
    public:
        FluidBase(/* args */) = default;
        virtual ~FluidBase() = default;
        virtual std::vector<double> Rho(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> Z(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> Pr(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> H(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> Gamma(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> Cp(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> k(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> Cv(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> mu(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> a(std::vector<double> P, std::vector<double> T) const { assert(false); }
        virtual std::vector<double> Lambda(std::vector<double> P, std::vector<double> T) const { assert(false); }

        // virtual std::valarray<double> Rho(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> Z(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> Pr(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> H(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> Gamma(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> Cp(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> k(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> Cv(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> mu(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> a(std::valarray<double> P, std::valarray<double> T) const { assert(false); }
        // virtual std::valarray<double> Lambda(std::valarray<double> P, std::valarray<double> T) const { assert(false); }

        virtual double Rho(double P, double T) const { assert(false); }
        virtual double Z(double P, double T) const { assert(false); }
        virtual double Pr(double P, double T) const { assert(false); }
        virtual double H(double P, double T) const { assert(false); }
        virtual double Gamma(double P, double T) const { assert(false); }
        virtual double Cp(double P, double T) const { assert(false); }
        virtual double k(double P, double T) const { assert(false); }
        virtual double Cv(double P, double T) const { assert(false); }
        virtual double mu(double P, double T) const { assert(false); }
        virtual double a(double P, double T) const { assert(false); }
        virtual double Lambda(double P, double T) const { assert(false); }

        // FIXME: add vector versions
        virtual Downstream GetShockDownStream(const Upstream &upsteam) const { assert(false); };
        virtual EdgeProp GetEdgeProperties(const Upstream &upsteam, double cp) const { assert(false); };
    };

    // TODO: Separate out transport properties
    class PerfectGas : public FluidBase
    {
        const double m_M = 29.6;
        const double m_R = 8.314 / m_M * 1000;

    public:
        virtual std::vector<double> Rho(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Z(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Pr(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> H(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Gamma(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Cp(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> k(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Cv(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> mu(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> a(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Lambda(std::vector<double> P, std::vector<double> T) const override;

        // virtual std::valarray<double> Rho(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Z(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Pr(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> H(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Gamma(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Cp(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> k(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Cv(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> mu(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> a(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Lambda(std::valarray<double> P, std::valarray<double> T) const override;

        virtual double Rho(double P, double T) const override;
        virtual double Z(double P, double T) const override;
        virtual double Pr(double P, double T) const override;
        virtual double H(double P, double T) const override;
        virtual double Gamma(double P, double T) const override;
        virtual double Cp(double P, double T) const override;
        virtual double k(double P, double T) const override;
        virtual double Cv(double P, double T) const override;
        virtual double mu(double P, double T) const override;
        virtual double a(double P, double T) const override;
        virtual double Lambda(double P, double T) const override;

        virtual Downstream GetShockDownStream(const Upstream &upsteam) const override;
        virtual EdgeProp GetEdgeProperties(const Upstream &upsteam, double cp) const override;

        PerfectGas(double M = 29.6) : m_M(M){};
    };

    class Hansen : public FluidBase
    {
    private:
        /* data */
        // static inline const std::vector<double> m_logPData = {2, 1, 0, -1, -2, -3, -4};
        static inline const std::vector<double> m_logPData = {-4, -3, -2, -1, 0, 1, 2};

        static inline const std::vector<double> m_TData = {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000};

        static inline const std::vector<double> m_ZData = {1, 1, 1, 1, 1, 1, 1,
                                                           1, 1, 1, 1, 1, 1, 1,
                                                           1, 1, 1, 1, 1, 1, 1,
                                                           1.016, 1.005, 1.002, 1.001, 1, 1, 1,
                                                           1.163, 1.088, 1.033, 1.011, 1.004, 1.001, 1,
                                                           1.2, 1.192, 1.149, 1.072, 1.026, 1.009, 1.003,
                                                           1.211, 1.203, 1.197, 1.167, 1.092, 1.035, 1.012,
                                                           1.287, 1.228, 1.208, 1.198, 1.165, 1.089, 1.033,
                                                           1.577, 1.337, 1.245, 1.213, 1.196, 1.149, 1.071,
                                                           1.91, 1.622, 1.359, 1.252, 1.214, 1.186, 1.118,
                                                           1.99, 1.898, 1.599, 1.348, 1.248, 1.208, 1.159,
                                                           2.008, 1.983, 1.849, 1.529, 1.316, 1.235, 1.189,
                                                           2.032, 2.006, 1.961, 1.752, 1.437, 1.279, 1.214,
                                                           2.088, 2.027, 1.997, 1.904, 1.607, 1.351, 1.243,
                                                           2.21, 2.067, 2.017, 1.971, 1.778, 1.457, 1.284,
                                                           2.446, 2.144, 2.044, 2.001, 1.896, 1.59, 1.341,
                                                           2.826, 2.284, 2.09, 2.023, 1.959, 1.727, 1.418,
                                                           3.282, 2.51, 2.166, 2.05, 1.993, 1.838, 1.512,
                                                           3.645, 2.832, 2.286, 2.09, 2.018, 1.914, 1.616,
                                                           3.843, 3.202, 2.462, 2.149, 2.042, 1.962, 1.718,
                                                           3.932, 3.526, 2.7, 2.234, 2.071, 1.993, 1.807,
                                                           3.969, 3.745, 2.983, 2.351, 2.111, 2.018, 1.876,
                                                           3.985, 3.867, 3.272, 2.505, 2.163, 2.042, 1.927,
                                                           3.993, 3.931, 3.52, 2.694, 2.232, 2.067, 1.965,
                                                           3.996, 3.963, 3.7, 2.91, 2.318, 2.098, 1.993,
                                                           3.998, 3.979, 3.818, 3.135, 2.426, 2.135, 2.017,
                                                           3.999, 3.988, 3.889, 3.347, 2.553, 2.18, 2.039,
                                                           3.999, 3.993, 3.932, 3.527, 2.7, 2.233, 2.062,
                                                           4, 3.996, 3.957, 3.667, 2.861, 2.297, 2.086,
                                                           4, 3.997, 3.973, 3.769, 3.028, 2.372, 2.113};

        static inline const std::vector<double> m_PrData = {0.738, 0.738, 0.738, 0.738, 0.738, 0.738, 0.738,
                                                            0.756, 0.756, 0.756, 0.756, 0.756, 0.756, 0.756,
                                                            0.767, 0.767, 0.767, 0.767, 0.767, 0.767, 0.767,
                                                            0.614, 0.668, 0.724, 0.766, 0.773, 0.773, 0.773,
                                                            0.771, 0.654, 0.611, 0.645, 0.696, 0.751, 0.762,
                                                            0.714, 0.745, 0.74, 0.636, 0.627, 0.68, 0.74,
                                                            0.606, 0.658, 0.737, 0.744, 0.66, 0.631, 0.678,
                                                            0.587, 0.58, 0.619, 0.759, 0.762, 0.662, 0.64,
                                                            0.764, 0.611, 0.578, 0.61, 0.752, 0.743, 0.654,
                                                            0.993, 0.799, 0.624, 0.581, 0.611, 0.767, 0.702,
                                                            0.871, 0.989, 0.785, 0.617, 0.583, 0.62, 0.748,
                                                            0.455, 0.891, 0.969, 0.736, 0.602, 0.592, 0.763,
                                                            0.392, 0.464, 0.955, 0.906, 0.673, 0.592, 0.61,
                                                            0.361, 0.404, 0.83, 0.986, 0.796, 0.62, 0.593,
                                                            0.342, 0.371, 0.424, 0.969, 0.927, 0.688, 0.595,
                                                            0.322, 0.351, 0.387, 0.648, 0.983, 0.788, 0.62,
                                                            0.279, 0.355, 0.363, 0.411, 0.943, 0.891, 0.666,
                                                            0.2, 0.316, 0.348, 0.382, 0.807, 0.961, 730,
                                                            0.114, 0.279, 0.336, 0.364, 0.497, 0.966, 0.806,
                                                            0.0576, 0.216, 0.319, 0.348, 0.429, 0.872, 0.886,
                                                            0.0314, 0.145, 0.295, 0.339, 0.404, 0.532, 0.937,
                                                            0.0213, 0.0877, 0.254, 0.327, 0.382, 0.463, 0.955,
                                                            0.0167, 0.0524, 0.201, 0.312, 0.369, 0.434, 947,
                                                            0.0143, 0.0346, 0.146, 0.292, 0.355, 0.412, 0.908,
                                                            0.0129, 0.0238, 0.101, 0.263, 0.343, 0.396, 0.728,
                                                            0.0121, 0.019, 0.0688, 0.227, 0.333, 0.383, 0.525,
                                                            0.011, 0.0162, 0.047, 0.185, 0.319, 0.369, 0.438,
                                                            0.0108, 0.0149, 0.0345, 0.144, 0.302, 0.36, 0.421,
                                                            0.0109, 0.013, 0.0245, 0.0986, 0.277, 0.349, 0.401,
                                                            0.011, 0.012, 0.0129, 0.0819, 0.253, 0.341, 0.394};

        static inline const std::vector<double> m_HData = {505000, 505000, 505000, 505000, 505000, 505000, 505000,
                                                           1050000, 1050000, 1050000, 1050000, 1050000, 1050000, 1050000,
                                                           1640000, 1640000, 1640000, 1640000, 1640000, 1640000, 1640000,
                                                           2530000, 2340000, 2280000, 2260000, 2250000, 2250000, 2250000,
                                                           5750000, 4420000, 3450000, 3060000, 2940000, 2890000, 2880000,
                                                           7050000, 6910000, 6140000, 4780000, 3970000, 3660000, 3560000,
                                                           8070000, 7810000, 7650000, 7110000, 5780000, 4770000, 4360000,
                                                           11300000, 9290000, 8640000, 8360000, 7740000, 6380000, 5400000,
                                                           21700000, 13600000, 10500000, 9470000, 9020000, 8120000, 6720000,
                                                           33700000, 23900000, 15000000, 11400000, 10200000, 9500000, 8220000,
                                                           37200000, 34100000, 23900000, 15400000, 12000000, 10700000, 9680000,
                                                           38800000, 37800000, 33200000, 22300000, 15000000, 12200000, 11000000,
                                                           40900000, 39600000, 38000000, 30700000, 19900000, 14400000, 12400000,
                                                           44600000, 41600000, 40200000, 36900000, 26500000, 17600000, 14000000,
                                                           52200000, 44600000, 42100000, 40200000, 33300000, 22000000, 16000000,
                                                           65800000, 49700000, 44400000, 42300000, 38400000, 27500000, 18700000,
                                                           87200000, 58200000, 47800000, 44300000, 41700000, 33300000, 22200000,
                                                           113000000, 71500000, 52900000, 46700000, 44000000, 38200000, 26400000,
                                                           134000000, 90000000, 60500000, 49900000, 46100000, 42000000, 31000000,
                                                           146000000, 111000000, 71100000, 54100000, 48300000, 44800000, 35700000,
                                                           153000000, 130000000, 85200000, 59800000, 50900000, 47100000, 39900000,
                                                           156000000, 144000000, 102000000, 67300000, 54100000, 49300000, 43500000,
                                                           159000000, 152000000, 119000000, 76800000, 58000000, 51400000, 46500000,
                                                           161000000, 157000000, 134000000, 88500000, 62800000, 53800000, 49000000,
                                                           162000000, 161000000, 146000000, 102000000, 68700000, 56500000, 51200000,
                                                           164000000, 163000000, 154000000, 116000000, 75800000, 59600000, 53400000,
                                                           166000000, 165000000, 160000000, 129000000, 84200000, 63100000, 55400000,
                                                           167000000, 167000000, 164000000, 141000000, 93600000, 67100000, 57600000,
                                                           169000000, 169000000, 167000000, 150000000, 104000000, 71800000, 59900000,
                                                           171000000, 170000000, 169000000, 157000000, 115000000, 77200000, 62400000};

        static inline const std::vector<double> m_GammaData = {1.387, 1.387, 1.387, 1.387, 1.387, 1.387, 1.387,
                                                               1.337, 1.337, 1.337, 1.337, 1.337, 1.337, 1.337,
                                                               1.306, 1.31, 1.312, 1.312, 1.312, 1.312, 1.312,
                                                               1.153, 1.209, 1.26, 1.286, 1.296, 1.299, 1.3,
                                                               1.157, 1.152, 1.161, 1.202, 1.249, 1.277, 1.288,
                                                               1.304, 1.239, 1.181, 1.178, 1.195, 1.235, 1.266,
                                                               1.176, 1.252, 1.27, 1.212, 1.202, 1.211, 1.241,
                                                               1.133, 1.15, 1.213, 1.26, 1.23, 1.223, 1.23,
                                                               1.19, 1.155, 1.154, 1.204, 1.251, 1.243, 1.24,
                                                               1.168, 1.203, 1.172, 1.166, 1.212, 1.252, 1.256,
                                                               1.257, 1.183, 1.214, 1.182, 1.183, 1.231, 1.262,
                                                               1.266, 1.237, 1.202, 1.221, 1.19, 1.206, 1.253,
                                                               1.188, 1.265, 1.217, 1.228, 1.22, 1.201, 1.235,
                                                               1.155, 1.21, 1.258, 1.216, 1.246, 1.217, 1.223,
                                                               1.164, 1.173, 1.237, 1.237, 1.244, 1.243, 1.223,
                                                               1.201, 1.168, 1.201, 1.252, 1.235, 1.264, 1.235,
                                                               1.242, 1.188, 1.183, 1.235, 1.243, 1.267, 1.255,
                                                               1.244, 1.224, 1.185, 1.213, 1.252, 1.26, 1.275,
                                                               1.216, 1.256, 1.203, 1.201, 1.248, 1.255, 1.288,
                                                               1.211, 1.263, 1.232, 1.201, 1.236, 1.259, 1.291,
                                                               1.256, 1.244, 1.263, 1.213, 1.226, 1.262, 1.287,
                                                               1.339, 1.23, 1.281, 1.233, 1.222, 1.261, 1.282,
                                                               1.427, 1.243, 1.28, 1.258, 1.226, 1.256, 1.28,
                                                               1.491, 1.288, 1.267, 1.283, 1.236, 1.252, 1.282,
                                                               1.528, 1.352, 1.257, 1.301, 1.251, 1.25, 1.284,
                                                               1.548, 1.419, 1.263, 1.307, 1.271, 1.251, 1.285,
                                                               1.558, 1.472, 1.288, 1.303, 1.291, 1.257, 1.284,
                                                               1.563, 1.509, 1.329, 1.295, 1.311, 1.266, 1.284,
                                                               1.565, 1.532, 1.377, 1.29, 1.326, 1.278, 1.284,
                                                               1.567, 1.547, 1.425, 1.293, 1.336, 1.293, 1.286};

        static inline const std::vector<double> m_CpData = {1030.3, 1030.3, 1030.3, 1030.3, 1030.3, 1030.3, 1030.3,
                                                            1136.5, 1136.5, 1136.5, 1136.5, 1136.5, 1136.5, 1136.5,
                                                            1231.2, 1214, 1208.3, 1205.4, 1205.4, 1205.4, 1205.4,
                                                            3320.6, 1931.5, 1463.7, 1311.6, 1265.7, 1248.4, 1245.6,
                                                            5748.6, 6888, 3897.5, 2189.8, 1567, 1360.4, 1294.4,
                                                            1552.7, 2261.6, 5068.4, 4775.7, 2760.9, 1788, 1443.6,
                                                            3068, 1905.7, 1868.4, 3533, 4276.3, 2726.5, 1799.5,
                                                            11703.9, 4724, 2436.6, 1983.2, 3171.4, 3578.9, 2364.9,
                                                            29127.6, 14080.2, 5648.2, 2758.1, 2189.8, 3174.2, 2904.4,
                                                            13629.6, 25020.7, 13265.1, 5447.3, 2749.5, 2459.6, 3022.1,
                                                            3679.3, 12794.5, 20864.9, 10687.9, 4548.9, 2600.2, 2743.7,
                                                            3280.4, 4350.9, 14169.2, 16531.2, 7734.7, 3553.1, 2597.4,
                                                            5490.3, 3375.1, 5912.2, 15604.2, 11804.3, 5255, 2884.4,
                                                            10567.3, 4701.1, 3699.4, 8948.7, 14235.2, 7605.5, 3587.5,
                                                            20095.7, 7754.7, 4015.1, 4970.8, 12252, 10099.5, 4729.8,
                                                            34899.2, 13092.9, 5567.8, 3926.2, 8113.5, 11574.7, 6147.5,
                                                            49754.3, 21329.8, 8279.9, 4279.2, 5332.5, 10946.2, 7662.9,
                                                            49742.8, 31989, 12372.6, 5424.3, 4256.2, 8719.1, 8911.4,
                                                            32993.5, 41310.8, 17992, 7249.6, 4204.6, 6457.5, 9442.3,
                                                            17303.2, 42217.7, 24805.4, 9789.6, 4744.1, 5011, 9020.4,
                                                            9181.1, 32936.1, 31380.6, 13067.1, 5699.8, 4359.5, 7881,
                                                            5717, 21094.5, 34950.9, 17004.8, 7020, 4256.2, 9419.3,
                                                            4284.9, 12679.7, 33142.8, 21266.7, 8687.5, 4508.8, 5441.5,
                                                            3679.3, 8010.2, 26800.1, 25132.6, 10685, 5008.2, 4703.9,
                                                            3412.4, 5645.3, 19332.3, 27554.9, 12963.8, 5705.6, 4307.9,
                                                            3291.9, 4465.7, 13302.5, 27595, 15411.9, 6569.4, 4181.6,
                                                            3231.6, 3871.6, 9270.1, 25138.3, 17828.4, 7582.5, 4247.6,
                                                            3202.9, 3561.7, 6807.6, 21071.5, 19920.7, 8727.7, 4454.2,
                                                            3185.7, 3395.2, 5355.4, 16733.8, 21341.3, 9979, 4769.9,
                                                            3180, 3303.4, 4505.9, 12814.6, 21797.7, 11299.2, 5171.7};

        static inline const std::vector<double> m_kData = {0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372,
                                                           0.0624, 0.0624, 0.0624, 0.0624, 0.0624, 0.0624, 0.0624,
                                                           0.0826, 0.0826, 0.0826, 0.0826, 0.0826, 0.0826, 0.0826,
                                                           0.337, 0.1765, 0.1233, 0.1057, 0.0994, 0.0994, 0.0994,
                                                           0.5249, 0.732, 0.4419, 0.2386, 0.1545, 0.1257, 0.1199,
                                                           0.1543, 0.2306, 0.6339, 0.577, 0.3369, 0.203, 0.1496,
                                                           0.4241, 0.2412, 0.1965, 0.4527, 0.5396, 0.3601, 0.2219,
                                                           1.844, 0.741, 0.357, 0.1963, 0.3668, 0.4834, 0.33,
                                                           3.981, 2.304, 0.958, 0.4333, 0.2237, 0.3994, 0.4203,
                                                           1.633, 3.558, 2.294, 0.9903, 0.4538, 0.2754, 0.4234,
                                                           0.5131, 1.652, 3.218, 1.987, 0.8682, 0.4218, 0.3566,
                                                           0.9311, 0.6005, 1.985, 2.842, 1.545, 0.6869, 0.2927,
                                                           2.053, 0.9357, 0.8677, 2.432, 2.292, 1.103, 0.5295,
                                                           4.417, 1.798, 0.5386, 1.366, 2.577, 1.637, 0.7702,
                                                           8.831, 3.398, 1.467, 1.007, 2.083, 2.124, 1.074,
                                                           14.92, 6.104, 2.461, 0.6016, 1.372, 2.32, 1.444,
                                                           19.53, 10.22, 4.01, 1.764, 0.9255, 2.096, 1.789,
                                                           17.53, 14.72, 6.372, 2.622, 0.6092, 1.642, 2.037,
                                                           10.55, 17.79, 9.432, 3.9, 1.279, 1.226, 2.09,
                                                           5.049, 16.72, 12.78, 5.561, 2.159, 0.8026, 1.946,
                                                           2.315, 12.11, 15.49, 7.602, 2.952, 1.314, 1.684,
                                                           1.203, 7.288, 16.27, 10, 3.969, 1.754, 1.416,
                                                           0.7201, 4.032, 14.53, 12.33, 5.189, 2.22, 1.184,
                                                           0.4956, 2.327, 11.17, 14.15, 6.688, 2.805, 1.037,
                                                           0.3993, 1.346, 7.69, 14.87, 8.176, 3.436, 0.7381,
                                                           0.3656, 0.9197, 5.072, 14.29, 9.761, 4.231, 1.303,
                                                           0.3357, 0.7307, 3.281, 12.52, 11.19, 5.124, 2.286,
                                                           0.3713, 0.6461, 2.242, 10.18, 12.24, 6.063, 2.678,
                                                           0.3979, 0.579, 1.489, 7.244, 12.92, 7.482, 3.272,
                                                           0.4252, 0.5648, 0.737, 5.891, 12.7, 8.194, 3.624};

        static inline const std::vector<double> m_muData = {2.67e-05, 2.67e-05, 2.67e-05, 2.67e-05, 2.67e-05, 2.67e-05, 2.67e-05,
                                                            4.16e-05, 4.16e-05, 4.16e-05, 4.16e-05, 4.16e-05, 4.16e-05, 4.16e-05,
                                                            5.27e-05, 5.27e-05, 5.27e-05, 5.27e-05, 5.27e-05, 5.27e-05, 5.27e-05,
                                                            6.19e-05, 6.19e-05, 6.19e-05, 6.19e-05, 6.19e-05, 6.19e-05, 6.19e-05,
                                                            7e-05, 7e-05, 7e-05, 7e-05, 7e-05, 7e-05, 7e-05,
                                                            7.72e-05, 7.72e-05, 7.72e-05, 7.72e-05, 7.72e-05, 7.72e-05, 7.72e-05,
                                                            8.47e-05, 8.47e-05, 8.47e-05, 8.43e-05, 8.41e-05, 8.39e-05, 8.38e-05,
                                                            9.28e-05, 9.21e-05, 9.19e-05, 9.18e-05, 9.14e-05, 9.07e-05, 9.02e-05,
                                                            0.000105, 0.000101, 9.93e-05, 9.89e-05, 9.85e-05, 9.78e-05, 9.67e-05,
                                                            0.000119, 0.000102, 0.000109, 0.000106, 0.000106, 0.000105, 0.000103,
                                                            0.00013, 0.000129, 0.000122, 0.000115, 0.000113, 0.000112, 0.00011,
                                                            0.00014, 0.00014, 0.000137, 0.000128, 0.000121, 0.000119, 0.000117,
                                                            0.000147, 0.000149, 0.000148, 0.000144, 0.000132, 0.000126, 0.000124,
                                                            0.000152, 0.000157, 0.000159, 0.000156, 0.000145, 0.000135, 0.000131,
                                                            0.000151, 0.000163, 0.000167, 0.000166, 0.00016, 0.000147, 0.000139,
                                                            0.000138, 0.000165, 0.000174, 0.000177, 0.000173, 0.00016, 0.000147,
                                                            0.00011, 0.000161, 0.000179, 0.000184, 0.000184, 0.000174, 0.000158,
                                                            7.08e-05, 0.000146, 0.00018, 0.000191, 0.000195, 0.000187, 0.00017,
                                                            3.68e-05, 0.00012, 0.000176, 0.000196, 0.000203, 0.0002, 0.000183,
                                                            1.71e-05, 8.6e-05, 0.000165, 0.000199, 0.000209, 0.000212, 0.000197,
                                                            8.15e-06, 5.35e-05, 0.000146, 0.000198, 0.000215, 0.000222, 0.00021,
                                                            4.4e-06, 3.04e-05, 0.000119, 0.000192, 0.000219, 0.000228, 0.000223,
                                                            2.8e-06, 1.68e-05, 8.87e-05, 0.000181, 0.000221, 0.000235, 0.000234,
                                                            1.9e-06, 1e-05, 6.14e-05, 0.000165, 0.000221, 0.000241, 0.000246,
                                                            1.46e-06, 5.83e-06, 4.03e-05, 0.000143, 0.000217, 0.000246, 0.000256,
                                                            1.32e-06, 3.97e-06, 2.61e-05, 0.000118, 0.000211, 0.000249, 0.000261,
                                                            1.18e-06, 3.03e-06, 1.69e-05, 9.22e-05, 0.0002, 0.000251, 0.000269,
                                                            1.2e-06, 2.57e-06, 1.15e-05, 7e-05, 0.000186, 0.000252, 0.000274,
                                                            1.4e-06, 2.27e-06, 7.34e-06, 4.68e-05, 0.000164, 0.000247, 0.00028,
                                                            1.42e-06, 2.13e-06, 2.84e-06, 3.77e-05, 0.000147, 0.000247, 0.000285};

        static inline const std::vector<double> m_aData = {1.387, 1.387, 1.387, 1.387, 1.387, 1.387, 1.387,
                                                           1.337, 1.337, 1.337, 1.337, 1.337, 1.337, 1.337,
                                                           1.306, 1.31, 1.312, 1.312, 1.312, 1.312, 1.312,
                                                           1.144, 1.206, 1.259, 1.285, 1.296, 1.299, 1.3,
                                                           1.132, 1.119, 1.144, 1.196, 1.247, 1.276, 1.288,
                                                           1.302, 1.232, 1.15, 1.147, 1.181, 1.229, 1.265,
                                                           1.171, 1.25, 1.265, 1.187, 1.166, 1.192, 1.234,
                                                           1.097, 1.137, 1.208, 1.254, 1.204, 1.187, 1.212,
                                                           1.092, 1.101, 1.133, 1.196, 1.241, 1.21, 1.208,
                                                           1.124, 1.103, 1.111, 1.143, 1.202, 1.233, 1.217,
                                                           1.249, 1.133, 1.113, 1.124, 1.161, 1.217, 1.229,
                                                           1.263, 1.225, 1.135, 1.124, 1.141, 1.186, 1.23,
                                                           1.183, 1.26, 1.193, 1.136, 1.137, 1.165, 1.214,
                                                           1.14, 1.205, 1.249, 1.167, 1.142, 1.155, 1.195,
                                                           1.128, 1.162, 1.231, 1.216, 1.156, 1.154, 1.182,
                                                           1.13, 1.143, 1.193, 1.242, 1.181, 1.159, 1.175,
                                                           1.136, 1.14, 1.168, 1.228, 1.214, 1.169, 1.174,
                                                           1.145, 1.144, 1.157, 1.203, 1.237, 1.185, 1.176,
                                                           1.157, 1.151, 1.155, 1.185, 1.237, 1.206, 1.181,
                                                           1.185, 1.159, 1.159, 1.176, 1.225, 1.228, 1.19,
                                                           1.244, 1.169, 1.165, 1.173, 1.212, 1.242, 1.201,
                                                           1.334, 1.187, 1.173, 1.175, 1.202, 1.246, 1.216,
                                                           1.424, 1.221, 1.181, 1.18, 1.198, 1.242, 1.232,
                                                           1.489, 1.276, 1.192, 1.187, 1.197, 1.236, 1.247,
                                                           1.528, 1.346, 1.207, 1.194, 1.198, 1.231, 1.257,
                                                           1.548, 1.415, 1.232, 1.202, 1.202, 1.227, 1.263,
                                                           1.558, 1.47, 1.269, 1.21, 1.208, 1.226, 1.265,
                                                           1.563, 1.508, 1.317, 1.221, 1.214, 1.226, 1.265,
                                                           1.565, 1.532, 1.37, 1.235, 1.221, 1.228, 1.264,
                                                           1.567, 1.546, 1.42, 1.254, 1.228, 1.232, 1.263};

    public:
        virtual std::vector<double> Rho(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Z(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Pr(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> H(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Gamma(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Cp(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> k(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Cv(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> mu(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> a(std::vector<double> P, std::vector<double> T) const override;
        virtual std::vector<double> Lambda(std::vector<double> P, std::vector<double> T) const override;

        // virtual std::valarray<double> Rho(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Z(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Pr(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> H(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Gamma(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Cp(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> k(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Cv(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> mu(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> a(std::valarray<double> P, std::valarray<double> T) const override;
        // virtual std::valarray<double> Lambda(std::valarray<double> P, std::valarray<double> T) const override;

        virtual double Rho(double P, double T) const override;
        virtual double Z(double P, double T) const override;
        virtual double Pr(double P, double T) const override;
        virtual double H(double P, double T) const override;
        virtual double Gamma(double P, double T) const override;
        virtual double Cp(double P, double T) const override;
        virtual double k(double P, double T) const override;
        virtual double Cv(double P, double T) const override;
        virtual double mu(double P, double T) const override;
        virtual double a(double P, double T) const override;
        virtual double Lambda(double P, double T) const override;

        virtual Downstream GetShockDownStream(const Upstream &upsteam) const override;
        virtual EdgeProp GetEdgeProperties(const Upstream &upsteam, double cp) const override;
    };

    template <typename... Args>
    std::shared_ptr<FluidBase> make_fluid(FluidType T, Args &&...args)
    {
        switch (T)
        {
        case FluidType::Hansen_Air:
            return safe_make_shared<Hansen>(std::forward<Args>(args)...);
        case FluidType::Perfect:
            return safe_make_shared<PerfectGas>(std::forward<Args>(args)...);
        default:
            assert(false);
            return nullptr;
        }
    }

    // class ShockBase
    // {
    // public:
    //     virtual Downstream GetDownstream(const Upstream &upstream) const { assert(false); };
    // };

    // class Ideal_Shock : public ShockBase
    // {

    // public:
    //     virtual Downstream GetDownstream(const Upstream &upstream) const override;
    // };

    // class Real_Shock : public ShockBase
    // {
    // private:
    //     Hansen m_hansen;

    // public:
    //     virtual Downstream GetDownstream(const Upstream &upstream) const override;
    // };

    // alias

    using Fluid = std::shared_ptr<FluidBase>;
} // namespace OTAP
