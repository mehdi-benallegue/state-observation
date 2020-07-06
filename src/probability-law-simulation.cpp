#include <state-observation/tools/probability-law-simulation.hpp>


namespace stateObservation
{
    namespace tools
    {

        boost::lagged_fibonacci1279 ProbabilityLawSimulation::gen_;

        double ProbabilityLawSimulation::getGaussianScalar(double std , double bias)
        {
            boost::normal_distribution<> g(0, 1);
            return g(gen_)*std+bias;

        }

        Matrix ProbabilityLawSimulation::getGaussianVector(const Matrix& std,
                                          const Matrix& bias,size_t rows, size_t cols)
        {
            boost::normal_distribution<> g(0, 1);
            Matrix ret= Matrix::Zero(rows,cols);
            for (size_t i=0;i<rows;++i)
            {
                for (size_t j=0; j<cols; ++j)
                    ret(i,j)=g(gen_);
            }
            ret=std*ret+bias;

            return ret;
        }
    }

}
