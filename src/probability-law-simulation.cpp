#include <state-observation/tools/probability-law-simulation.hpp>


namespace stateObservation
{
    namespace tools
    {

        boost::lagged_fibonacci1279 ProbabilityLawSimulation::gen_;

        Matrix ProbabilityLawSimulation::getGaussianVector(const Matrix& std,
                                          const Matrix& bias,unsigned rows, unsigned cols)
        {
            boost::normal_distribution<> g(0, 1);
            Matrix ret= Matrix::Zero(rows,cols);
            for (unsigned i=0;i<rows;++i)
            {
                for (unsigned j=0; j<cols; ++j)
                    ret(i,j)=g(gen_);
            }
            ret=std*ret+bias;

            return ret;
        }
    }

}
