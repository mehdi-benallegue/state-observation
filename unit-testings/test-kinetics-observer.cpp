#include <iostream>

#include <state-observation/dynamics-estimators/kinetics-observer.hpp>

using namespace stateObservation;
int main()
{
    KineticsObserver o;

    Vector x0(o.getStateSize());
    x0.setZero();
    Vector xf(x0);
    Vector xs(x0);

      


    std::cout <<"successfully built" <<std::endl;

    return 0;
}