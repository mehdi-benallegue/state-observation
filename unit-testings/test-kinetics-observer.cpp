#include <iostream>

#include <state-observation/dynamics-estimators/kinetics-observer.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>

using namespace stateObservation;
int main()
{
    KineticsObserver o(7);

    
    std::cout << o.getStateSize()<<std::endl;
    
    Vector x0(o.getStateSize());
    x0.setZero();
    Vector xf(x0);
    Vector xs(x0);

    GaussianWhiteNoise random4(4);

    


    kine::Orientation ori;

    
    Vector4 q_i;

    q_i.setZero();
    
    q_i = random4.getNoisy(q_i);

    std::cout << "q_i " <<  std::endl;
    std::cout << q_i.transpose() << std::endl;

    q_i.normalize();

    std::cout << "q_i " <<  std::endl;
    std::cout << q_i.transpose() << std::endl;

    Quaternion q(q_i);
    AngleAxis aa(q);
    Matrix3 M(q.toRotationMatrix());

    std::cout << "q " <<  std::endl;
    std::cout << q.coeffs().transpose() << std::endl;
    std::cout << "M " << std::endl;
    std::cout << M << std::endl;
  


    ori = kine::Orientation(q);

    Matrix3 M2 = ori;

    std::cout << "q " <<  std::endl;
    std::cout << ori.toVector().transpose() << std::endl;
    std::cout << "M " << std::endl;
    std::cout << ori.matrix3() << std::endl;


    ori = q;

    
      


    std::cout <<"successfully built" <<std::endl;

    return 0;
}