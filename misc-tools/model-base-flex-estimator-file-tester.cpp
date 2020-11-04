#include <fstream>
#include <iostream>

#include <stdlib.h>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/examples/offline-model-base-flex-estimation.hpp>

#include <time.h>

using namespace stateObservation;

int main(int argc, char * argv[])
{
  if(argc != 6)
  {
    std::cout << "Parameters : measurements inputs numberOfContacts dt m" << std::endl;
    std::cout << "Leaving" << std::endl;
  }
  else
  {
    IndexedVectorArray y, u, numberOfContacts;

    std::cout << "Read Sensors filename: " << argv[1] << std::endl;
    y.readVectorsFromFile(argv[1]);
    std::cout << "Sensors loaded, size:" << y.size() << std::endl;
    std::cout << "Read Inputs, filename: " << argv[2] << std::endl;
    u.readVectorsFromFile(argv[2]);
    std::cout << "Inputs loaded, size:" << u.size() << std::endl;

    std::cout << "Read contact number, filename: " << argv[3] << std::endl;
    numberOfContacts.readVectorsFromFile(argv[3]);
    std::cout << "Contact numbers loaded, size:" << numberOfContacts.size() << std::endl;

    // std::cout <<"numberOfContacts " <<numberOfContacts.getLastIndex()<< " "
    //           << numberOfContacts[numberOfContacts.getLastIndex()].size() << " "
    //           << numberOfContacts[numberOfContacts.getLastIndex()] << " "
    //           << std::endl;//<< " "<< numberOfContacts[2]<< " "<< numberOfContacts[3]<< " ";

    Matrix xh0 = Vector::Zero(flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::size, 1);

    double dt = atof(argv[4]);
    double mass = atof(argv[5]);

    std::cout << "Rebuiding state" << std::endl;
    IndexedVectorArray xhat = examples::offlineModelBaseFlexEstimation(
        y, u, xh0, numberOfContacts, dt, mass, true, IndexedMatrixArray(), IndexedMatrixArray(), Matrix3::Zero(),
        Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero(), 0x0, 0x0, 0x0, 0x0, 1);
    std::cout << "State rebuilt" << std::endl;

    xhat.writeInFile("xhat");
  }
  return 0;
}
