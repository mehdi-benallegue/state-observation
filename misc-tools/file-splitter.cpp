#include <fstream>
#include <iostream>
#include <stdexcept> // std::invalid_argument

#include <stdlib.h>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/examples/offline-model-base-flex-estimation.hpp>

#include <time.h>

using namespace stateObservation::flexibilityEstimation;
typedef IMUElasticLocalFrameDynamicalSystem::state state;
typedef IMUElasticLocalFrameDynamicalSystem::input input;

using namespace stateObservation;

/// ///////////////////////////////////////////////////////////////

int main(int argc, char * argv[])

{
  try
  {
    if(argc % 2 == 0) throw std::invalid_argument("argc");

    for(int i = 1; i < argc; i = i + 2) /// i is incremented twice
    {
      std::string type, filename;
      type = argv[i];
      filename = argv[i + 1];
      if(type == "-u")
      {

        ///-------------- input----------------------------------

        IndexedVectorArray u;
        std::cout << "Input " << std::endl;
        std::cout << " filename: " << filename << std::endl;
        u.readVectorsFromFile(filename.c_str());
        std::cout << "Inputs loaded, size:" << u.size() << std::endl;

        Vector ui;

        std::cout << "Splitting" << std::endl;

        IndexedVectorArray posCom, velCom, accCom, inertia, angMoment, dotInertia, dotAngMoment, posIMU, oriIMU,
            linVelIMU, angVelIMU, linAccIMU, additionalForces, contact1, contact2;
        for(TimeIndex i = u.getFirstIndex(); i < u.getNextIndex(); ++i)
        {
          ui = u[i];
          posCom.setValue(ui.segment<3>(input::posCom), i);
          velCom.setValue(ui.segment<3>(input::velCom), i);
          accCom.setValue(ui.segment<3>(input::accCom), i);
          inertia.setValue(ui.segment<6>(input::inertia), i);
          angMoment.setValue(ui.segment<3>(input::angMoment), i);
          dotInertia.setValue(ui.segment<6>(input::dotInertia), i);
          dotAngMoment.setValue(ui.segment<3>(input::dotAngMoment), i);
          posIMU.setValue(ui.segment<3>(input::posIMU), i);
          oriIMU.setValue(ui.segment<3>(input::oriIMU), i);
          linVelIMU.setValue(ui.segment<3>(input::linVelIMU), i);
          angVelIMU.setValue(ui.segment<3>(input::angVelIMU), i);
          linAccIMU.setValue(ui.segment<3>(input::linAccIMU), i);
          additionalForces.setValue(ui.segment<6>(input::additionalForces), i);

          Vector c1, c2;
          unsigned numberOfContacts = 0;

          if(u.size() > input::sizeBase)
          {
            c1.resize(13);
            numberOfContacts++;
            c1 << numberOfContacts, ui.segment<12>(input::contacts);
          }
          else
          {
            c1.resize(1);
            c1 << numberOfContacts;
            contact1.setValue(c1, i);
          }

          if(u[i].size() > input::sizeBase + 12)
          {
            c2.resize(13);
            numberOfContacts++;
            c1(0) = numberOfContacts;
            c2 << numberOfContacts, ui.segment<12>(input::contacts + 12);
          }
          else
          {
            c2.resize(1);
            c2 << numberOfContacts;
          }

          contact1.setValue(c1, i);
          contact2.setValue(c2, i);
        }

        std::cout << "Splitting finished" << std::endl;

        std::cout << "Writing files" << std::endl;

        posCom.writeInFile(filename + "-00-posCom");
        velCom.writeInFile(filename + "-01-velCom");
        accCom.writeInFile(filename + "-02-accCom");
        inertia.writeInFile(filename + "-03-inertia");
        angMoment.writeInFile(filename + "-04-angMoment");
        dotInertia.writeInFile(filename + "-05-dotInertia");
        dotAngMoment.writeInFile(filename + "-06-dotAngMoment");
        posIMU.writeInFile(filename + "-07-posIMU");
        oriIMU.writeInFile(filename + "-08-oriIMU");
        linVelIMU.writeInFile(filename + "-09-linVelIMU");
        angVelIMU.writeInFile(filename + "-10-angVelIMU");
        linAccIMU.writeInFile(filename + "-11-linAccIMU");
        additionalForces.writeInFile(filename + "-12-additionalForces");
        contact1.writeInFile(filename + "-13-contact1");
        contact2.writeInFile(filename + "-14-contact2");

        std::cout << "Done " << std::endl;
      }
      else if(type == "-y")
      {

        IndexedVectorArray y;
        std::cout << "Measurements" << std::endl;
        std::cout << " filename: " << filename << std::endl;
        y.readVectorsFromFile(filename.c_str());
        std::cout << "Measurements loaded, size:" << y.size() << std::endl;

        IndexedVectorArray accelerometer, gyrometer, force1, force2;

        Vector yi;

        std::cout << "Splitting" << std::endl;

        for(TimeIndex i = y.getFirstIndex(); i < y.getNextIndex(); ++i)
        {
          yi = y[i];

          if(yi.size() > 0)
          {
            accelerometer.setValue(yi.segment<3>(0), i);
            gyrometer.setValue(yi.segment<3>(3), i);

            Vector f1, f2;
            int numberofContact = 0;
            if(yi.size() > 6)
            {
              numberofContact++;
              f1.resize(7);
              f1 << numberofContact, yi.segment<6>(6);
            }
            else
            {
              f1.resize(1);
              f1 << numberofContact;
            }
            if(yi.size() > 12)
            {
              f2.resize(7);
              numberofContact++;
              f1(0) = numberofContact;
              f2 << numberofContact, yi.segment<6>(12);
            }
            else
            {
              f2.resize(1);
              f2 << numberofContact;
            }

            force1.setValue(f1, i);
            force2.setValue(f2, i);
          }
        }

        std::cout << "Splitting finished" << std::endl;

        std::cout << "Writing files" << std::endl;

        accelerometer.writeInFile(filename + "-00-accelerometer");
        gyrometer.writeInFile(filename + "-01-gyrometer");
        force1.writeInFile(filename + "-02-force1");
        force2.writeInFile(filename + "-03-force2");

        std::cout << "Done " << std::endl;
      }

      else if(type == "-x")
      {
        IndexedVectorArray x;
        std::cout << "State " << std::endl;
        std::cout << " filename: " << filename << std::endl;
        x.readVectorsFromFile(filename.c_str());
        std::cout << "State loaded, size:" << x.size() << std::endl;

        std::cout << "Splitting" << std::endl;

        ///------------- state------------------------
        Vector xi;
        IndexedVectorArray pos, ori, linVel, angVel, fc1, fc2, unmodeledForces, comBias, drift;

        for(TimeIndex i = x.getFirstIndex(); i < x.getNextIndex(); ++i)
        {
          xi = x[i];
          pos.setValue(xi.segment<3>(state::pos), i);
          ori.setValue(xi.segment<3>(state::ori), i);
          linVel.setValue(xi.segment<3>(state::linVel), i);
          angVel.setValue(xi.segment<3>(state::angVel), i);
          fc1.setValue(xi.segment<6>(state::fc), i);
          fc2.setValue(xi.segment<6>(state::fc + 6), i);
          unmodeledForces.setValue(xi.segment<6>(state::unmodeledForces), i);
          comBias.setValue(xi.segment<3>(state::comBias), i);
          drift.setValue(xi.segment<3>(state::drift), i);
        }

        std::cout << "Splitting finished" << std::endl;

        std::cout << "Writing files" << std::endl;

        pos.writeInFile(filename + "-00-pos");
        ori.writeInFile(filename + "-01-ori");
        linVel.writeInFile(filename + "-02-linVel");
        angVel.writeInFile(filename + "-03-angVel");
        fc1.writeInFile(filename + "-04-fc1");
        fc2.writeInFile(filename + "-05-fc2");
        unmodeledForces.writeInFile(filename + "-06-unmodeledForces");
        comBias.writeInFile(filename + "-07-comBias");
        drift.writeInFile(filename + "-08-drift");

        std::cout << "Done " << std::endl;
      }
      else
      {
        std::cout << std::endl << "ERROR: \"" << type << "\" wrong parameter!" << std::endl << std::endl;
        throw std::invalid_argument("argc");
      }
    }
  }

  catch(const std::invalid_argument & ia)
  {
    std::cout << "USAGE: argument" << std::endl;

    std::cout << "    -u <filname>" << std::endl;
    std::cout << "    -y <filname>" << std::endl;
    std::cout << "    -x <filname>" << std::endl;
    std::cout << std::endl;
    std::cout << "Leaving" << std::endl;

    return 1;
  }

  std::cout << "Successfully finished" << std::endl;

  return 0;
}
