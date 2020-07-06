/**
 * \file      offline-ekf-flexibility-estimation.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief      Gives an implementation of flexibility estimation for IMU reconstruction
 *             with or without given input using a Kalman filter. The source is in a file
 *             imu-attitude-trajectory-reconstruction.hxx
 *
 * \details
 *
 *
 */

#ifndef FLEXIBILITYESTIMATION_OFFLINEEKFFLEXIBILITYESTIMATION_H
#define FLEXIBILITYESTIMATION_OFFLINEEKFFLEXIBILITYESTIMATION_H

#include <state-observation/flexibility-estimation/fixed-contact-ekf-flex-estimator-imu.hpp>
#include <vector>

namespace stateObservation
{
    namespace examples
    {

        /*! \fn IndexedVectorArray offlineEKFFlexibilityEstimation(
         *   const stateObservation::IndexedVectorArray & y,
         *   const stateObservation::IndexedVectorArray & u,
         *   const Matrix & xh0,
         *   unsigned numberOfContacts,
         *   const std::vector<Vector3> & contactsPositions,
         *   double dt);
         *
         *  \brief Provides the estimation of the flexibility of a robot
         *         given the measurements of an IMU and the input
         *         (which provides the reference trajectory of the IMU). This method uses
         *         extended Kalman filtering, we need to provide it with an
         *         initial guess, the number of contacts and their positions, and
         *         the time sampling period.
         *
         *
         *  \param y IMU measurements
         *  \param u the inputs of the dynamical system
         *  \param xh0 an initial guess of the state
         *  \param numberOfContacts the number of contacts
         *  \param contactsPositions a vector of positions of the vector
         *  \param dt the time discretization period
         */
        stateObservation::IndexedVectorArray offlineEKFFlexibilityEstimation(
            const stateObservation::IndexedVectorArray & y,
            const stateObservation::IndexedVectorArray & u,
            const Vector & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3, Eigen::aligned_allocator<Vector3> > & contactsPositions,
            double dt,
            IndexedVectorArray * ino=0x0,
            IndexedVectorArray * premea = 0x0);


        /*! \fn IndexedVectorArray offlineEKFFlexibilityEstimation(
         *   const stateObservation::IndexedVectorArray & y,
         *   const stateObservation::IndexedVectorArray & u,
         *   const Matrix & xh0,
         *   unsigned numberOfContacts,
         *   const std::vector<Vector3> & contactsPositions,
         *   double dt);
         *
         *  \brief Provides the estimation of the flexibility of a robot
         *         given the measurements of an IMU. This method uses
         *         extended Kalman filtering, we need to provide it with an
         *         initial guess, the number of contacts and their positions, and
         *         the time sampling period.
         *
         *
         *  \param y IMU measurements
         *  \param xh0 an initial guess of the state
         *  \param numberOfContacts the number of contacts
         *  \param contactsPositions a vector of positions of the vector
         *  \param dt the time discretization period
         */
        stateObservation::IndexedVectorArray offlineEKFFlexibilityEstimation(
            const stateObservation::IndexedVectorArray & y,
            const Vector & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3, Eigen::aligned_allocator<Vector3> > & contactsPositions,
            double dt);



# include <state-observation/examples/offline-ekf-flexibility-estimation.hxx>

    }

}

#endif // FLEXIBILITYESTIMATION_OFFLINEEKFFLEXIBILITYESTIMATION_H
