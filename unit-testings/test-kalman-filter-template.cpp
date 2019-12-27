#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/random.hpp>

#include <algorithm>

#include <bitset>
#include <iomanip>

#include <state-observation/observer/compile-time/compile-time-kalman-filter.hpp>
#include <state-observation/observer/compile-time/compile-time-extended-kalman-filter.hpp>

#include <boost/utility/binary.hpp>

#include <state-observation/tools/probability-law-simulation.hpp>


double testExtendedKalmanFilter()
{
    ///the number of samples
    const static unsigned kmax=1000;

    ///define the type of the extended Kalman filter
    typedef stateObservation::compileTime::ExtendedKalmanFilter<4,3,1> ekf;

    ///The functor that describes the dynamics of the state
    ///and the measurement
    class KalmanFunctor:
                public ekf::DynamicsFunctorBase
    {

    public:
        ///Constructor
        KalmanFunctor()
        {
            s_=ekf::StateVector::Random()*0.1;
            m_=ekf::MeasureVector::Random()*0.1;
            t_=ekf::StateVector::Random();
            n_=ekf::MeasureVector::Random();
            a_=ekf::Amatrix::Random()*0.6;
            c_=ekf::Cmatrix::Random();
        }

        ///The dynamics of the state xk1=f(xk,u,k)
        virtual ekf::StateVector stateDynamics(const ekf::StateVector& xk, const ekf::InputVector& u, unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::StateVector xk1;
            xk1=a_*xk+cos(10*(xk.transpose()*xk)[0])*s_+t_+(u*u.transpose())[0]*s_;
            return xk1;
        }

        ///The dynamics of the state yk=h(xk,u,k)
        virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& xk, const ekf::InputVector &u,unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::MeasureVector yk;
            yk=c_*xk+cos(10*(xk.transpose()*xk)[0])*m_+n_;
            return yk;
        }

    private:
        ///Special instructions to have a static-sized eigen vector as a member
        enum { NeedsToAlign = ((sizeof(ekf::StateVector)%16)==0)||
                              ((sizeof(ekf::MeasureVector)%16)==0)||
                              ((sizeof(ekf::Amatrix)%16)==0)||
                              ((sizeof(ekf::Cmatrix)%16)==0)
             };
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        ///containers for the vectors and matrices
        ekf::StateVector s_;
        ekf::StateVector t_;
        ekf::MeasureVector m_;
        ekf::MeasureVector n_;

        ekf::Amatrix a_;
        ekf::Cmatrix c_;
    };

    ///containers for the state, the measurements and the input
    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector uk[kmax+1];

     ///the standard deviation matrix to generate the gaussian noise
     ekf::Rmatrix r1=ekf::Rmatrix::Random()*0.01;
     ekf::Qmatrix q1=ekf::Qmatrix::Random()*0.01;

     ///instanciate the functor
     KalmanFunctor func;

    { /// Construction of the sequence of states measurements and inputs

        ///initializations
        ekf::StateVector x=ekf::StateVector::Zero();
        xk[0]=x;
        uk[0]=ekf::InputVector::Random();

        for (unsigned k=1; k<=kmax; ++k)
        {
            ///generation of random inputs
            uk[k]=ekf::InputVector::Random();

            ///generation of Gaussian white noises
            ekf::StateVector v= stateObservation::tools::ProbabilityLawSimulation::
                                getGaussianVector(q1,ekf::StateVector::Zero(),ekf::stateSize);

            ekf::MeasureVector w= stateObservation::tools::ProbabilityLawSimulation::
                                  getGaussianVector(r1,ekf::MeasureVector::Zero(),ekf::measureSize);

            ///the dynamics is executed here
            xk[k]=x=func.stateDynamics(x,uk[k-1],k-1)+v;
            yk[k-1]=func.measureDynamics(x,uk[k],k)+w;
        }
    }

    ///instanciation of the extended Kalman filter
    ekf f;

    ///set the functor of the extended Kalman filter
    f.setFunctor(&func);

    ///generation of a random initial estimation of the state
    ekf::StateVector xh=ekf::StateVector::Random();

    ///set the initial state of the estimator
    f.setState(xk[0],0);

    ///set the covariance matrix of the initial estimation error
    ekf::Pmatrix p=ekf::Pmatrix::Zero();
    for (unsigned i=0;i<ekf::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();
    f.setStateCovariance(p);

    ///the covariance matrices for the process noise and the measurements noise
    ekf::Rmatrix r(r1*r1.transpose());
    ekf::Qmatrix q(q1*q1.transpose());

    ///set the covariance matrices for the extended Kalman filter
    f.setR(r);
    f.setQ(q);

    ///set initial input
    f.setInput(uk[0],0);

    ///set the derivation step for the finite difference method
    ekf::StateVector dx=ekf::StateVector::Constant(1)*1e-8;


    unsigned i;
    for (i=1;i<=kmax;++i)
    {
        ///give the measurements and the inputs at instant i to the ekf
        f.setMeasurement(yk[i-1],i);
        f.setInput(uk[i],i);

        ///obtain jacobians by finite differences method
        ekf::Amatrix a=f.getAMatrixFD(dx);
        ekf::Cmatrix c= f.getCMatrixFD(dx);

        ///set the jacobians to the ekf
        f.setA(a);
        f.setC(c);

        ///get the estimation of the state at instant i;
        xh=f.getEstimatedState(i);
     }

    ekf::StateVector error=xh-xk[kmax];

    return error.norm();
}

double testExtendedKalmanFilterLTV()
{
    const static unsigned kmax=1000;

    typedef stateObservation::compileTime::ExtendedKalmanFilter<4,3,1> ekf;


    struct KalmanFunctorLTV:
                public ekf::DynamicsFunctorBase
    {

public:
        KalmanFunctorLTV()
        {
            s_=ekf::StateVector::Random();
            n_=ekf::MeasureVector::Random();

            for (unsigned i=0;i<=kmax;++i)
            {
                a.push_back( ekf::Amatrix::Random()*0.5 );
                c.push_back( ekf::Cmatrix::Random() );
            }

        }

        virtual ekf::StateVector stateDynamics(const ekf::StateVector& x, const ekf::InputVector& u, unsigned k)
        {
            ekf::StateVector xk1;
            unsigned kk=std::min(k,kmax);
            xk1=a[kk]*x+(u*u.transpose())[0]*s_;

            return xk1;
        }

        virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& x, const ekf::InputVector &u,unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::MeasureVector yk;
            unsigned kk=std::min(k,kmax);
            yk=c[kk]*x+(u*u.transpose())[0]*n_;
            return yk;
        }

        std::vector<ekf::Amatrix,Eigen::aligned_allocator<ekf::Amatrix> > a;
        std::vector<ekf::Cmatrix,Eigen::aligned_allocator<ekf::Cmatrix> > c;
private:
        ekf::StateVector s_;
        ekf::MeasureVector n_;
    };


    ekf f;

    ekf::Amatrix a;
    KalmanFunctorLTV func;

    f.setFunctor(&func);



    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector uk[kmax+1];

    ekf::StateVector x=ekf::StateVector::Zero();



    boost::lagged_fibonacci1279 gen_;

    ekf::Rmatrix r1=ekf::Rmatrix::Random()*0.01;

    ekf::Qmatrix q1=ekf::Qmatrix::Random()*0.01;

    ekf::Rmatrix r(r1*r1.transpose());
    ekf::Qmatrix q(q1*q1.transpose());

    xk[0]=x;
    uk[0]=ekf::InputVector::Random();

    for (unsigned k=1; k<=kmax; ++k)
    {
        ekf::StateVector v;
        for (unsigned i=0;i<ekf::stateSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            v[i]=g(gen_);
        }
        v=q1*v;

        ekf::MeasureVector w;
        for (unsigned i=0;i<ekf::measureSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            w[i]=g(gen_);
        }
        w=r1*w;

        uk[k]=ekf::InputVector::Random();

        x=func.stateDynamics(x,uk[k-1],k-1)+v;

        xk[k]=x;
        yk[k-1]=func.measureDynamics(x,uk[k],k)+w;

    }

    ekf::StateVector xh=ekf::StateVector::Random();

    f.setState(xh,0);

    ekf::Pmatrix p=ekf::Pmatrix::Zero();



    f.setStateCovariance(p);

    f.setR(r);
    f.setQ(q);

    f.setInput(uk[0],0);

    for (unsigned i=0;i<ekf::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    ekf::StateVector dx=ekf::StateVector::Constant(1)*1e-8;

    unsigned i;
    for (i=1;i<=kmax;++i)
    {
        f.setMeasurement(yk[i-1],i);
        f.setInput(uk[i],i);


        ekf::Amatrix a=f.getAMatrixFD(dx);
        ekf::Cmatrix c= f.getCMatrixFD(dx);

        f.setA(a);
        f.setC(c);

        xh=f.getEstimatedState(i);
    }

    ekf::StateVector error(xh-xk[kmax]);

    return error.norm();

}

double testExtendedKalmanFilterZeroInput()
{
    const static unsigned kmax=1000;


    typedef stateObservation::compileTime::ExtendedKalmanFilter<4,3> ekf;


    class KalmanFunctor:
                public ekf::DynamicsFunctorBase
    {

    public:
        KalmanFunctor()
        {
            m_=ekf::MeasureVector::Random()*0.1;
            s_=ekf::StateVector::Random()*0.1;
            t_=ekf::StateVector::Random();
            n_=ekf::MeasureVector::Random();
        }

        virtual ekf::StateVector stateDynamics(const ekf::StateVector& x, const ekf::InputVector& u, unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::StateVector xk1;
            xk1=a_*x+cos(10*(x.transpose()*x)[0])*s_+t_;
            return xk1;
        }

        virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& x, const ekf::InputVector &u,unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::MeasureVector yk;
            yk=c_*x+cos(10*(x.transpose()*x)[0])*m_+n_;
            return yk;
        }

        void setA(const ekf::Amatrix& a)
        {
            a_=a;
        }

        void setC(const ekf::Cmatrix& c)
        {
            c_=c;
        }

    private:
        ///Special instructions to have a static-sized eigen vector as a member
        enum { NeedsToAlign = ((sizeof(ekf::StateVector)%16)==0)||
                              ((sizeof(ekf::MeasureVector)%16)==0)||
                              ((sizeof(ekf::Amatrix)%16)==0)||
                              ((sizeof(ekf::Cmatrix)%16)==0)
             };
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        ekf::StateVector s_;
        ekf::StateVector t_;
        ekf::MeasureVector m_;
        ekf::MeasureVector n_;

        ekf::Amatrix a_;
        ekf::Cmatrix c_;

    };


    ekf f;

    KalmanFunctor func;

    func.setA(ekf::Amatrix::Random()*0.5);
    func.setC(ekf::Cmatrix::Random());



    f.setFunctor(&func);


    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector u= ekf::InputVector::Zero();

    ekf::StateVector x=ekf::StateVector::Zero();

    boost::lagged_fibonacci1279 gen_;

    ekf::Rmatrix r1=ekf::Rmatrix::Random()*0.01;

    ekf::Qmatrix q1=ekf::Qmatrix::Random()*0.01;

    ekf::Rmatrix r(r1*r1.transpose());
    ekf::Qmatrix q(q1*q1.transpose());

    xk[0]=x;


    for (unsigned k=1; k<=kmax; ++k)
    {
        ekf::StateVector v;
        for (unsigned i=0;i<ekf::stateSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            v[i]=g(gen_);
        }
        v=q1*v;

        ekf::MeasureVector w;
        for (unsigned i=0;i<ekf::measureSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            w[i]=g(gen_);
        }
        w=r1*w;

        x=func.stateDynamics(x,u,k-1)+v;

        xk[k]=x;
        yk[k-1]=func.measureDynamics(x,u,k)+w;

    }

    ekf::StateVector xh=ekf::StateVector::Random();

    f.setState(xk[0],0);

    ekf::Pmatrix p=ekf::Pmatrix::Zero();

    f.setStateCovariance(p);

    f.setR(r);
    f.setQ(q);

    for (unsigned i=0;i<ekf::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    ekf::StateVector dx=ekf::StateVector::Constant(1)*1e-8;

    unsigned i;
    for (i=1;i<=kmax;++i)
    {
        f.setMeasurement(yk[i-1],i);


        ekf::Amatrix a=f.getAMatrixFD(dx);
        ekf::Cmatrix c= f.getCMatrixFD(dx);

        f.setA(a);
        f.setC(c);
        xh=f.getEstimatedState(i);

    }

    ekf::StateVector error=xh-xk[kmax];

    return error.norm();

}


double testKalmanFilter()
{

    typedef stateObservation::compileTime::KalmanFilter<4,3,2> filter;

    filter f;
    filter::Amatrix a;

    a<<	-0.6785714 ,0.1156463 ,	0.4392517 ,	0.2863946 ,
    0.0865306  ,-0.0273810,	0.3355102 ,	0.0184150 ,
    - 0.4172789,-0.2036735,	-0.4434014,	-0.2666667,
    0.4200680  ,0.5387075 ,	0.4883673 , -0.6598639;

    filter::Bmatrix b=filter::Bmatrix::Random();
    filter::Cmatrix c=filter::Cmatrix::Random();
    filter::Dmatrix d=filter::Dmatrix::Random();

    const unsigned kmax=1000;

    filter::StateVector xk[kmax+1];
    filter::MeasureVector yk[kmax];
    filter::InputVector uk[kmax+1];

    filter::StateVector x=filter::StateVector::Zero();

    xk[0]=x;



    filter::Rmatrix r1=filter::Rmatrix::Random()*0.01;

    filter::Qmatrix q1=filter::Qmatrix::Random()*0.01;

    filter::Rmatrix r(r1*r1.transpose());
    filter::Qmatrix q(q1*q1.transpose());

    uk[0]=filter::InputVector::Random();

    boost::lagged_fibonacci1279 gen_;
    for (unsigned k=1; k<=kmax; ++k)
    {

        filter::StateVector v;
        for (unsigned i=0;i<filter::stateSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            v[i]=g(gen_);
        }
        v=q1*v;

        filter::MeasureVector w;
        for (unsigned i=0;i<filter::measureSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            w[i]=g(gen_);
        }
        w=r1*w;

        uk[k]=filter::InputVector::Random();




        xk[k]=x=a*x+v+b*uk[k-1];
        yk[k-1]=c*x+w+d*uk[k];
    }

    filter::StateVector xh=filter::StateVector::Random();

    f.setState(xh,0);

    filter::Pmatrix p=filter::Pmatrix::Zero();

    for (unsigned i=0;i<filter::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    f.setStateCovariance(p);

    f.setA(a);
    f.setB(b);
    f.setC(c);
    f.setD(d);

    f.setR(r);
    f.setQ(q);

    f.setInput(uk[0],0);

    unsigned i;
    for (i=1;i<=kmax;++i)
    {
        f.setMeasurement(yk[i-1],i);
        f.setInput(uk[i],i);
    }

    filter::StateVector error=f.getEstimatedState(kmax)-xk[kmax];

    return error.norm();
}

double testKalmanFilterZeroInput()
{
    typedef stateObservation::compileTime::KalmanFilter<4,3> filter;

    filter f;
    filter::Amatrix a;

    a<<	-0.6785714 ,0.1156463 ,	0.4392517 ,	0.2863946 ,
    0.0865306  ,-0.0273810,	0.3355102 ,	0.0184150 ,
    - 0.4172789,-0.2036735,	-0.4434014,	-0.2666667,
    0.4200680  ,0.5387075 ,	0.4883673 , -0.598639;

    filter::Bmatrix b=filter::Bmatrix::Random();
    filter::Cmatrix c=filter::Cmatrix::Random();
    filter::Dmatrix d=filter::Dmatrix::Random();

    const unsigned kmax=1000;

    filter::StateVector xk[kmax+1];
    filter::MeasureVector yk[kmax];

    filter::StateVector x=filter::StateVector::Zero();

    xk[0]=x;

    boost::lagged_fibonacci1279 gen_;

    filter::Rmatrix r1=filter::Rmatrix::Random()*0.01;

    filter::Qmatrix q1=filter::Qmatrix::Random()*0.01;

    filter::Rmatrix r(r1*r1.transpose());
    filter::Qmatrix q(q1*q1.transpose());

    for (unsigned k=1; k<=kmax; ++k)
    {

        filter::StateVector v;
        for (unsigned i=0;i<filter::stateSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            v[i]=g(gen_);
        }
        v=q1*v;

        filter::MeasureVector w;
        for (unsigned i=0;i<filter::measureSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            w[i]=g(gen_);
        }
        w=r1*w;

        xk[k]=x=a*x+v;
        yk[k-1]=c*x+w;
    }

    filter::StateVector xh=filter::StateVector::Random();

    f.setState(xh,0);

    filter::Pmatrix p=filter::Pmatrix::Zero();

    for (unsigned i=0;i<filter::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    f.setStateCovariance(p);

    f.setA(a);
    f.setB(b);
    f.setC(c);
    f.setD(d);

    f.setR(r);
    f.setQ(q);

    unsigned i;

    for (i=1;i<=kmax;++i)
    {
        f.setMeasurement(yk[i-1],i);
    }

    filter::StateVector error=f.getEstimatedState(kmax)-xk[kmax];

    return error.norm();
}

int main()
{
    short exit=0;
    double error;
    std::cout<<"Starting"<<std::endl;

    if ((error=testKalmanFilter())<0.1)
    {
        std::cout<<"Test Kalman filter SUCCEEDED: estimationError = "<<error<<std::endl;
    }
    else
    {
        exit=exit | BOOST_BINARY( 1 );
        std::cout<<"Test Kalman filter FAILED: estimationError = "<<error<<std::endl;
    }
    if ((error=testKalmanFilterZeroInput())<0.1)
    {
        std::cout<<"Test Kalman filter (zero input) SUCCEEDED: estimationError = "<<error<<std::endl;
    }
    else
    {
        exit=exit | BOOST_BINARY( 10 );
        std::cout<<"Test Kalman filter (zero input) FAILED: estimationError = "<<error<<std::endl;
    }
    if ((error=testExtendedKalmanFilter())<0.1)
    {
        std::cout<<"Test extended Kalman filter SUCCEEDED: estimationError = "<<error<<std::endl;
    }
    else
    {
        exit=exit | BOOST_BINARY( 100 );
        std::cout<<"Test extended Kalman filter FAILED: estimationError = "<<error<<std::endl;
    }
    if ((error=testExtendedKalmanFilterLTV())<0.1)
    {
        std::cout<<"Test extended Kalman filter (LTV) SUCCEEDED: estimationError = "<<error<<std::endl;
    }
    else
    {
        exit=exit | BOOST_BINARY( 1000 );
        std::cout<<"Test extended Kalman filter (LTV) FAILED: estimationError = "<<error<<std::endl;
    }
    if ((error=testExtendedKalmanFilterZeroInput())<0.1)
    {
        std::cout<<"Test extended Kalman filter (zero input) SUCCEEDED: estimationError = "<<error<<std::endl;
    }
    else
    {
        exit=exit | BOOST_BINARY( 10000 );
        std::cout<<"Test extended Kalman filter (zero input) FAILED: estimationError = "<<error<<std::endl<<std::endl;
    }


    std::cout<<"Test exit code "<< std::bitset< 16 >(exit) <<std::endl;

    return exit;

}
