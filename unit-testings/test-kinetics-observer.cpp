#include <iostream>

#include <state-observation/dynamics-estimators/kinetics-observer.hpp>
#include <state-observation/tools/probability-law-simulation.hpp>

using namespace stateObservation;
using namespace kine;

int testOrientation(int errcode)
{   

    

    
    
    Vector4 q_i;

    tools::ProbabilityLawSimulation s;
    q_i = s.getGaussianVector(Matrix4::Identity(),Vector4::Zero(),4);
    q_i.normalize();

    Quaternion q(q_i);
    AngleAxis aa(q);
    Matrix3 M(q.toRotationMatrix());

    double err=0.;

 
    {
        kine::Orientation ori1;
        ori1 = q;
        Matrix3 M1 = ori1;

        err +=  (Quaternion(ori1.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()  
                 + (M1 - M).norm();
    }

    {
        kine::Orientation ori2;
        ori2 = M;
        
        err +=  (Quaternion(ori2.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() ;
                + (ori2.matrix3() - M).norm();
    }

    {
        kine::Orientation ori3;
        ori3 = aa;

        err +=   (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() 
                 + (ori3.matrix3()-M).norm() ;
    }

    {
        kine::Orientation ori4;
        ori4 = Vector3(aa.angle()*aa.axis());

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()  
                + (ori4.matrix3()-M).norm() ;
    }


    {
        kine::Orientation ori1(q);
        Matrix3 M1 = ori1;

        err +=  (Quaternion(ori1.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()  
                 + (M1 - M).norm();
    }

    {
        kine::Orientation ori2(M);
        
        err +=  (Quaternion(ori2.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() ;
                + (ori2.matrix3() - M).norm();
    }

    {
        kine::Orientation ori3(aa);

        err +=   (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() 
                 + (ori3.matrix3()-M).norm() ;
    }

    {
        kine::Orientation ori4(Vector3(aa.angle()*aa.axis()));

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
    }

    {
        kine::Orientation ori4(q,M);

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
    }

    {   
        kine::Orientation ori4;
        ori4.setRandom();
        ori4 = M;

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
    }

    {   
        kine::Orientation ori4;
        ori4.setRandom();
        ori4 = q;

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
    }



    {
        kine::Orientation ori1(q);
        
        kine::Orientation ori2 = ori1.inverse();
        kine::Orientation ori3 = ori2.inverse();

        Matrix3 M3 = ori1;

        err +=  (Quaternion(ori3.toVector4()).toRotationMatrix() - M).norm()  
                 + (M3 - M).norm();
    }


    {

        kine::Orientation ori0(q);
        kine::Orientation ori1(M);
        ori1=ori1.inverse();
        const kine::Orientation & ori2=ori0;
        const kine::Orientation & ori3=ori1;
        
        kine::Orientation ori00 = ori0*ori1;
        kine::Orientation ori01 = ori1*ori0;
        kine::Orientation ori02 = ori0*ori3; 
        kine::Orientation ori03 = ori3*ori0;
        kine::Orientation ori04 = ori2*ori1;
        kine::Orientation ori05 = ori1*ori2;
        kine::Orientation ori06 = ori2*ori3;
        kine::Orientation ori07 = ori3*ori2;

        kine::Orientation ori10 (ori0,ori1);
        kine::Orientation ori11 (ori1,ori0);
        kine::Orientation ori12 (ori0,ori3); 
        kine::Orientation ori13 (ori3,ori0);
        kine::Orientation ori14 (ori2,ori1);
        kine::Orientation ori15 (ori1,ori2);
        kine::Orientation ori16 (ori2,ori3);
        kine::Orientation ori17 (ori3,ori2);

        err +=  (ori00.matrix3() - Matrix3::Identity()).norm()+
                (ori01.matrix3() - Matrix3::Identity()).norm()+
                (ori02.matrix3() - Matrix3::Identity()).norm()+
                (ori03.matrix3() - Matrix3::Identity()).norm()+
                (ori04.matrix3() - Matrix3::Identity()).norm()+
                (ori05.matrix3() - Matrix3::Identity()).norm()+
                (ori06.matrix3() - Matrix3::Identity()).norm()+
                (ori07.matrix3() - Matrix3::Identity()).norm()+
                (ori10.matrix3() - Matrix3::Identity()).norm()+
                (ori11.matrix3() - Matrix3::Identity()).norm()+
                (ori12.matrix3() - Matrix3::Identity()).norm()+
                (ori13.matrix3() - Matrix3::Identity()).norm()+
                (ori14.matrix3() - Matrix3::Identity()).norm()+
                (ori15.matrix3() - Matrix3::Identity()).norm()+
                (ori16.matrix3() - Matrix3::Identity()).norm()+
                (ori17.matrix3() - Matrix3::Identity()).norm();
    }

    {
        kine::Orientation ori1(q);
        kine::Orientation ori2(M);
 
        kine::Orientation ori3 = ori2.inverse()*ori1;

        err +=  (Quaternion(ori3.toVector4()).toRotationMatrix() - Matrix3::Identity()).norm();
    }

    {
        Vector3 v1 = Vector3::Random()*0.1;
        kine::Orientation ori1(M);
        kine::Orientation ori2(ori1);
        ori2.integrate(v1);

        Vector3 v2 = ori1.differentiate(ori2) ;
 
        err +=  (v1- v2).norm();
    }

    if (err>1e-14)
    {
        return errcode;
    }
    return 0;
}

int testKinematics (int errcode)
{
    typedef kine::Kinematics::Flags Flags;


    kine::Kinematics k0;
    kine::Kinematics k1;

 
    

    Flags::Byte flag0 = BOOST_BINARY(000000);
    Flags::Byte flag1 = BOOST_BINARY(000000);


    kine::Kinematics k,l,k2,k3;

    
    int count = pow(2,6)*pow(2,6);
    double err=0;
    double threshold=1e-30*count;

    for (int i = 0; i < count; i++)
    {

        Vector3 pos0 = Vector3::Random();
        kine::Orientation ori0 = kine::Orientation::randomRotation();
        Vector3 linvel0 = Vector3::Random();
        Vector3 angvel0 = Vector3::Random();
        Vector3 linacc0 = Vector3::Random();
        Vector3 angacc0 = Vector3::Random();


        Vector3 pos1 = Vector3::Random();
        kine::Orientation ori1= kine::Orientation::randomRotation();
        Vector3 linvel1 = Vector3::Random();
        Vector3 angvel1 = Vector3::Random();
        Vector3 linacc1 = Vector3::Random();
        Vector3 angacc1 = Vector3::Random();

        k0.reset();
        k1.reset();

        if (flag0 & Flags::position)
        {
            k0.position = pos0;
        }
        if (true) ///the orientation has to be set 
        {
            k0.orientation = ori0;
        }
        if (flag0 & Flags::linVel)
        {
            k0.linVel = linvel0;
        }
        if (flag0 & Flags::angVel)
        {
            k0.angVel = angvel0;
        }
        if (flag0 & Flags::linAcc)
        {
            k0.linAcc = linacc0;
        }
        if (flag0 & Flags::angAcc)
        {
            k0.angAcc = angacc0;
        }

        if (flag1 & Flags::position)
        {
            k1.position = pos1;
        }
        if (! k1.position.isSet() || flag1 & Flags::orientation)///the orientation has to be set 
        {
            k1.orientation = ori1;
        }
        if (flag1 & Flags::linVel)
        {
            k1.linVel = linvel1;
        }
        if (flag1 & Flags::angVel)
        {
            k1.angVel = angvel1;
        }
        if (flag1 & Flags::linAcc)
        {
            k1.linAcc = linacc1;
        }
        if (flag1 & Flags::angAcc)
        {
            k1.angAcc = angacc1;
        }

        if ((flag0 = (flag0+1) & Flags::all)==0) ///update the flags to span all the possibilties
            flag1 = (flag1+1) & Flags::all;

        k2=k1;
        if (!k1.orientation.isSet())
        {
            k2.orientation.setZeroRotation();
        }

        k = ((k0 
            * k2) 
            * k2.getInverse() )
            * k0.getInverse();

        if (k.position.isSet())
        {
            err += k.position().squaredNorm();
        }
        if (k.orientation.isSet())
        {
            err += k.orientation.toRotationVector().squaredNorm();
        }
        if (k.linVel.isSet())
        {
            err +=  k.linVel().squaredNorm();
        }
        if (k.angVel.isSet())
        {
            err += k.angVel().squaredNorm();
        }
        if (k.linAcc.isSet())
        {
            err += k.linAcc().squaredNorm();
        }
        if (k.angAcc.isSet())
        {
            err += k.angAcc().squaredNorm();
        }
        
    }

     if (err>threshold )
    {
        std::cout << "Error too large : " << err <<std::endl;
        return errcode;
    }
    
    err =0;

    threshold = 1e-19*count;



    flag0 = BOOST_BINARY(000000);
    flag1 = BOOST_BINARY(000000);

     
    for (int i = 0; i < count; i++)
    {
        Vector3 pos0 = Vector3::Random();
        kine::Orientation ori0 = kine::Orientation::randomRotation();
        Vector3 linvel0 = Vector3::Random();
        Vector3 angvel0 = Vector3::Random();
        Vector3 linacc0 = Vector3::Random();
        Vector3 angacc0 = Vector3::Random();

        k.reset();
        k.position = pos0;
        k.orientation = ori0;
        k.linVel = linvel0;
        k.angVel = angvel0;
        k.linAcc = linacc0;
        k.angAcc = angacc0;

        l = k;

        double dt = 0.001;
        k.integrate(dt);


        k0.reset();
        k1.reset();

        if (i==-1)
        {
            std::cout <<"debug me ";
        }

        if (flag0 & Flags::position)
        {
            k0.position = l.position;
        }
        if (flag0 & Flags::orientation)  
        {
            k0.orientation = l.orientation;
        }
        if (flag0 & Flags::linVel)
        {
            k0.linVel = l.linVel;
        }
        if (flag0 & Flags::angVel)
        {
            k0.angVel = l.angVel;
        }
        if (flag0 & Flags::linAcc)
        {
            k0.linAcc = l.linAcc;
        }
        if (flag0 & Flags::angAcc)
        {
            k0.angAcc = l.angAcc;
        }

        if (flag1 & Flags::position)
        {
            k1.position = k.position;
        }
        if (flag1 & Flags::orientation)
        {
            k1.orientation = k.orientation;
        }
        if (flag1 & Flags::linVel)
        {
            k1.linVel = k.linVel;
        }
        if (flag1 & Flags::angVel)
        {
            k1.angVel = k.angVel;
        }
        if (flag1 & Flags::linAcc)
        {
            k1.linAcc = k.linAcc;
        }
        if (flag1 & Flags::angAcc)
        {
            k1.angAcc = k.angAcc;
        }


        if ((flag0 = (flag0+1) & Flags::all)==0) ///update the flags to span all the possibilties
            flag1 = (flag1+1) & Flags::all;

        kine::Kinematics::Flags::Byte flag =BOOST_BINARY(000000);

        if (k1.position.isSet() || 
            (k0.position.isSet() && k0.linVel.isSet() && k0.linAcc.isSet())
           )
        {
            flag = flag | kine::Kinematics::Flags::position;
        }
        if (k1.orientation.isSet() || 
            (k0.orientation.isSet() && k0.angVel.isSet() && k0.angAcc.isSet())
            )
        {
            flag = flag | kine::Kinematics::Flags::orientation;
        }
        if (k1.linVel.isSet() || 
            (k0.position.isSet() && k1.position.isSet()) || 
            (k0.linAcc.isSet() && k0.linVel.isSet())
            )
        {
            flag = flag | kine::Kinematics::Flags::linVel;
        }
        if (k1.angVel.isSet() ||
            (k0.orientation.isSet() && k1.orientation.isSet()) || 
            (k0.angVel.isSet() && k0.angAcc.isSet())
           )
        {
            flag = flag | kine::Kinematics::Flags::angVel;
        }
        if (k1.linAcc.isSet() ||
             (k0.linVel.isSet() && k1.linVel.isSet())||
             (k0.linVel.isSet() && k0.position.isSet() && k1.position.isSet())
            )
        {
            flag = flag | kine::Kinematics::Flags::linAcc;
        }
        if (k1.angAcc.isSet()||
            (k0.angVel.isSet() && k1.angVel.isSet()) ||
            (k0.angVel.isSet() && k0.orientation.isSet() && k1.orientation.isSet())
            )
        {
            flag = flag | kine::Kinematics::Flags::angAcc;
        }

        k3 = k0;
        k0.update(k1,dt,flag);

        if (k0.position.isSet())
        {
            err += (k.position()-k0.position()).squaredNorm();
        }
        if (k0.orientation.isSet())
        {
            err += (k.orientation.differentiate(k0.orientation)).squaredNorm();
        }
        if (k0.linVel.isSet())
        {
            if ((k.linVel()-k0.linVel()).squaredNorm()<1e-8)
            {
                err +=  (k.linVel()-k0.linVel()).squaredNorm();
            }
            else
            {
                err += ((k.position()-l.position())/dt-k0.linVel()).squaredNorm();
            }
            
            
        }
        if (k0.angVel.isSet())
        {
            if ((k.angVel()-k0.angVel()).squaredNorm()<1e-9)
            {
                err +=  (k.angVel()-k0.angVel()).squaredNorm();
            }
            else
            {
                err += (l.orientation.differentiate(k.orientation)/dt-k0.angVel()).squaredNorm();
            }
        }
        if (k0.linAcc.isSet())
        {
            if ((k.linAcc()-k0.linAcc()).squaredNorm() < 1e-5)
            {
                err +=(k.linAcc()-k0.linAcc()).squaredNorm() ;
            }
            else
            {
                err +=(k.linAcc()-2*k0.linAcc()).squaredNorm() ;
            }            
        }
        if (k0.angAcc.isSet())
        {            
            if ((k.angAcc()-k0.angAcc()).squaredNorm()< 1e-5)
            {
                err +=(k.angAcc()-k0.angAcc()).squaredNorm() ;
            }
            else
            {
                 err +=(k.angAcc()-2*k0.angAcc()).squaredNorm() ;
            }            
        }

//        std::cout<< i<<" "<<err << std::endl;

        if (err>threshold )
        {
            break;
        }
        
    }

    
    if (err>threshold )
    {
        std::cout << "Error too large : " << err <<std::endl;
        return errcode;
    }
    
    return 0;


    

    

}

int testKineticsObserver(int code)
{
    double dt = 0.001;
    KineticsObserver o(6);

    Vector x0(o.getStateSize());
    x0.setZero();
    Vector xf(x0);
    Vector xs(x0);

    o.setSamplingTime(dt);

    Kinematics stateKine;

    stateKine.position.set() << 0,0,0.7;
    stateKine.orientation.setZeroRotation();
    stateKine.linVel.set().setZero();
    stateKine.angVel.set().setZero();


    o.setStateKinematics(stateKine);

    Vector x = o.getStateVector();
    stateObservation::TimeIndex index = o.getStateVectorSampleTime();

    std::cout << index << " " << x.transpose() << std::endl;

    o.update();

    Kinematics k = o.getKinematics();

    

    std::cout << k;

    Kinematics l = o.getKinematicsOf(k);

    std::cout << l;

    

     


    return 0;
}

int main()
{
    int returnVal;

    if ((returnVal = testOrientation(1))) /// it is not an equality check
    {
        std::cout<< "Orientation test failed, code : 1" <<std::endl;
        return returnVal;
    }
    else
    {
        std::cout << "Orientation test succeeded" << std::endl;
    }
    


    if ((returnVal = testKinematics(2))) /// it is not an equality check
    {
        std::cout<< "Kinematics test failed, code : 1" <<std::endl;
        return returnVal;
    }
    else
    {
        std::cout << "Kinematics test succeeded" << std::endl;
    }

    if ((returnVal = testKineticsObserver(3)))
    {
        std::cout<< "Kinetics Observer test failed, code : 3" <<std::endl;
        return returnVal;
    }
    else
    {
        std::cout << "Kinetics Observer test succeeded" << std::endl;
    }


    

    std::cout<< "test succeeded" <<std::endl;
    return 0;
}