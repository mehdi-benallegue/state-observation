#include <iostream>

#include <state-observation/dynamics-estimators/kinetics-observer.hpp>
#include <state-observation/tools/probability-law-simulation.hpp>

using namespace stateObservation;

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

    unsigned testNum =0;

  
    {
        kine::Orientation ori1;
        ori1 = q;
        Matrix3 M1 = ori1;

        err +=  (Quaternion(ori1.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()  
                 + (M1 - M).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori2;
        ori2 = M;
        
        err +=  (Quaternion(ori2.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() ;
                + (ori2.matrix3() - M).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori3;
        ori3 = aa;

        err +=   (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() 
                 + (ori3.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori4;
        ori4 = Vector3(aa.angle()*aa.axis());

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()  
                + (ori4.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }


    {
        kine::Orientation ori1(q);
        Matrix3 M1 = ori1;

        err +=  (Quaternion(ori1.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()  
                 + (M1 - M).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori2(M);
        
        err +=  (Quaternion(ori2.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() ;
                + (ori2.matrix3() - M).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori3(aa);

        err +=   (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() 
                 + (ori3.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori4(Vector3(aa.angle()*aa.axis()));

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori4(q,M);

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {   
        kine::Orientation ori4;
        ori4.setRandom();
        ori4 = M;

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {   
        kine::Orientation ori4;
        ori4.setRandom();
        ori4 = q;

        err +=  (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm()  
                + (ori4.matrix3()-M).norm() ;
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }



    {
        kine::Orientation ori1(q);
        
        kine::Orientation ori2 = ori1.inverse();
        kine::Orientation ori3 = ori2.inverse();

        Matrix3 M3 = ori1;

        err +=  (Quaternion(ori3.toVector4()).toRotationMatrix() - M).norm()  
                 + (M3 - M).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
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
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        kine::Orientation ori1(q);
        kine::Orientation ori2(M);
 
        kine::Orientation ori3 = ori2.inverse()*ori1;

        err +=  (Quaternion(ori3.toVector4()).toRotationMatrix() - Matrix3::Identity()).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    {
        Vector3 v1 = Vector3::Random()*0.1;
        kine::Orientation ori1(M);
        kine::Orientation ori2(ori1);
        ori2.integrate(v1);

        Vector3 v2 = ori1.differentiate(ori2) ;
 
        err +=  (v1- v2).norm();
        std::cout << "orientation error " << ++testNum << " " << err <<std::endl;
    }

    if (err>1e-14)
    {
        return errcode;
    }
    return 0;
}

int testKinematics (int errcode)
{
    double err;
    double threshold;

    kine::Kinematics k0();

    Vector3 pos0 = Vector3::Random();

    kine::Orientation ori0;
    ori.setRandom();

    Vector3 linvel0 = Vector3::Random();
    Vector3 angvel0 = Vector3::Random();
    Vector3 linacc0 = Vector3::Random();
    Vector3 angacc0 = Vector3::Random();

    kine::Kinematics k1();

    Vector3 pos1 = Vector3::Random();

    kine::Orientation ori1;
    ori.setRandom();

    Vector3 linvel1 = Vector3::Random();
    Vector3 angvel1 = Vector3::Random();
    Vector3 linacc1 = Vector3::Random();
    Vector3 angacc1 = Vector3::Random();

    
    if (err>threshold )
    {
        return errcode;
    }
    
    return 0;


    

    

}

int main()
{

     KineticsObserver o(7);

    
    std::cout << o.getStateSize()<<std::endl;
    
    Vector x0(o.getStateSize());
    x0.setZero();
    Vector xf(x0);
    Vector xs(x0);

    int returnVal;

    if (returnVal = testOrientation(1)) /// it is not an equality check
    {
        std::cout<< "test failed, code : 1" <<std::endl;
        return returnVal;
    }

    

    std::cout<< "test succeeded" <<std::endl;
    return 0;
}