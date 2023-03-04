/// \file test_tmp.cpp
/// \brief A simple test file for surfactant module
/// \author LSEC: Song Lu

#include<iostream>
#include<surfactant/sfpde.h>
#include "phg.h"
using namespace std;
//void lsFun(double x, double y, double z, double *value)
//{
//    *value = x * x + y * y + z * z - 1.0;
//}
//
//void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
//{
//    grad[0] = x + x;
//    grad[1] = y + y;
//    grad[2] = z + z;
//}
//void IntFun(double x, double y, double z, double *ff)
//{
// *ff = 1;
//}
//int main()
//{
//    /*
//        double a[3] = {1,3,4};
//        double b[3] = {5,2,0};
//        double r[3];
//        double r1;
//        vecMinus(a,b,r);
//        for(int i=0; i<3; i++)
//        {
//            cout<<r[i]<<endl;
//        }
//        crossMul(a,b,r);
//        for(int i=0; i<3; i++)
//        {
//            cout<<r[i]<<endl;
//        }
//        r1 = dotP3(a,b);
//        cout<<r1<<endl;
//        double tet[4][3] = {{0,0,0},
//            {1,0,0},{0,1,0},{0,0,1}
//        };
//        double i = 0;
//        double x = 0;
//        double y = 0;
//        double z = 0.2;
//        double r2;
//        r2 = getBaryCoord(tet,i,x,y,z);
//        cout << r2 <<endl;
//    */
//
//
//    double tet[4][3] = {{0,0,0},
//        {0,0,1},{0,1,1},{-1,0,0}
//    };
//    double res;
//
//    int n = phgQuadInterface2(
//                lsFun,		/* the level set function */
//                2,		/* polynomial order of the level set function */
//                lsGrad,	/* the gradient of the level set function */
//                tet,		/* coordinates of the vertices of the tetra */
//                IntFun,		/* the integrand */
//                1,		/* dimension of the integrand */
//                DOF_PROJ_NONE,	/* projection type for surface integral */
//                0,		/* integration type (-1, 0, 1) */
//                10,		/* order of the 1D Gaussian quadrature */
//                &res,		/* the computed integral */
//                NULL		/* pointer returning the computed rule */
//            );
//    cout<<res<<endl;
//
//
//    exit(0);
//}





//int main()
//{
//
//    typedef DROPS::SSORPcCL SurfPcT;
////     typedef DROPS::JACPcCL SurfPcT;
//    SurfPcT surfpc;
//    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
//    SurfSolverT surfsolver( surfpc, 10000, 1e-10, true);
//    DROPS::VecDescCL xp2( &ifacep2idx);
//    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
//    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
//    DROPS::WriteFEToFile( xp2, mg, "xp2_iface.txt", /*binary=*/ false);
//    exit(0);
//}








//int main()
//{
//exit(0);
//}
