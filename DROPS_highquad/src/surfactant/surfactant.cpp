/// \file
/// \brief Solve a non-stationary convection-diffusion-equation on a moving interface
/// \author LNM RWTH Aachen: Joerg Grande, LSEC: Song Lu

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "levelset/levelsetmapper.h"
#include "out/output.h"
#include "out/vtkOut.h"
#include "misc/dynamicload.h"
#include "misc/funcmap.h"
#include "misc/omp_variable.h"
#include "geom/subtriangulation.h"
#include "num/gradient_recovery.h"
#include "surfactant/sfpde.h"
/*
#include "surfphasesep/separation.
*/
#include <cmath>
#include <fstream>
#include <string>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <surfactant/sfpde.h>
#include <surfactant/femP3.h>
//#include <phg.h>
#include "num/directsolver.h"

using namespace DROPS;

DROPS::ParamCL P;

std::unique_ptr<VTKOutCL> vtkwriter;

DROPS::InVecMap& invecmap= DROPS::InVecMap::getInstance();
DROPS::InScaMap& inscamap= DROPS::InScaMap::getInstance();

instat_vector_fun_ptr the_wind_fun;
instat_scalar_fun_ptr the_lset_fun;
instat_vector_fun_ptr the_normal_fun;//only for surfactant extension method

instat_scalar_fun_ptr the_rhs_fun;
instat_scalar_fun_ptr the_sol_fun;
instat_vector_fun_ptr the_sol_grad_fun;

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc_wind[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
bnd_val_fun bf_wind[6];

instat_scalar_fun_ptr sigma( 0);
SurfaceTensionCL sf( sigma, 0);
DROPS::LsetBndDataCL lsbnd( 6);


// Surface divergence of a vector field w
inline double div_gamma_wind (const Point3DCL& n, const SMatrixCL<3,3>& dw)
{
    return trace( dw) - inner_prod( n, dw*n);
}

// laplace-beltrami of a function u
inline double laplace_beltrami_u (const Point3DCL& n,      const SMatrixCL<3,3>& dn,
                                  const Point3DCL& grad_u, const SMatrixCL<3,3>& Hess_u)
{
    const double tr_PHessu= trace( Hess_u) - inner_prod( n, Hess_u*n),
                 tr_Pdn= trace( dn) - inner_prod( n, dn*n),
                 ngradu= inner_prod( n, grad_u);
    return tr_PHessu - tr_Pdn*ngradu;
}


DROPS::Point3DCL WindVelocity;
DROPS::Point3DCL constant_wind (const DROPS::Point3DCL&, double)
{
    return WindVelocity;
}
static RegisterVectorFunction regvec_constant_wind( "ConstantWind", constant_wind);

#if 1
DROPS::Point3DCL RadDrop;
DROPS::Point3DCL PosDrop;
#endif

double ellipsoid (const DROPS::Point3DCL& p, double)
{
    const DROPS::Point3DCL x= (p - PosDrop)/RadDrop;
    return x.norm_sq() - 1.;
}
static RegisterScalarFunction regsca_ellipsoid_lset( "Ellipsoid", ellipsoid);


DROPS::Point2DCL RadTorus; // R= RadTorus[0], r= RadTorus[1]; R > r for tori with a hole.
double torus (const Point3DCL& p, double)
{
    return std::sqrt( std::pow( RadTorus[0] - std::sqrt(p[0]*p[0] + p[1]*p[1]), 2) + std::pow( p[2], 2)) - RadTorus[1];
}
static RegisterScalarFunction regsca_torus_lset( "Torus", torus);

// ==non-stationary test case "HeatConduction"
// "ConstantWind" 0
// "Ellipsoid" pos:[0,0,0], rad:[1,1,1]
double heat_conduction_u0 (const Point3DCL& p)
{
    //  return 1+p[0]*p[1]*p[2]/pow(p.norm(),3);//Ex0 and new Ex00
    //  return 2*p[0]*p[1]*p[2]/pow(p.norm(),3);//Ex0-1
    //     return 1; //Ex0-ode
    //   return 1;//new Ex0SurfTransp
    return 1+p[0]+p[1]+p[2]; //new Ex000
//	double a=100;
//	return 1*a+(p[0]+p[1]+p[2])/p.norm()/a;//new Ex000-1; new Ex000-2
}

double heat_conduction_rhs (const Point3DCL& p, double t)
{
    return 0;//new Ex000//new Ex000-constant wind
//	double a=100;
//	return (p[0]+p[1]+p[2])*exp(-t)/p.norm()/a;//new Ex000-1
//	return 2*(p[0]+p[1]+p[2])/p.norm()/a-a*exp(-t); //new Ex000-2
    //   return 11*p[0]*p[1]*p[2]*exp(-t);//Ex0
    //     return 11*p[0]*p[1]*p[2]*exp(-t)+12*p[0]*p[1]*p[2];//Ex0-1
//	return -exp(-t);//Ex0-ode

//	double D=P.get<double>("SurfTransp.Visc");
//	double m=10;
//	return m*cos(m*t);//Ex0-ode1
//	return p[0]/p.norm()*p[1]/p.norm()*p[2]/p.norm()*(12*D*sin(m*t)+m*cos(m*t));//new Ex0
//        return p[0]*p[1]*p[2]*(12*D+m*cos(m*t))*std::exp(sin(m*t));//new Ex00
}
static RegisterScalarFunction regsca_heat_conduction_rhs( "HeatConductionRhs", heat_conduction_rhs);

double heat_conduction_sol (const Point3DCL& p, double t)
{
    double D=P.get<double>("SurfTransp.Visc");
    // return std::exp( -6.*t)*heat_conduction_u0( p);
    //    return 1+std::exp(-t)*p[0]*p[1]*p[2];//Ex0
//	return (1+std::exp(-t))*p[0]*p[1]*p[2];//Ex0-1
//	return exp(-t);//Ex0-odeSurfTransp
    //      return std::exp(-2.0*t)*(p[0]+p[1]+p[2])+1;//new Ex000
//	double a=100;
//	return exp(-t)*(p[0]+p[1]+p[2])/p.norm()/a+1;//new Ex000-1
//	return a*exp(-t)+(p[0]+p[1]+p[2])/p.norm()/a;//new Ex000-2
//	double m=10;
//	return 1+sin(m*t);//Ex0-ode1
//	return 1+p[0]/p.norm()*p[1]/p.norm()*p[2]/p.norm()*sin(m*t);//new Ex0
//        return 1+p[0]*p[1]*p[2]*std::exp(sin(m*t));//mew Ex00
    DROPS::Point3DCL x= (p - PosDrop - t*constant_wind(p, t));
    return std::exp(-2.0*D*t)*(x[0]+x[1]+x[2])+1;//new Ex000-- constant wind
}
static RegisterScalarFunction regsca_heat_conduction_sol( "HeatConductionSol", heat_conduction_sol);

DROPS::Point3DCL heat_conduction_surfgradsol (const DROPS::Point3DCL& p, double t)
{
    Point3DCL surfgrad;
    double D=P.get<double>("SurfTransp.Visc");
//	double a=100;
    /*	double b=(p[0]+p[1]+p[2])/p.norm();
    	surfgrad[0]=1-b*p[0]/p.norm();
    	surfgrad[1]=1-b*p[1]/p.norm();
    	surfgrad[2]=1-b*p[2]/p.norm();
            surfgrad*=std::exp(-2.*t);*///new Ex000
//	surfgrad*=exp(-t)/a;  //new Ex000-1
//	surfgrad/=a;
//	double nm=p.norm();
//	double a=3*p[0]/nm*p[1]/nm*p[2]/nm;
//	surfgrad[0]=p[1]/nm*p[2]/nm-a*p[0]/nm;
//	surfgrad[1]=p[0]/nm*p[2]/nm-a*p[1]/nm;
//	surfgrad[2]=p[0]/nm*p[1]/nm-a*p[2]/nm;//*/
//	double m=10;
//	surfgrad*=sin(m*t); //new Ex0
//        surfgrad*=std::exp(sin(m*t));//new Ex00
//	surfgrad*=std::exp(-t); //Ex0
//	surfgrad*=(1+std::exp(-t)); //Ex0-1
//	surfgrad*=0; //Ex0-ode

    DROPS::Point3DCL x= (p - PosDrop - t*constant_wind(p, t));
    double b=(x[0]+x[1]+x[2])/x.norm();
    surfgrad[0]=1-b*x[0]/x.norm();
    surfgrad[1]=1-b*x[1]/x.norm();
    surfgrad[2]=1-b*x[2]/x.norm();
    surfgrad*=std::exp(-2.*D*t);///new Ex000--constant wind
    return surfgrad;
}
static RegisterVectorFunction regvec_heat_conduction_surfgradsol( "HeatConductionSurfGradSol", heat_conduction_surfgradsol);

// ==non-stationary test-case "ToroidalFlow"==
// Level set: "torus" with RadTorus

double angular_velocity (const Point3DCL& p, double)
{
    return 1. + p[2];
}

DROPS::SMatrixCL<3,3> rotation_matrix (const DROPS::Point3DCL& p, double t)
{
    const double omega= angular_velocity( p, 0.);
    SMatrixCL<3,3> m;
    m(0,0)= std::cos(omega*t);
    m(0,1)= -std::sin(omega*t);
    m(1,0)= std::sin(omega*t);
    m(1,1)=  std::cos(omega*t);
    m(2,2)= 1.;
    return m;
}

DROPS::Point3DCL toroidal_flow (const DROPS::Point3DCL& p, double t)
{
    return rotation_matrix( p, t)*p;
}

DROPS::Point3DCL toroidal_flow_wind (const DROPS::Point3DCL& p, double)
{
    const double omega= angular_velocity( p, 0.);
    return MakePoint3D( -p[1]*omega, p[0]*omega, 0.);
}
static RegisterVectorFunction regvec_toroidal_flow_wind( "ToroidalFlowWind", toroidal_flow_wind);

double toroidal_flow_sol (const Point3DCL& p, double t)
{
    const Point3DCL q( toroidal_flow( p, -t));
    return std::exp( -t)*q[0]*q[1];
}
static RegisterScalarFunction regsca_toroidal_flow_sol( "ToroidalFlowSol", toroidal_flow_sol);

double toroidal_flow_rhs (const Point3DCL& p, double t)
{
    const Point3DCL w( toroidal_flow_wind( p, t));
    const double omega= angular_velocity( p, 0.);
    SMatrixCL<3,3> dw;
    dw(0,1)= -omega;
    dw(0,2)= -p[1];
    dw(1,0)=  omega;
    dw(1,2)=  p[0];

    const Point2DCL xhat( MakePoint2D( p[0], p[1]));
    const double norm_xhat= xhat.norm();
    const double l= torus( p, t) + RadTorus[1];
    const Point2DCL tt= (norm_xhat - RadTorus[0])/(l*norm_xhat)*xhat;
    const Point3DCL n( MakePoint3D( tt[0], tt[1], p[2]/l));
    SMatrixCL<3,3> dn;
    dn= eye<3,3>() - outer_product( n, n);
    SMatrixCL<2,2> dnhat= RadTorus[0]/norm_xhat*(eye<2,2>() - outer_product( xhat/norm_xhat, xhat/norm_xhat));
    dn(0,0)-= dnhat(0,0);
    dn(0,1)-= dnhat(0,1);
    dn(1,0)-= dnhat(1,0);
    dn(1,1)-= dnhat(1,1);
    dn*= 1./l;

    const double c= std::cos( omega*t),
                 s= std::sin( omega*t);
    const Point3DCL z( toroidal_flow( p, -t));
    const Point3DCL dz0( MakePoint3D(  c, s, t*(-s*p[0] + c*p[1]))),
          dz1( MakePoint3D( -s, c, t*(-c*p[0] - s*p[1])));
    const Point3DCL grad_u= std::exp( -t)*(z[1]*dz0 + z[0]*dz1);
    SMatrixCL<3,3> Hess_u;
    Hess_u(0,2)= -c*z[0] - s*z[1];
    Hess_u(1,2)= -s*z[0] + c*z[1];
    Hess_u(2,0)= -c*z[0] - s*z[1];
    Hess_u(2,1)= -s*z[0] + c*z[1];
    Hess_u(2,2)= t*(z[0]*(s*p[0] - c*p[1]) + z[1]*(-c*p[0] - s*p[1]));
    Hess_u*= t;
    Hess_u+= outer_product( dz0, dz1) + outer_product( dz1, dz0);
    Hess_u*= std::exp( -t);

    const double u= toroidal_flow_sol( p, t),
                 mat_der=  -u,
                 reaction= div_gamma_wind( n, dw)*u,
                 diffusion= -laplace_beltrami_u( n, dn, grad_u, Hess_u);

    return mat_der + reaction + diffusion;
}
static RegisterScalarFunction regsca_toroidal_flow_rhs( "ToroidalFlowRhs", toroidal_flow_rhs);


// ==non-stationary test case "AxisScaling"==

double axis_scaling (double t)
{
    return 1. + 0.25*std::sin( t);
}

DROPS::Point3DCL axis_scaling_wind (const DROPS::Point3DCL& p, double t)
{
    return MakePoint3D( 0.25*std::cos( t)/axis_scaling( t)*p[0], 0., 0.);
}
static RegisterVectorFunction regvec_axis_scaling_wind( "AxisScalingWind", axis_scaling_wind);

double axis_scaling_lset (const Point3DCL& p, double t)
{
    return std::pow( p[0]/axis_scaling( t), 2) + std::pow( p[1], 2) + std::pow( p[2], 2) - 1.;
}
static RegisterScalarFunction regsca_axis_scaling_lset( "AxisScalingLset", axis_scaling_lset);

double axis_scaling_lset_ini (const Point3DCL& p, double)
{
    static const double t_end= P.get<double>( "Time.FinalTime");
    const double tout= t_end <= M_PI/2. ? t_end : M_PI/2.,
                 tin= t_end >= M_PI*3./2. ? M_PI*3./2. : (t_end >= M_PI ? t_end : 0.),
                 lout= axis_scaling_lset( p, tout),
                 lin= axis_scaling_lset( p, tin);
    return lout >= 0. ? lout :
           lin  <= 0. ? lin  :
           0.;
}

double axis_scaling_sol (const Point3DCL& p, double t)
{
    return std::exp( -0.5*t)*heat_conduction_u0( p);
}
static RegisterScalarFunction regsca_axis_scaling_sol( "AxisScalingSol", axis_scaling_sol);

double axis_scaling_rhs (const Point3DCL& p, double t)
{
    const double a= axis_scaling( t);
    const double bf4= (p/MakePoint3D( std::pow( a, 2), 1., 1.)).norm_sq();
    const double bf6= (p/MakePoint3D( std::pow( a, 3), 1., 1.)).norm_sq();

    const double mat_der=  0.25*std::cos( t)/a - 0.5;
    const double reaction= 0.25*std::cos( t)/a/bf4*( std::pow( p[1], 2) + std::pow( p[2], 2));
    const double diffusion= (-2./std::pow( a, 2) - (1. + 1./std::pow( a, 2))*( 2. + 1./std::pow( a, 2) - bf6/bf4))/bf4;

//     const Point3DCL tt( p/MakePoint3D( std::pow( a, 2), 1., 1.));
//     const double l= tt.norm();
//     const Point3DCL n( tt/l);
//     SMatrixCL<3,3> dn( (eye<3,3>() - outer_product( n, n))/l );
//     dn(0,0)/= std::pow( a, 2); dn(1,0)/= std::pow( a, 2); dn(2,0)/= std::pow( a, 2);
//
//     const Point3DCL w( MakePoint3D( 0.25*p[0]*std::cos( t)/a, 0., 0.));
//     SMatrixCL<3,3> dw;
//     dw(0,0)= 0.25*std::cos( t)/a;
//
//     const Point3DCL grad_u( std::exp( -0.5*t)*MakePoint3D( p[1], p[0], 0.));
//     SMatrixCL<3,3> Hess_u;
//     Hess_u(0,1)= 1.;
//     Hess_u(1,0)= 1.;
//     Hess_u*= std::exp( -0.5*t);
//
//     const double err= div_gamma_wind( n, dw) - reaction,
//                errlb= laplace_beltrami_u( n, dn, grad_u, Hess_u) - diffusion*axis_scaling_sol( p, t);
//     if (std::fabs( err) > 1e-12 || std::fabs( errlb) > 1e-12) {
//         std::cerr << err   << " " << div_gamma_wind( n, dw) << " " << reaction << "\n"
//                   << errlb << " " << laplace_beltrami_u( n, dn, grad_u, Hess_u) << " " << diffusion*axis_scaling_sol( p, t) << "\n";
//         exit( 1);
//     }
//     return (mat_der + div_gamma_wind( n, dw))*axis_scaling_sol( p, t) - laplace_beltrami_u( n, dn, grad_u, Hess_u);

    return (mat_der + reaction - diffusion)*axis_scaling_sol( p, t);
}
static RegisterScalarFunction regsca_axis_scaling_rhs( "AxisScalingRhs", axis_scaling_rhs);


// ==non-stationary test case "Collision"==

const double collision_p= 3.0;

// sphere 1: c1= (-1.5 0 0) + v*t
Point3DCL collision_center_1 (double t)
{
    return MakePoint3D( -1.5, 0., 0.) + WindVelocity*t;
}
// sphere 2: c2= ( 1.5 0 0) - v*t = -c1
Point3DCL collision_center_2 (double t)
{
    return -collision_center_1( t);
}

double collision_Dt_lset (const DROPS::Point3DCL& x, double t)
{
    Point3DCL x1= x - collision_center_1( t),
              x2= x - collision_center_2( t);
    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p + 2.),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p + 2.);
    return collision_p*(DROPS::inner_prod( WindVelocity, x1)/n1
                        - DROPS::inner_prod( WindVelocity, x2)/n2);
}

DROPS::Point3DCL collision_Dx_lset (const DROPS::Point3DCL& x, double t)
{
    Point3DCL x1= x - collision_center_1( t),
              x2= x - collision_center_2( t);
    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p + 2.),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p + 2.);
    return collision_p*(x1/n1 + x2/n2);
}

DROPS::Point3DCL collision_wind (const DROPS::Point3DCL& x, double t)
{
    const Point3DCL Dphi= collision_Dx_lset( x, t);
    const double Dtphi=   collision_Dt_lset( x, t);
    const double n= Dphi.norm_sq();
    return n < 1e-6 ? DROPS::Point3DCL() : (Dtphi/n)*Dphi;
}
static RegisterVectorFunction regvec_collision_wind( "CollisionWind", collision_wind);

double collision_lset (const Point3DCL& x, double t)
{
    const Point3DCL c1= collision_center_1( t);
    const Point3DCL x1= x - c1;
    const Point3DCL c2= collision_center_2( t);
    const Point3DCL x2= x - c2;

    const double n1= x1.norm() < 1e-3 ? 1e-3 : std::pow( x1.norm(), collision_p),
                 n2= x2.norm() < 1e-3 ? 1e-3 : std::pow( x2.norm(), collision_p);
    return RadDrop[0] - 1./n1  - 1./n2;
}
static RegisterScalarFunction regsca_collision_lset( "CollisionLset", collision_lset);

double collision_sol (const Point3DCL& x, double)
{
    return 2.*std::cos( x[0])*std::cos(M_PI*x[1]);
}
static RegisterScalarFunction regsca_collision_sol( "CollisionSol", collision_sol);

double collision_sol2 (const Point3DCL& x, double)
{
    return x[0] < 0. ? 0 : 3. - x[0];
}
static RegisterScalarFunction regsca_collision_sol2( "CollisionSol2", collision_sol2);

double collision_rhs (const Point3DCL&, double)
{
    return 0;
}
static RegisterScalarFunction regsca_collision_rhs( "CollisionRhs", collision_rhs);

double zero_function (const Point3DCL&, double)
{
    return 0;
}
static RegisterScalarFunction regsca_zero_sca_fun( "ZeroScalarFun", zero_function);

DROPS::Point3DCL mergedrop_wind (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL  c1,c2;
    c1[0]=1.5*(t-1),c1[1]=0,c1[2]=0;
    c2[0]=-1.5*(t-1),c2[1]=0,c2[2]=0;

    DROPS::Point3DCL x1( p - c1);
    DROPS::Point3DCL x2( p - c2);

    Point3DCL wd;
    /* if(x1.norm()!=0&&x2.norm()!=0)
     {
         double a=-4.5*x1[0]/std::pow(x1.norm(),5) + 4.5*x2[0]/std::pow(x2.norm(),5);
         wd=3/std::pow(x1.norm(),5)*x1+3/std::pow(x2.norm(),5)*x2;
         if(wd.norm_sq()!=0)
             wd=-a/wd.norm_sq()*wd;
         else
            wd=-a*wd;
     }*/
    double a=-4.5*x1[0]*std::pow(x2.norm(),5) + 4.5*x2[0]*std::pow(x1.norm(),5);
    wd=3*std::pow(x2.norm(),5)*x1+3*std::pow(x1.norm(),5)*x2;
    if(wd.norm_sq()!=0)
        wd=-a/wd.norm_sq()*wd;
    else
        wd=-a*wd;
    return wd;
}
static RegisterVectorFunction regvec_mergedrop_wind( "MergeDropWind", mergedrop_wind);

double mergedrop_lset (const Point3DCL& p, double t)
{
    // std::cout<<t<<" ";
    DROPS::Point3DCL  c1,c2;
    c1[0]=1.5*(t-1),c1[1]=0,c1[2]=0;
    c2[0]=-1.5*(t-1),c2[1]=0,c2[2]=0;

    DROPS::Point3DCL x1( p - c1);
    DROPS::Point3DCL x2( p - c2);
    // std::cout<<1.0-1./std::pow(x1.norm(),3)<<std::endl;
    if(x1.norm()!=0&&x2.norm()!=0)
        return  1- 1./std::pow(x1.norm(),3)-1.0/std::pow(x2.norm(),3);
    else
        return   -10*(std::pow(x1.norm(),3)*std::pow(x2.norm(),3)-std::pow(x1.norm(),3)-std::pow(x2.norm(),3));//-1./std::pow(x2.norm(),3);

}
static RegisterScalarFunction regsca_mergedrop_lset( "MergeDropLset", mergedrop_lset);

DROPS::Point3DCL normal_mergedrop (const Point3DCL& p, double t)
{
    // std::cout<<t<<" ";
    DROPS::Point3DCL  c1,c2;
    c1[0]=1.5*(t-1),c1[1]=0,c1[2]=0;
    c2[0]=-1.5*(t-1),c2[1]=0,c2[2]=0;

    DROPS::Point3DCL x1( p - c1);
    DROPS::Point3DCL x2( p - c2);

    DROPS::Point3DCL grad;
    if(x1.norm()!=0&&x2.norm()!=0)
    {
        grad=-3*(std::pow(x1.norm(),5)*x1+std::pow(x2.norm(),5)*x2);
        return  grad.norm()>0.000001?grad/grad.norm():grad;
    }
    else
        return   0*grad;//-1./std::pow(x2.norm(),3);

}
static RegisterVectorFunction regsca_normal_mergedrop( "NormalMergeDrop", normal_mergedrop);


double mergedrop_sol (const Point3DCL& p, double t)
{
    if(p[0]>=0)
        return   3-p[0];
    else
        return 0;
}
static RegisterScalarFunction regsca_mergedrop_sol( "MergeDropSol", mergedrop_sol);

DROPS::Point3DCL mergedrop_surfgradsol (const DROPS::Point3DCL& p, double t)
{
    Point3DCL surfgrad;
    return surfgrad;
}
static RegisterVectorFunction regvec_mergedrop_surfgradsol( "MergeDropSurfGradSol", mergedrop_surfgradsol);


double mergedrop_rhs (const Point3DCL& p, double t)
{
    return 0; //Ex3
}
static RegisterScalarFunction regsca_mergedrop_rhs( "MergeDropRhs", mergedrop_rhs);
//////////////////////////////////////////////////////////////////////////////////////////////
// ==Some spheres==

double sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    //double value=0;
    //lsFun(x[0], x[1], x[2], &value);
    return x.norm() - RadDrop[0];
    //return value;
}
static RegisterScalarFunction regsca_sphere_dist_lset( "SphereDist", sphere_dist);


//static RegisterScalarFunction regsca_sphere_dist_lset( "SphereDist", sphere_dist);


DROPS::Point3DCL d_sphere_dist (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - PosDrop);
    return x/x.norm();
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (PosDrop + t*constant_wind(p, t)));
    return x.norm() - RadDrop[0];
}
static RegisterScalarFunction regsca_movell_lset( "MovingEllipsoid", sphere_2move);


DROPS::Point3DCL normal_sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x= (p - PosDrop - t*constant_wind(p, t));//--constant wind
    if(x.norm()>0.000001)
        return x/x.norm();
    else
        return x;
}
static RegisterVectorFunction regvec_moving_sphere_normal( "NormalMovingEllipsoid", normal_sphere_2move);

SMatrixCL<3,3> dp_sphere (const DROPS::Point3DCL& x, double)
{
    const double normx= x.norm();
    return normx == 0. ? SMatrixCL<3,3>() : RadDrop[0]/normx*(eye<3,3>() - outer_product( x/normx, x/normx));
}

double abs_det_sphere (const TetraCL& tet, const BaryCoordCL& xb, const SurfacePatchCL& p)
{
    if (p.empty())
        return 0.;

    // Compute the jacobian of p.
    SMatrixCL<3,3> dp= dp_sphere( GetWorldCoord( tet, xb), 0.);

    const Bary2WorldCoordCL b2w( tet);
    QRDecompCL<3,2> qr;
    SMatrixCL<3,2>& M= qr.GetMatrix();
    M.col(0, b2w( p.vertex_begin()[1]) - b2w( p.vertex_begin()[0]));
    M.col(1, b2w( p.vertex_begin()[2]) - b2w( p.vertex_begin()[0]));
    qr.prepare_solve();
    SMatrixCL<3,2> U;
    Point3DCL tmp;
    for (Uint i= 0; i < 2; ++i)
    {
        tmp= std_basis<3>( i + 1);
        qr.apply_Q( tmp);
        U.col( i, tmp);
    }
    const SMatrixCL<2,2> Gram= GramMatrix( dp*U);
    return std::sqrt( Gram(0,0)*Gram(1,1) - Gram(0,1)*Gram(1,0));
}

// ==stationary test case "LaplaceBeltramixyz", created by ls==
// Sphere around 0, RadDrop 1, wind == 0
// "Levelset": "SphereDist"
// f = 2*(x+y+z) u=x+y+z for problem -\Delta u + u = f
//double xyz_rhs (const DROPS::Point3DCL& p, double)
//{//my test case，f = 3*(x+y+z) for problem -\Delta u + u = f
//    return 3*(p[0]+p[1]+p[2])/p.norm();
//}




// ==stationary test case "LaplaceBeltrami0"==
// Sphere around 0, RadDrop 1, wind == 0
// "Levelset": "SphereDist"
// A right hand side from C.J. Heine...
// const double a( -13./8.*std::sqrt( 35./M_PI));
const double a( 12.);
double laplace_beltrami_0_rhs (const DROPS::Point3DCL& p, double)
{
    return a/std::pow( p.norm(), 3.)*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_rhs( "LaplaceBeltrami0Rhs", laplace_beltrami_0_rhs);

// ...and the corresponding solution (extended)
double laplace_beltrami_0_sol (const DROPS::Point3DCL& p, double)
{
    return p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
//     return 1./12.*laplace_beltrami_0_rhs( p, 0.);
//    return 1. + p.norm_sq()/(12. + p.norm_sq())*laplace_beltrami_0_rhs( p, 0.);
}
static RegisterScalarFunction regsca_laplace_beltrami_0_sol( "LaplaceBeltrami0Sol", laplace_beltrami_0_sol);

// The tangential gradient of laplace_beltrami_0_sol with respect to the exact sphere.
DROPS::Point3DCL laplace_beltrami_0_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)*( MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) - (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
    return tmp; // This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
}

double sol0t (const DROPS::Point3DCL& p, double t)
{
    const Point3DCL q( p - (PosDrop + t*constant_wind(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

//    return q.norm_sq()/(12. + q.norm_sq())*val;
    return 1. + q.norm_sq()/(12. + q.norm_sq())*val;
}

// ==stationary test case "LaplaceBeltrami1"==
// Torus with R= RadTorus[0]= 1., r= RadTorus[1]= 0.6, wind == 0
// "Levelset": "Torus"

// angle from positive x-axis to (x,y,0)
double t_angle (const Point3DCL& p, double)
{
    return std::atan2( p[1], p[0]);
}

// distance from the circle in the x-y-plane around 0 with radius R
double rho (const Point3DCL& p, double)
{
    return std::sqrt( p[2]*p[2] + std::pow( std::sqrt( p[0]*p[0] + p[1]*p[1]) - RadTorus[0], 2));
}

// angle from positive (x,y,0)-direction to p.
double p_angle (const Point3DCL& p, double)
{
    return std::atan2( p[2], std::sqrt( p[0]*p[0] + p[1]*p[1]) - RadTorus[0]);
}

double laplace_beltrami_1_rhs (const Point3DCL& p, double)
{
    const double pa= p_angle( p, 0.);
    const double ta= t_angle( p, 0.);

    using std::sin;
    using std::cos;
    using std::pow;
    const double t0= (9.*sin( 3.*ta)*cos( 3.*pa + ta))/(RadTorus[1]*RadTorus[1]);
    const double t1= -(-10.*sin( 3.*ta)*cos(3.*pa + ta) - 6.*cos( 3.*ta)*sin( 3.*pa + ta))
                     /pow(RadTorus[0] + RadTorus[1]*cos( pa), 2);
    const double t2= -(3.*sin( pa)*sin( 3.*ta)*sin( 3.*pa + ta))/(RadTorus[1]*(RadTorus[0] + RadTorus[1]*cos( pa)));
    return t0 + t1 + t2;
}
static RegisterScalarFunction regsca_laplace_beltrami_1_rhs( "LaplaceBeltrami1Rhs", laplace_beltrami_1_rhs);

double laplace_beltrami_1_sol (const Point3DCL& p, double)
{
    return std::sin(3.*t_angle( p, 0.))*std::cos( 3.*p_angle( p, 0.) + t_angle( p, 0.));
}
static RegisterScalarFunction regsca_laplace_beltrami_1_sol( "LaplaceBeltrami1Sol", laplace_beltrami_1_sol);

template<class DiscP1FunType>
double L2_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);//store reference tetra
    const double t= discsol.GetTime();//discsol is piecewise linear function
    QuadDomain2DCL qdom;//2D domain quadrature class
    std::valarray<double> qsol,//q is value on qdom
        qdiscsol;

    double d( 0.);
    DROPS_FOR_TRIANG_CONST_TETRA( discsol.GetMG(), ls.GetLevel(), it)//go through tetra
    {
        //   LocalP1CL<P1EvalCL<double, DROPS::NoBndDataCL<double>,
        //           DROPS::VecDescBaseCL<DROPS::VectorBaseCL<double> > >> local_p1_f(*it,discsol);
        //typename PEvalT::LocalFET( tet, f)
//        const BaryCoordCL& test {1,2,3,4};
//        LocalP1CL<double> local_p1_f(*it,discsol);
//        std::cout<<"here:"<<local_p1_f(test)<<std::endl;
        make_CompositeQuad5Domain2D (qdom, *it, lat, ls, lsbnd);// get quadrature domain
        resize_and_evaluate_on_vertexes ( discsol, *it, qdom, qdiscsol);//discsol.getDoF, get calculated value on qdom
        resize_and_evaluate_on_vertexes ( extsol,  *it, qdom, t, qsol);//get true value on qdom
        d+= quad_2D( std::pow( qdiscsol - qsol, 2), qdom);//lenght of q, 63 35 28
        //  std::cout<<d<<std::endl;
    }
    return std::sqrt( d);
}
//#ifdef Debug

//const DROPS::P1EvalCL<double, DROPS::NoBndDataCL<double>,
//     DROPS::VecDescBaseCL<DROPS::VectorBaseCL<double> > > *discsol_ref;
LocalP1CL<double> local_p1_f;
const P1EvalCL<double, NoBndDataCL<double>,
      VecDescBaseCL<VectorBaseCL<double> > > *discsol_ref;
LocalP1CL<double> *local_p1_f_ptr;
void errL2IntFun(double x, double y, double z, double *ff)
{
    double BaryCoordArr[4];
    for(int i=0; i<4; i++)
    {
        BaryCoordArr[i] = getBaryCoord(tet,i, x,y,z);
    }
    const BaryCoordCL& BaryCoord{BaryCoordArr[0],BaryCoordArr[1],BaryCoordArr[2],BaryCoordArr[3]};
    double qdiscsol = (*local_p1_f_ptr)(BaryCoord);
    const DROPS::Point3DCL& p{x,y,z};
    double qsol = laplace_beltrami_xyz_sol (p, 0.);
    *ff = std::pow( qdiscsol - qsol, 2);
}

template<class DiscP1FunType>
double L2_error_high_quad (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                           const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    const double t= discsol.GetTime();
    double d( 0.);
    double res(0.);
    DROPS_FOR_TRIANG_CONST_TETRA( discsol.GetMG(), ls.GetLevel(), it)
    {
        for (Uint i= 0; i < 4; ++i)
        {
            auto vtx = (*it).GetVertex(i);
            auto coord = vtx->GetCoord();
            for(int j=0; j<4; j++)
            {
                tet[i][j] = coord[j];
                //  std::cout<<tet[i][j]<<std::endl;
            }
        }
        LocalP1CL<double> local_p1_f(*it,discsol);
        local_p1_f_ptr = &local_p1_f;
        int n = phgQuadInterface2(
                    lsFun,		/* the level set function */
                    2,		/* polynomial order of the level set function */
                    lsGrad,	/* the gradient of the level set function */
                    tet,		/* coordinates of the vertices of the tetra */
                    errL2IntFun,		/* the integrand */
                    1,		/* dimension of the integrand */
                    DOF_PROJ_NONE,	/* projection type for surface integral */
                    0,		/* integration type (-1, 0, 1) */
                    orderG,		/* order of the 1D Gaussian quadrature */
                    &res,		/* the computed integral */
                    NULL		/* pointer returning the computed rule */
                );
        d += res;
    }
    return std::sqrt( d);
}
//#endif

/// The nodal interpolant of extsol on the interface-FE-space ist computed first.
/// The H1-error is then computed between the interpolant and the numerical solution.
template<class DiscP1FunType>
double H1_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    IdxDescCL* idx= const_cast<IdxDescCL*>( discsol.GetSolution()->RowIdx);
    MatDescCL A( idx, idx);
    SetupLBP1( discsol.GetMG(), &A, ls, lsbnd, 1.);
    VecDescCL sol_vec( idx);
    P1Init (extsol, sol_vec, discsol.GetMG(), discsol.GetTime());
    // sol_vec.t= discsol.GetTime();
    // SetupInterfaceRhsP1( discsol.GetMG(), &sol_vec, ls, lsbnd, extsol);
    const VectorCL diff( discsol.GetSolution()->Data - sol_vec.Data);
    return std::sqrt( dot(A.Data*diff, diff));
}



template<class DiscP1FunType>
double H1_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol,DROPS::IdxDescCL &ifaceidx)
{
    //IdxDescCL* idx= const_cast<IdxDescCL*>( discsol.GetSolution()->RowIdx);
    IdxDescCL* idx = &ifaceidx;
    MatDescCL A( idx, idx);
    SetupLBP1( discsol.GetMG(), &A, ls, lsbnd, 1.);
    VecDescCL sol_vec( idx);
    P1Init (extsol, sol_vec, discsol.GetMG(), discsol.GetTime());
    // sol_vec.t= discsol.GetTime();
    // SetupInterfaceRhsP1( discsol.GetMG(), &sol_vec, ls, lsbnd, extsol);
    const VectorCL diff( discsol.GetSolution()->Data - sol_vec.Data);
    return std::sqrt( dot(A.Data*diff, diff));
}


template<class DiscP1FunType>
double H1_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                 const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol,DROPS::IdxDescCL &ifaceidx, MatDescCL &A)
{
    //IdxDescCL* idx= const_cast<IdxDescCL*>( discsol.GetSolution()->RowIdx);
    IdxDescCL* idx = &ifaceidx;
    //MatDescCL A( idx, idx);
    //SetupLBP1( discsol.GetMG(), &A, ls, lsbnd, 1.);
    VecDescCL sol_vec( idx);
    P1Init (extsol, sol_vec, discsol.GetMG(), discsol.GetTime());
    // sol_vec.t= discsol.GetTime();
    // SetupInterfaceRhsP1( discsol.GetMG(), &sol_vec, ls, lsbnd, extsol);
    const VectorCL diff( discsol.GetSolution()->Data - sol_vec.Data);
    return std::sqrt( dot(A.Data*diff, diff));
}

template<class DiscP1FunType>
double H1_error_p2 (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                    const DiscP1FunType& discsol, VecDescCL &sol_vec,DROPS::IdxDescCL &ifaceidx, MatDescCL &A)
{
    //IdxDescCL* idx= const_cast<IdxDescCL*>( discsol.GetSolution()->RowIdx);
    //IdxDescCL* idx = &ifaceidx;
    //MatDescCL A( idx, idx);
    //SetupLBP1( discsol.GetMG(), &A, ls, lsbnd, 1.);
    //VecDescCL sol_vec( idx);
    //P2Init (extsol, sol_vec, discsol.GetMG(), discsol.GetTime());
    // sol_vec.t= discsol.GetTime();
    // SetupInterfaceRhsP1( discsol.GetMG(), &sol_vec, ls, lsbnd, extsol);
    const VectorCL diff( discsol.GetSolution()->Data - sol_vec.Data);
    return std::sqrt( dot(A.Data*diff, diff));
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
                DROPS::instat_scalar_fun_ptr extsol)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);
    const double t= ls.t;
    QuadDomain2DCL qdom;
    std::valarray<double> qsol;

    double d( 0.);
    DROPS_FOR_TRIANG_CONST_TETRA( mg, ls.GetLevel(), it)
    {
        make_CompositeQuad5Domain2D (qdom, *it, lat, ls, lsbnd);
        resize_and_evaluate_on_vertexes ( extsol,  *it, qdom, t, qsol);
        d+= quad_2D( qsol*qsol, qdom);
    }
    return std::sqrt( d);
}

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::instat_scalar_fun_ptr d, double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
    ls.Data[it->Unknowns( idx)]= ls.Data[it->Unknowns( idx)]=
                                     0.5*(ls.Data[it->GetVertex( 0)->Unknowns( idx)] + ls.Data[it->GetVertex( 1)->Unknowns( idx)]);
    ls.t= t;
}

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
    ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
    ls.t= t;
}

void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, instat_vector_fun_ptr LsgVel, double t= 0.)
{
    VectorCL& lsgvel= vec->Data;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, sit)
    {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                                                   LsgVel(sit->GetCoord(), t));
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, sit)
    {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                                                   LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t));
    }
    vec->t= t;
}


SurfactantP1BaseCL* make_surfactant_timedisc( MultiGridCL& mg, LevelsetP2CL& lset,
        VecDescCL& v, const BndDataCL<Point3DCL>& Bnd_v,
        const ParamCL& P, const double & dist=0)
{
    SurfactantP1BaseCL* ret= 0;
    const std::string method= P.get<std::string>( "SurfTransp.Method");

    if (method == std::string( "cGcG"))
        ret= new SurfactantcGP1CL( mg,
                                   P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
                                   &v, Bnd_v, lset.Phi, lset.GetBndData(),
                                   P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
                                   P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "spacetime-cGdG"))
        ret= new SurfactantSTP1CL( mg,
                                   P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
                                   &v, Bnd_v, lset.Phi, lset.GetBndData(),
                                   /* cG_in_t_ */ false, /* use_mass_div */ P.get<bool>( "SurfTransp.UseMassDiv"),
                                   P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
                                   P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "spacetime-cGcG"))
        ret= new SurfactantSTP1CL( mg,
                                   P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
                                   &v, Bnd_v, lset.Phi, lset.GetBndData(),
                                   /* cG_in_t_ */ true, /* use_mass_div */ P.get<bool>( "SurfTransp.UseMassDiv"),
                                   P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
                                   P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "characteristic-transport"))
        ret= new SurfactantCharTransportP1CL ( mg,
                                               P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
                                               &v, Bnd_v, lset.Phi, lset.GetBndData(),
                                               P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"),
                                               P.get<double>("SurfTransp.XFEMReduced"));
    else if (method == std::string( "NarrowBandStabilization"))
    {
        std::cout<<"Test Method NarrowBand :"<<dist<<std::endl;
        ret= new SurfactantNarrowBandStblP1CL ( mg,
                                                P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"),
                                                &v, Bnd_v, lset.Phi, lset.GetBndData(),the_normal_fun, dist,
                                                1./(P.get<DROPS::Point3DCL>("Mesh.E1")[0]/(P.get<double>("Mesh.N1"))),
                                                P.get<int>("SurfTransp.Solver.Iter"),P.get<double>("SurfTransp.Solver.Tol"),
                                                P.get<double>("SurfTransp.XFEMReduced"));/**/
    }
    /* else if (method == std::string( "CahnHilliardNarrowBand"))
     {
         std::cout<<"Test Method NarrowBand :"<<dist<<std::endl;
         ret= new CahnHilliardNarrowBandStblP1CL ( mg,
                                                 P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"), 0.1,
                                                 &v, Bnd_v, lset.Phi, lset.GetBndData(),the_normal_fun, dist,
                                                 1./(P.get<DROPS::Point3DCL>("Mesh.E1")[0]/(P.get<double>("Mesh.N1"))), 1,
                                                 P.get<int>("SurfTransp.Solver.Iter"),P.get<double>("SurfTransp.Solver.Tol")
                                                 );

     }/**/
    else
        throw DROPSErrCL( std::string( "make_surfactant_timedisc: Unknown method '") + method + std::string( "'.\n"));

    return ret;
}


void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;

    if (P.get<std::string>("SurfTransp.Exp.Levelset") == std::string( "AxisScalingLset"))
        dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( axis_scaling_lset_ini);

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun, 0.);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    //DROPS::LevelsetP2CL& lset2( *LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")) );
    //lset2.idx.CreateNumbering( mg.GetLastLevel(), mg);
    //lset2.Phi.SetIdx( &lset2.idx);
    //LSInit( mg, lset2.Phi, the_lset_fun, 0.);

    const double Vol= lset.GetVolume();
    lset.InitVolume( Vol);
    std::cout << "droplet volume: " << Vol << std::endl;

    BndDataCL<Point3DCL> Bnd_v( 6, bc_wind, bf_wind);
    IdxDescCL vidx( vecP2_FE);
    vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( mg, &v, Bnd_v, the_wind_fun, 0.);

    //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps"));
    double dist=2*P.get<DROPS::Point3DCL>("SurfTransp.Exp.Velocity").norm()*P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps")
                +2*P.get<DROPS::Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")/pow(2,P.get<int>("Mesh.AdaptRef.FinestLevel")+1);

    std::unique_ptr<SurfactantP1BaseCL> timediscp( make_surfactant_timedisc( mg, lset, v, Bnd_v, P,dist));
    SurfactantP1BaseCL& timedisc= *timediscp;
    timedisc.SetRhs( the_rhs_fun);

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    InterfaceP1RepairCL ic_repair( mg, lset.Phi, lset.GetBndData(), timedisc.ic);
    adap.push_back( &ic_repair);
    //LevelsetRepairCL lset2repair( lset2);
    //adap.push_back( &lset2repair);

    // Init Interface-Sol
    timedisc.idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData(), dist);
    std::cout << "NumUnknowns: " << timedisc.idx.NumUnknowns() << std::endl;
    timedisc.ic.SetIdx( &timedisc.idx);
    timedisc.SetInitialValue( the_sol_fun, 0.);
    timedisc.iface.SetIdx( &timedisc.idx);
    timedisc.iface_old.SetIdx( &timedisc.idx);

    BndDataCL<> nobnd( 0);
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    if (vtkwriter.get() != 0)
    {
        vtkwriter->Register( make_VTKScalar(      lset.GetSolution(),              "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.ic,                 "InterfaceSol"));
        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.iface,                 "interface_mesh"));
        vtkwriter->Register( make_VTKIfaceScalar( mg, timedisc.iface_old,                 "old_interface_mesh"));
        vtkwriter->Register( make_VTKVector(      make_P2Eval( mg, Bnd_v, v),      "Velocity"));
        //vtkwriter->Register( make_VTKScalar(      lset2.GetSolution(),             "Levelset2"));
        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Write( 0.);
    }
    //if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
    //    DROPS::WriteFEToFile( timedisc.ic, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path"), P.get<bool>( "SolutionOutput.Binary"));

    const double dt= P.get<double>("Time.FinalTime")/P.get<double>("Time.NumSteps");
    double L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "L_2x-error: " << L_2x_err
              << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
              << std::endl;
    double L_inftL_2x_err= L_2x_err;
    std::cout << "L_inftL_2x-error: " <<  L_inftL_2x_err << std::endl;
    double H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
    std::cout << "H_1x-error: " << H_1x_err << std::endl;
    double L_2tH_1x_err_sq= 0.5*dt*std::pow( H_1x_err, 2);
    BndDataCL<> ifbnd( 0);
    std::cout << "initial surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';

    dynamic_cast<DistMarkingStrategyCL*>( adap.get_marking_strategy())->SetDistFct( lset);

    double total_mass=Integral_Gamma(mg, lset.Phi,lset.GetBndData(),make_P1Eval(mg, ifbnd, timedisc.ic));
    //In the first step, we use some smaller time step to solve the problem
    // double discrete_mass=0;
    int N1=(int)1./dt;
    for (int step= 1; step <= N1; ++step)
    {
        const double cur_time= step*dt/N1;
        timedisc.InitTimeStep();

        LSInit( mg, lset.Phi, the_lset_fun, cur_time);
        InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
        timedisc.DoStep0( cur_time);
        total_mass=Integral_Gamma(mg, lset.Phi,lset.GetBndData(),make_P1Eval(mg, ifbnd, timedisc.ic));
        //  if(step==1)
        //   	discrete_mass=total_mass;
        //L_2tL_2x_err_sq+= (step > 1 ? 0.5 : 0.)*dt/10*std::pow( L_2x_err, 2);
        L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "L_2x-error: " << L_2x_err
                  << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
                  << std::endl;
        L_inftL_2x_err= std::max( L_inftL_2x_err, L_2x_err);
        std::cout << "L_inftL_2x-eerror: " << L_inftL_2x_err << std::endl;
        //L_2tL_2x_err_sq+= 0.5*dt/10*std::pow( L_2x_err, 2);
        // std::cout << "L_2tL_2x-error: " << std::sqrt( L_2tL_2x_err_sq) << std::endl;
        L_2tH_1x_err_sq+= (step > 1 ? 0.5 : 0.)*dt/N1*std::pow( H_1x_err, 2);
        //    H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "H_1x-error: " << H_1x_err << std::endl;
        L_2tH_1x_err_sq+= 0.5*dt/N1*std::pow( H_1x_err, 2);
        //    L_2tH_1x_err_sq+= dt*std::pow( H_1x_err, 2);
        std::cout << "L_2tH_1x-error: " << std::sqrt( L_2tH_1x_err_sq) << std::endl;
    }

    for (int step= 2; step <= P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "======================================================== step " << step << ":\n";
        ScopeTimerCL timer( "Strategy: Time-loop");
        const double cur_time= step*dt;
        // Assumes (as the rest of Drops), that the current triangulation is acceptable to perform the time-step.
        // If dt is large and AdapRef.Width is small, this may not be true.
        // Watch for large differences in numbers of old and new dof.
        timedisc.InitTimeStep();
        LSInit( mg, lset.Phi, the_lset_fun, cur_time);
        InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
        timedisc.DoStep( cur_time);
        std::cout << "surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';
        L_2x_err= L2_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "L_2x-error: " << L_2x_err
                  << "\nnorm of true solution: " << L2_norm( mg, lset.Phi, lset.GetBndData(), the_sol_fun)
                  << std::endl;
        L_inftL_2x_err= std::max( L_inftL_2x_err, L_2x_err);
        std::cout << "L_inftL_2x-error: " << L_inftL_2x_err << std::endl;
        L_2tH_1x_err_sq+= (step > 1 ? 0.5 : 0.)*dt*std::pow( H_1x_err, 2);
        H_1x_err= H1_error( lset.Phi, lset.GetBndData(), timedisc.GetSolution(), the_sol_fun);
        std::cout << "H_1x-error: " << H_1x_err << std::endl;
        L_2tH_1x_err_sq+= 0.5*dt*std::pow( H_1x_err, 2);
        std::cout << "L_2tH_1x-error: " << std::sqrt( L_2tH_1x_err_sq) << std::endl;
        if (vtkwriter.get() != 0 && step % P.get<int>( "VTK.Freq") == 0)
        {
            LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ cur_time);
            vtkwriter->Write( cur_time);
        }
        if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0 && step % P.get<int>( "SurfTransp.SolutionOutput.Freq") == 0)
        {
            std::ostringstream os1,
                os2;
            os1 << P.get<int>( "Time.NumSteps");
            os2 << P.get<std::string>( "SurfTransp.SolutionOutput.Path") << std::setw( os1.str().size()) << step;
            DROPS::WriteFEToFile( timedisc.ic, mg, os2.str(), P.get<bool>( "SurfTransp.SolutionOutput.Binary"));
        }
//        lset2.DoStep();
//        VectorCL rhs( lset2.Phi.Data.size());
//        lset2.ComputeRhs( rhs);
//        lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, cur_time));
//        lset2.SetTimeStep( dt);
//        lset2.DoStep( rhs);

//         std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
//         if (P.get("Levelset.VolCorr", 0)) {
//             double dphi= lset.AdjustVolume( Vol, 1e-9);
//             std::cout << "volume correction is " << dphi << std::endl;
//             lset.Phi.Data+= dphi;
//             std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
//         }
        //if (C.rpm_Freq && step%C.rpm_Freq==0) { // reparam levelset function
        // lset.ReparamFastMarching( C.rpm_Method);
        const bool doGridMod= P.get<int>("Mesh.AdaptRef.Freq") && step%P.get<int>("Mesh.AdaptRef.Freq") == 0;
        const bool gridChanged= doGridMod ? adap.UpdateTriang() : false;
        if (gridChanged)
        {
            std::cout << "Triangulation changed.\n";
            vidx.DeleteNumbering( mg);
            vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
            v.SetIdx( &vidx);
            InitVel( mg, &v, Bnd_v, the_wind_fun, cur_time);
            LSInit( mg, lset.Phi, the_lset_fun, cur_time);
            the_sol_vd.SetIdx( &lset.idx);
            LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ cur_time);
            // timedisc.Update(); // Called unconditionally in DoStep.

            //lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), dt);

            std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            lset.AdjustVolume();
            lset.GetVolumeAdjuster()->DebugOutput( std::cout);
        }
    }
    std::cout << std::endl;
    //delete &lset2;
}

#define HighQuadP1
#define P1_SURFAPPROX
void StationaryStrategyP1 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    adap.MakeInitialTriang();

    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    DROPS::IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifaceidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    //only consider the dofs on tet which interact with surface

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
#ifndef HighQuadP1
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
#else
    DROPS::SetupInterfaceMassP1HighQuad( mg, &M, lset.Phi, lset.GetBndData());
#endif
    std::cout << "M is set up.\n";


    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
#ifndef HighQuadP1
    DROPS::SetupLBP1( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfTransp.Visc"));
#else
    DROPS::SetupLBP1HighQuad( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfTransp.Visc"));
#endif
    std::cout << "A is set up.\n";

    DROPS::MatrixCL L;
    L.LinComb( 1.0, A.Data, 1.0, M.Data);
//   DROPS::MatrixCL& L= A.Data;
    DROPS::VecDescCL b( &ifaceidx);
#ifndef HighQuadP1
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, lset.GetBndData(), the_rhs_fun);
#else
    DROPS::SetupInterfaceRhsP1HighQuad(mg, &b, lset.Phi, lset.GetBndData(), the_rhs_fun);
#endif
    DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    DROPS::WriteFEToFile( b, mg, "rhs_iface.txt", /*binary=*/ false);


#if 1
    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data, x.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
#endif

#if 0
//define direct solver
    DROPS::VecDescCL x( &ifaceidx);
    DROPS::DirectSymmSolverCL dsolver(L);
    dsolver.Solve(L,x.Data,b.Data);
    //dsolver.Update(A);
    //dsolver.Solve(A,x,b);
    //DROPS::DirectNonSymmSolverCL dnsolver(A);
    //dnsolver.Solve(A,x,b);

#endif

    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( x, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path"), P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);//only consider the dof interact with exact surface
    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);//fixup the edge-dofs, much more, how to understand???
    DROPS::NoBndDataCL<> nobnd;
    if (vtkwriter.get() != 0)
    {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, x, "InterfaceSol"));
        vtkwriter->Write( 0.);
    }
#ifndef HighQuadP1
    double L2_err( L2_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), the_sol_fun));
    //double H1_err= H1_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, x), the_sol_fun,ifaceidx);
#else
    double L2_err( L2_error_high_quad( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), the_sol_fun));
#endif
    double H1_err= H1_error( lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, x), the_sol_fun,ifaceidx,A);
    std::cout << "L_2-error: " << L2_err<<std::endl;
    std::cout << "H_1-error: "<< H1_err <<std::endl;
}

/// \brief Accumulate L2-norms and errors on the higher order zero level.
/// Works for P1IF_FE, P2IF_FE, and C-functions. All functions are evaluated on the P2-levelset.
class InterfaceL2AccuP2CL : public TetraAccumulatorCL
{
private:
    const InterfaceCommonDataP2CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    NoBndDataCL<> nobnddata;
    const VecDescCL* fvd;

    instat_scalar_fun_ptr f;
    double f_time;

    LocalLaplaceBeltramiP2CL loc_lb;

    instat_vector_fun_ptr f_grad;
    double f_grad_time;

    InterfaceL2AccuP2CL* tid0p; // The object in OpenMP-thread 0, in which the following variables are updated.
    std::vector<double> f_grid_norm,
        f_grid_int,
        f_norm,
        f_int,
        err,
        area,
        f_grid_grad_norm,
        f_grad_norm,
        grad_err;

public:
    double f_grid_norm_acc,
           f_grid_int_acc,
           f_norm_acc,
           f_int_acc,
           err_acc,
           area_acc,
           f_grid_grad_norm_acc,
           f_grad_norm_acc,
           grad_err_acc;

    InterfaceL2AccuP2CL (const InterfaceCommonDataP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), f_grad( 0), f_grad_time( 0.) {}
    virtual ~InterfaceL2AccuP2CL () {}

    void set_name (const std::string& n)
    {
        name_= n;
    }
    void set_grid_function (const VecDescCL& fvdarg)
    {
        fvd= &fvdarg;
    }
    void set_function (const instat_scalar_fun_ptr farg, double f_time_arg= 0.)
    {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const instat_vector_fun_ptr farg, double f_time_arg= 0.)
    {
        f_grad= farg;
        f_grad_time= f_time_arg;
    }

    virtual void begin_accumulation ()
    {
        std::cout << "InterfaceL2AccuP2CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\".\n";

        tid0p= this;
        f_grid_norm.clear();
        f_grid_norm.resize( omp_get_max_threads(), 0.);
        f_grid_int.clear();
        f_grid_int.resize( omp_get_max_threads(), 0.);

        f_norm.clear();
        f_norm.resize( omp_get_max_threads(), 0.);
        f_int.clear();
        f_int.resize( omp_get_max_threads(), 0.);

        err.clear();
        err.resize( omp_get_max_threads(), 0.);
        area.clear();
        area.resize( omp_get_max_threads(), 0.);

        f_grid_grad_norm.clear();
        f_grid_grad_norm.resize( omp_get_max_threads(), 0.);
        f_grad_norm.clear();
        f_grad_norm.resize( omp_get_max_threads(), 0.);
        grad_err.clear();
        grad_err.resize( omp_get_max_threads(), 0.);

        f_grid_norm_acc= f_grid_int_acc= f_norm_acc= f_int_acc= err_acc= area_acc= 0.;
        f_grid_grad_norm_acc= f_grad_norm_acc= grad_err_acc= 0.;
    }

    virtual void finalize_accumulation()
    {
        std::cout << "InterfaceL2AccuP2CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";//"P2-solution"
        area_acc= std::accumulate( area.begin(), area.end(), 0.);
        std::cout << "\n\tarea: " << area_acc;//output area
        if (fvd != 0)//grid solution: norm and integral
        {
            f_grid_norm_acc=  std::sqrt( std::accumulate( f_grid_norm.begin(), f_grid_norm.end(), 0.));
            f_grid_int_acc= std::accumulate( f_grid_int.begin(), f_grid_int.end(), 0.);
            std::cout << "\n\t|| f_grid ||_L2: " << f_grid_norm_acc
                      << "\tintegral: " << f_grid_int_acc;
        }
        if (f != 0)//eaxct solution
        {
            f_norm_acc=  std::sqrt( std::accumulate( f_norm.begin(), f_norm.end(), 0.));
            f_int_acc= std::accumulate( f_int.begin(), f_int.end(), 0.);
            std::cout << "\n\t|| f ||_L2: " << f_norm_acc
                      << "\t integral: " << f_int_acc;
        }
        if (fvd != 0 && f != 0)//both
        {
            err_acc=  std::sqrt( std::accumulate( err.begin(), err.end(), 0.));//accumulate err
            std::cout << "\n\t|| f - f_grid ||_L2: " << err_acc;

            const double mvf_err= std::sqrt( std::pow( err_acc, 2) - std::pow( f_grid_int_acc - f_int_acc, 2)/area_acc);
            std:: cout << "\t|| f - c_f - (f_grid -c_{f_grid}) ||_L2: " << mvf_err;//some kind of "mean"
        }

        if (fvd != 0)//grid solution for gradient
        {
            f_grid_grad_norm_acc=  std::sqrt( std::accumulate( f_grid_grad_norm.begin(), f_grid_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grid_grad ||_L2: " << f_grid_grad_norm_acc;
        }
        if (f_grad != 0)//eaxct solution
        {
            f_grad_norm_acc=  std::sqrt( std::accumulate( f_grad_norm.begin(), f_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grad ||_L2: " << f_grad_norm_acc;
        }
        if (fvd != 0 && f_grad != 0)//both get err
        {
            grad_err_acc=  std::sqrt( std::accumulate( grad_err.begin(), grad_err.end(), 0.));
            std::cout << "\n\t|| f_grad - f_grid_grad ||_L2: " << grad_err_acc;
        }
        std::cout << std::endl;
    }

    virtual void visit (const TetraCL& t)//accumulator main function
    {
        const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();//close cdata
        if (cdata.empty())//empty cdata, finish
            return;

        const int tid= omp_get_thread_num();//multithread

        tid0p->area[tid]+= quad_2D( cdata.qdom_projected.absdets(), cdata.qdom);//for parallelization

        std::valarray<double> qfgrid,//to store solution's dicrete value on triangle
            qf;
        if (fvd != 0)
        {
            // XXX: Check, whether incomplete P2-Data exists locally (which is allowed for P2IF_FE, but not handled correctly by this class --> Extend fvd). Likewise for P1IF_FE...
            if (fvd->RowIdx->GetFE() == P2IF_FE)
                resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
                resize_and_evaluate_on_vertexes( make_P1Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
            else if (fvd->RowIdx->GetFE() == P3IF_FE)
                resize_and_evaluate_on_vertexes( make_P3Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
//             resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), cdata.qdom_projected, qfgrid);
            tid0p->f_grid_int[tid]+=  quad_2D( cdata.qdom_projected.absdets()*qfgrid,        cdata.qdom);
            tid0p->f_grid_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*qfgrid*qfgrid, cdata.qdom);
        }
        if (f != 0)
        {
            resize_and_evaluate_on_vertexes( f, cdata.qdom_projected.vertexes(), f_time, qf);//init qf
            tid0p->f_int[tid]+= quad_2D( cdata.qdom_projected.absdets()*qf, cdata.qdom);
            tid0p->f_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*qf*qf, cdata.qdom);
        }
        if (fvd != 0 && f != 0)
        {
            std::valarray<double> qerr= qfgrid - qf;//valarray-type variant can be substrate directly
            auto tmp = quad_2D( cdata.qdom_projected.absdets()*qerr*qerr, cdata.qdom);
            tid0p->err[tid]+= tmp;
            //cout2txt(tmp);
        }




        GridFunctionCL<Point3DCL> qfgradgrid,
                       qfgrad;
        if (fvd != 0)
        {
            qfgradgrid.resize( cdata.qdom.vertex_size());
            if (fvd->RowIdx->GetFE() == P2IF_FE)
            {
                LocalP2CL<> lp2( t, *fvd, nobnddata);
                loc_lb.setup( t, cdata);
                for (Uint i= 0; i < 10; ++i)
                    qfgradgrid+= lp2[i]*loc_lb.get_qgradp2( i);
            }
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
            {
// // XXX Implement this case.
//                 LocalP1CL<> lp1( t, *fvd, nobnddata);
//                 loc_lb_p1.setup( t, cdata);
//                 for (Uint i= 0; i < 4; ++i)
//                     qfgradgrid+= lp1[i]*loc_lb_p1.get_qgradp1( i);
            }
            tid0p->f_grid_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgradgrid, qfgradgrid), cdata.qdom);
        }
        if (f_grad != 0)
        {
            resize_and_evaluate_on_vertexes( f_grad, cdata.qdom_projected.vertexes(), f_grad_time, qfgrad);
            for (Uint i= 0; i < cdata.qdom_projected.vertexes().size(); ++i)
            {
                Point3DCL n= cdata.quaqua.local_ls_grad( *cdata.qdom_projected.vertexes()[i].first, cdata.qdom_projected.vertexes()[i].second);
                n/= n.norm();
                qfgrad[i]-= inner_prod( n, qfgrad[i])*n;
            }
            tid0p->f_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgrad, qfgrad), cdata.qdom);
        }
        if (fvd != 0 && f_grad != 0)
        {
            GridFunctionCL<Point3DCL> qerr( qfgradgrid - qfgrad);
            tid0p->grad_err[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qerr, qerr), cdata.qdom);
        }
    }

    virtual InterfaceL2AccuP2CL* clone (int /*clone_id*/)
    {
        return new InterfaceL2AccuP2CL( *this);
    }
};

static LocalP2CL<double> localP2RhsCp;
LocalP1CL<Point3DCL> dfG;
//err accumulator unit, high quad version
class InterfaceL2AccuP2CLHighQuad: public TetraAccumulatorCL
{
private:
    const InterfaceCommonDataP2CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    NoBndDataCL<> nobnddata;
    const VecDescCL* fvd;

    instat_scalar_fun_ptr f;
    double f_time;

    LocalLaplaceBeltramiP2CL loc_lb;

    instat_vector_fun_ptr f_grad;
    double f_grad_time;

    InterfaceL2AccuP2CLHighQuad* tid0p; // The object in OpenMP-thread 0, in which the following variables are updated.
    std::vector<double> f_grid_norm,
        f_grid_int,
        f_norm,
        f_int,
        err,
        area,
        f_grid_grad_norm,
        f_grad_norm,
        grad_err;
    double res;

public:
    double f_grid_norm_acc,
           f_grid_int_acc,
           f_norm_acc,
           f_int_acc,
           err_acc,
           area_acc,
           f_grid_grad_norm_acc,
           f_grad_norm_acc,
           grad_err_acc;

    InterfaceL2AccuP2CLHighQuad (const InterfaceCommonDataP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), f_grad( 0), f_grad_time( 0.) {}
    virtual ~InterfaceL2AccuP2CLHighQuad () {}

    void set_name (const std::string& n)
    {
        name_= n;
    }
    void set_grid_function (const VecDescCL& fvdarg)
    {
        fvd= &fvdarg;
    }
    void set_function (const instat_scalar_fun_ptr farg, double f_time_arg= 0.)
    {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const instat_vector_fun_ptr farg, double f_time_arg= 0.)
    {
        f_grad= farg;
        f_grad_time= f_time_arg;
    }

    virtual void begin_accumulation ()
    {
        std::cout << "InterfaceL2AccuP2CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\".\n";

        tid0p= this;
        f_grid_norm.clear();
        f_grid_norm.resize( omp_get_max_threads(), 0.);
        f_grid_int.clear();
        f_grid_int.resize( omp_get_max_threads(), 0.);

        f_norm.clear();
        f_norm.resize( omp_get_max_threads(), 0.);
        f_int.clear();
        f_int.resize( omp_get_max_threads(), 0.);

        err.clear();
        err.resize( omp_get_max_threads(), 0.);
        area.clear();
        area.resize( omp_get_max_threads(), 0.);

        f_grid_grad_norm.clear();
        f_grid_grad_norm.resize( omp_get_max_threads(), 0.);
        f_grad_norm.clear();
        f_grad_norm.resize( omp_get_max_threads(), 0.);
        grad_err.clear();
        grad_err.resize( omp_get_max_threads(), 0.);

        f_grid_norm_acc= f_grid_int_acc= f_norm_acc= f_int_acc= err_acc= area_acc= 0.;
        f_grid_grad_norm_acc= f_grad_norm_acc= grad_err_acc= 0.;
    }

    virtual void finalize_accumulation()
    {
        std::cout << "InterfaceL2AccuP2CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";//"P2-solution"
        area_acc= std::accumulate( area.begin(), area.end(), 0.);
        std::cout << "\n\tarea: " << area_acc;//output area
        if (fvd != 0)//grid solution: norm and integral
        {
            f_grid_norm_acc=  std::sqrt( std::accumulate( f_grid_norm.begin(), f_grid_norm.end(), 0.));
            f_grid_int_acc= std::accumulate( f_grid_int.begin(), f_grid_int.end(), 0.);
            std::cout << "\n\t|| f_grid ||_L2: " << f_grid_norm_acc
                      << "\tintegral: " << f_grid_int_acc;
        }
        if (f != 0)//eaxct solution
        {
            f_norm_acc=  std::sqrt( std::accumulate( f_norm.begin(), f_norm.end(), 0.));
            f_int_acc= std::accumulate( f_int.begin(), f_int.end(), 0.);
            std::cout << "\n\t|| f ||_L2: " << f_norm_acc
                      << "\t integral: " << f_int_acc;
        }
        if (fvd != 0 && f != 0)//both
        {
            err_acc=  std::sqrt( std::accumulate( err.begin(), err.end(), 0.));//accumulate err
            std::cout << "\n\t|| f - f_grid ||_L2: " << err_acc;

            const double mvf_err= std::sqrt( std::pow( err_acc, 2) - std::pow( f_grid_int_acc - f_int_acc, 2)/area_acc);
            std:: cout << "\t|| f - c_f - (f_grid -c_{f_grid}) ||_L2: " << mvf_err;//some kind of "mean"
        }

        if (fvd != 0)//grid solution for gradient
        {
            f_grid_grad_norm_acc=  std::sqrt( std::accumulate( f_grid_grad_norm.begin(), f_grid_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grid_grad ||_L2: " << f_grid_grad_norm_acc;
        }
        if (f_grad != 0)//eaxct solution
        {
            f_grad_norm_acc=  std::sqrt( std::accumulate( f_grad_norm.begin(), f_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grad ||_L2: " << f_grad_norm_acc;
        }
        if (fvd != 0 && f_grad != 0)//both get err
        {
            grad_err_acc=  std::sqrt( std::accumulate( grad_err.begin(), grad_err.end(), 0.));
            std::cout << "\n\t|| f_grad - f_grid_grad ||_L2: " << grad_err_acc;
        }
        std::cout << std::endl;
    }
    static void localErrIntFunP2(double x, double y, double z, double *ff)
    {
        auto BCs = getBaryCoords(tet,x,y,z);
        double fcal = localP2RhsCp(BCs);
        DROPS::Point3DCL p(x,y,z);
        double fext = laplace_beltrami_xyz_sol(p,0);
        *ff = (fcal-fext)*(fcal-fext);

        //     DROPS::BaryCoordCL tmp{0.051566846126417189,0.07044162180172904,0.1244933792871029,
        //   0.75349815278465082};
        //   std::cout<<localP2RhsCp(tmp)<<std::endl;
        //  *ff = 1;
    }

    static void localErrIntFunH1(double x, double y, double z, double *ff)
    {
        //qfgradgridG f_grad
        //auto BCs = getBaryCoords(tet,x,y,z);

        DROPS::Point3DCL p(x,y,z);
        auto BCs = getBaryCoords(tet,x,y,z);

        Point3DCL df = laplace_beltrami_xyz_sol_grad(p,0);
        Point3DCL dfgrid = dfG(BCs);

        double n[3];
        getSfNormalVec(x,y,z,n);
        SVectorCL<3> sf_grad_f;
        SVectorCL<3> sf_grad_fgrid;
        getSurfaceGradient(df,n,sf_grad_f);
        getSurfaceGradient(dfgrid,n,sf_grad_fgrid);
        //sf_grad_fgrid = dfgrid;
        double s = 0;
        for(int i=0;i<3;i++)
        {

            s = s+ (sf_grad_f[i]-sf_grad_fgrid[i])*(sf_grad_f[i]-sf_grad_fgrid[i]);
        }
        *ff = s;

    //    *ff = dotP3(sf_grad_f-sf_grad_fgrid,sf_grad_f-sf_grad_fgrid);
 //       *ff = 0;
//        double dfcal
//
//        auto BCs = getBaryCoords(tet,x,y,z);
//        double fcal = localP2RhsCp(BCs);
//        DROPS::Point3DCL p(x,y,z);
//        double fext = laplace_beltrami_xyz_sol(p,0);
//        *ff = (fcal-fext)*(fcal-fext);

        //     DROPS::BaryCoordCL tmp{0.051566846126417189,0.07044162180172904,0.1244933792871029,
        //   0.75349815278465082};
        //   std::cout<<localP2RhsCp(tmp)<<std::endl;
        //  *ff = 1;
    }

    virtual void visit (const TetraCL& t)//accumulator main function
    {
        GetTet2DArr(t,tet);//change tetra to C-style
        const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();//close cdata
        if (cdata.empty())//empty cdata, finish
            return;

        const int tid= omp_get_thread_num();//multithread

        tid0p->area[tid]+= quad_2D( cdata.qdom_projected.absdets(), cdata.qdom);//for parallelization

        std::valarray<double> qfgrid,//to store solution's dicrete value on triangle
            qf;
        if (fvd != 0)
        {
            // XXX: Check, whether incomplete P2-Data exists locally (which is allowed for P2IF_FE, but not handled correctly by this class --> Extend fvd). Likewise for P1IF_FE...
            if (fvd->RowIdx->GetFE() == P2IF_FE)
                resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
                resize_and_evaluate_on_vertexes( make_P1Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
//             resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), cdata.qdom_projected, qfgrid);
            tid0p->f_grid_int[tid]+=  quad_2D( cdata.qdom_projected.absdets()*qfgrid,        cdata.qdom);
            tid0p->f_grid_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*qfgrid*qfgrid, cdata.qdom);
        }
        if (f != 0)
        {
            resize_and_evaluate_on_vertexes( f, cdata.qdom_projected.vertexes(), f_time, qf);//init qf
            tid0p->f_int[tid]+= quad_2D( cdata.qdom_projected.absdets()*qf, cdata.qdom);
            tid0p->f_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*qf*qf, cdata.qdom);
        }
        if (fvd != 0 && f != 0)
        {
            //std::valarray<double> qerr= qfgrid - qf;//valarray-type variant can be substrate directly
            //tid0p->err[tid]+= quad_2D( cdata.qdom_projected.absdets()*qerr*qerr, cdata.qdom);
            LocalP2CL<double> localP2Rhs(t,make_P2Eval( mg, nobnddata, *fvd));
            //      LocalP2CL<double> localP2Rhs(t,*fvd,nobnddata);
            localP2RhsCp = localP2Rhs;
            int n = phgQuadInterface2(
                        lsFun,		/* the level set function */
                        2,		/* polynomial order of the level set function */
                        lsGrad,	/* the gradient of the level set function */
                        tet,		/* coordinates of the vertices of the tetra */
                        localErrIntFunP2,		/* the integrand */
                        1,		/* dimension of the integrand */
                        DOF_PROJ_NONE,	/* projection type for surface integral */
                        0,		/* integration type (-1, 0, 1) */
                        orderG,		/* order of the 1D Gaussian quadrature */
                        &res,		/* the computed integral */
                        NULL		/* pointer returning the computed rule */
                    );


            tid0p->err[tid]+= res;

            //cout2txt(res);

        }
        GridFunctionCL<Point3DCL> qfgradgrid;
        GridFunctionCL<Point3DCL> qfgrad;
        if (fvd != 0)
        {
            qfgradgrid.resize( cdata.qdom.vertex_size());
            if (fvd->RowIdx->GetFE() == P2IF_FE)
            {
                LocalP2CL<> lp2( t, *fvd, nobnddata);
                loc_lb.setup( t, cdata);
                for (Uint i= 0; i < 10; ++i)
                    qfgradgrid+= lp2[i]*loc_lb.get_qgradp2( i);
            }
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
            {
// // XXX Implement this case.
//                 LocalP1CL<> lp1( t, *fvd, nobnddata);
//                 loc_lb_p1.setup( t, cdata);
//                 for (Uint i= 0; i < 4; ++i)
//                     qfgradgrid+= lp1[i]*loc_lb_p1.get_qgradp1( i);
            }
            tid0p->f_grid_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgradgrid, qfgradgrid), cdata.qdom);
        }
        if (f_grad != 0)
        {
            resize_and_evaluate_on_vertexes( f_grad, cdata.qdom_projected.vertexes(), f_grad_time, qfgrad);
            for (Uint i= 0; i < cdata.qdom_projected.vertexes().size(); ++i)
            {
                Point3DCL n= cdata.quaqua.local_ls_grad( *cdata.qdom_projected.vertexes()[i].first, cdata.qdom_projected.vertexes()[i].second);
                n/= n.norm();
                qfgrad[i]-= inner_prod( n, qfgrad[i])*n;
            }
            tid0p->f_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgrad, qfgrad), cdata.qdom);
        }
        if (fvd != 0 && f_grad != 0)
        {
            LocalP2CL<> localP2Fcal(t,make_P2Eval( mg, nobnddata, *fvd));
            //LocalP2CL<> localP2Fcal( t, *fvd, nobnddata);
           // LocalP1CL<Point3DCL> df;
            loc_lb.setup( t, cdata);
            auto df = localP2Fcal[0]*loc_lb.get_gradp2( 0);
            for (Uint i= 1; i < 10; ++i)
            {
                df = df + localP2Fcal[i]*loc_lb.get_gradp2( i);
            }
            dfG = df;


            int n = phgQuadInterface2(
                        lsFun,		/* the level set function */
                        2,		/* polynomial order of the level set function */
                        lsGrad,	/* the gradient of the level set function */
                        tet,		/* coordinates of the vertices of the tetra */
                        localErrIntFunH1,		/* the integrand */
                        1,		/* dimension of the integrand */
                        DOF_PROJ_NONE,	/* projection type for surface integral */
                        0,		/* integration type (-1, 0, 1) */
                        orderG,		/* order of the 1D Gaussian quadrature */
                        &res,		/* the computed integral */
                        NULL		/* pointer returning the computed rule */
                    );


            tid0p->grad_err[tid] += res;

        }
    }

    virtual InterfaceL2AccuP2CLHighQuad* clone (int /*clone_id*/)
    {
        return new InterfaceL2AccuP2CLHighQuad( *this);
    }
};


/// \brief Accumulate different error measures for the approximation of the level
/// set function $\varphi$ by the piecewise linear approximation $\varphi_h$.
class InterfaceApproxErrorDeformAccuCL : public TetraAccumulatorCL
{
private:
    VecDescCL* yG_, ///< error-term for H^1 smooth \varphi
               * ydist_; ///< distance per QuaQuamapper

    InterfaceCommonDataDeformP2CL* cdata_;

    IdxT numry[4];
    double vec[4];

    instat_scalar_fun_ptr d_; // exact distance
    instat_vector_fun_ptr Dd_; // gradient of exact distance
    GridFunctionCL<> qd;
    GridFunctionCL<Point3DCL> qDderr;

    OpenMPVar_MinInit_Max_CL<double> max_h,
                             max_d,
                             max_Dderr;

public:
    InterfaceApproxErrorDeformAccuCL (InterfaceCommonDataDeformP2CL* cdataarg, VecDescCL* yg, VecDescCL* ydist)
        : yG_( yg), ydist_( ydist), cdata_ (cdataarg), d_ (0), Dd_ (0) {}
    virtual ~InterfaceApproxErrorDeformAccuCL () {}

    InterfaceApproxErrorDeformAccuCL& set_d  (instat_scalar_fun_ptr darg)
    {
        d_= darg;
        return *this;
    }
    InterfaceApproxErrorDeformAccuCL& set_Dd (instat_vector_fun_ptr Ddarg)
    {
        Dd_= Ddarg;
        return *this;
    }

    virtual void begin_accumulation ()
    {
        std::cout << "#InterfaceApproxErrorDeformAccuCL::begin_accumulation"
                  ": " << ydist_->RowIdx->NumUnknowns() << " rows.\n";
        max_h.scatter ();
        max_d.scatter ();
        max_Dderr.scatter ();
    }

    virtual void finalize_accumulation()
    {
        std::cout << "#InterfaceApproxErrorDeformAccuCL::finalize_accumulation: ";
        max_h.reduce();
        std::cout << "\n\tmax_h: " << max_h.value() << "\n";
        max_d.reduce();
        std::cout << "\tmax_d: " << max_d.value() << "\n";
        max_Dderr.reduce();
        std::cout << "\tmax_dDerr: " << max_Dderr.value() << "\n";
    }

    virtual void visit (const TetraCL& t)
    {
        InterfaceCommonDataDeformP2CL& cdata= cdata_->get_clone ();
        const int tid= omp_get_thread_num();
        if (cdata.empty ())
            return;

        resize_and_evaluate_on_vertexes( d_, t, cdata.qdom2d_full, 0., qd);
        qd= std::abs (qd);
        max_d.value (tid)= std::max (max_d.value (tid), *std::max_element (&qd[0], &qd[0] + cdata.qdom2d_full.vertex_size ()));
        resize_and_evaluate_on_vertexes( Dd_, t, cdata.qdom2d_full, 0., qDderr);
        const SurfacePatchCL::FacetT& facet= cdata.surf.facet_begin()[0];
        const BaryCoordCL verts[3]= { cdata.surf.vertex_begin()[facet[0]],
                                      cdata.surf.vertex_begin()[facet[1]],
                                      cdata.surf.vertex_begin()[facet[2]]
                                    };
        cdata.Phi.set_surface_patch (verts, cdata.pos_pt);
        for (Uint i= 0; i < cdata.qdom2d_full.vertex_size (); ++i)
        {
            cdata.Phi.set_point (cdata.qdom2d_only_weights.vertex_begin ()[i], true);
            Point3DCL n= cdata.Phi.dPhi(cdata.qdom2d_only_weights.vertex_begin ()[i])*cdata.Phi.w;
            n/=n.norm();
            qDderr[i]-= n;
            max_Dderr.value (tid)= std::max (max_Dderr.value (tid), qDderr[i].norm ());
        }

//         for (Uint d= 0; d < 10; ++d) {
//             G+= cdata.gradp2[d]( bc)*cdata.locp2_ls[d];
//         }
        max_h.value (tid)= std::max (max_h.value (tid), ::cbrt( std::abs( cdata.det_T)));
        GetLocalNumbP1NoBnd( numry, t, *ydist_->RowIdx);
        for (int i= 0; i < 4; ++i)
        {
            ydist_->Data[numry[i]]= std::max (ydist_->Data[numry[i]], max_d.value (tid));
//             yG_->Data[numry[i]]= std::max( yG_->Data[numry[i]], G.norm ());
        }
    }

    virtual InterfaceApproxErrorDeformAccuCL* clone (int /*clone_id*/)
    {
        InterfaceApproxErrorDeformAccuCL* p= new InterfaceApproxErrorDeformAccuCL( *this);
        p->max_h.make_reference_to (max_h);
        p->max_d.make_reference_to (max_d);
        p->max_Dderr.make_reference_to (max_Dderr);
        return p;
    }
};

/// \brief Accumulate L2-norms and errors on the deformed zero level.
/// Works for P1IF_FE, P2IF_FE, and C-functions. All functions are evaluated on the deformed levelset.
class InterfaceL2AccuDeformP2CL : public TetraAccumulatorCL
{
private:
    const InterfaceCommonDataDeformP2CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    NoBndDataCL<> nobnddata;
    const VecDescCL* fvd;

    instat_scalar_fun_ptr f;
    double f_time;

    LocalLaplaceBeltramiDeformP2CL loc_lb;
    LocalNormalLaplaceDeformP2CL loc_ngrad;

    instat_vector_fun_ptr f_grad;
    double f_grad_time;

    InterfaceL2AccuDeformP2CL* tid0p; // The object in OpenMP-thread 0, in which the following variables are updated.
    std::vector<double> f_grid_norm,
        f_grid_int,
        f_norm,
        f_int,
        err,
        area,
        f_grid_grad_norm,
        f_grad_norm,
        grad_err,
        normal_grad;

public:
    double f_grid_norm_acc,
           f_grid_int_acc,
           f_norm_acc,
           f_int_acc,
           err_acc,
           area_acc,
           f_grid_grad_norm_acc,
           f_grad_norm_acc,
           grad_err_acc,
           normal_grad_acc;

    InterfaceL2AccuDeformP2CL (const InterfaceCommonDataDeformP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), loc_ngrad (1.), f_grad( 0), f_grad_time( 0.) {}
    virtual ~InterfaceL2AccuDeformP2CL () {}

    void set_name (const std::string& n)
    {
        name_= n;
    }
    void set_grid_function (const VecDescCL& fvdarg)
    {
        fvd= &fvdarg;
    }
    void set_function (const instat_scalar_fun_ptr farg, double f_time_arg= 0.)
    {
        f= farg;
        f_time= f_time_arg;
    }
    void set_grad_function (const instat_vector_fun_ptr farg, double f_time_arg= 0.)
    {
        f_grad= farg;
        f_grad_time= f_time_arg;
    }

    virtual void begin_accumulation ()
    {
        std::cout << "InterfaceL2AccuDeformP2CL::begin_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\".\n";

        tid0p= this;
        f_grid_norm.clear();
        f_grid_norm.resize( omp_get_max_threads(), 0.);
        f_grid_int.clear();
        f_grid_int.resize( omp_get_max_threads(), 0.);

        f_norm.clear();
        f_norm.resize( omp_get_max_threads(), 0.);
        f_int.clear();
        f_int.resize( omp_get_max_threads(), 0.);

        err.clear();
        err.resize( omp_get_max_threads(), 0.);
        area.clear();
        area.resize( omp_get_max_threads(), 0.);

        f_grid_grad_norm.clear();
        f_grid_grad_norm.resize( omp_get_max_threads(), 0.);
        f_grad_norm.clear();
        f_grad_norm.resize( omp_get_max_threads(), 0.);
        grad_err.clear();
        grad_err.resize( omp_get_max_threads(), 0.);
        normal_grad.clear();
        normal_grad.resize( omp_get_max_threads(), 0.);

        f_grid_norm_acc= f_grid_int_acc= f_norm_acc= f_int_acc= err_acc= area_acc= 0.;
        f_grid_grad_norm_acc= f_grad_norm_acc= grad_err_acc= 0.;
    }

    virtual void finalize_accumulation()
    {
        std::cout << "InterfaceL2AccuDeformP2CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";
        area_acc= std::accumulate( area.begin(), area.end(), 0.);
        std::cout << "\n\tarea: " << area_acc;
        if (fvd != 0)
        {
            f_grid_norm_acc=  std::sqrt( std::accumulate( f_grid_norm.begin(), f_grid_norm.end(), 0.));
            f_grid_int_acc= std::accumulate( f_grid_int.begin(), f_grid_int.end(), 0.);
            std::cout << "\n\t|| f_grid ||_L2: " << f_grid_norm_acc
                      << "\tintegral: " << f_grid_int_acc;
        }
        if (f != 0)
        {
            f_norm_acc=  std::sqrt( std::accumulate( f_norm.begin(), f_norm.end(), 0.));
            f_int_acc= std::accumulate( f_int.begin(), f_int.end(), 0.);
            std::cout << "\n\t|| f ||_L2: " << f_norm_acc
                      << "\t integral: " << f_int_acc;
        }
        if (fvd != 0 && f != 0)
        {
            err_acc=  std::sqrt( std::accumulate( err.begin(), err.end(), 0.));
            std::cout << "\n\t|| f - f_grid ||_L2: " << err_acc;

            const double mvf_err= std::sqrt( std::pow( err_acc, 2) - std::pow( f_grid_int_acc - f_int_acc, 2)/area_acc);
            std:: cout << "\t|| f - c_f - (f_grid -c_{f_grid}) ||_L2: " << mvf_err;
        }

        if (fvd != 0)
        {
            f_grid_grad_norm_acc=  std::sqrt( std::accumulate( f_grid_grad_norm.begin(), f_grid_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grid_grad ||_L2: " << f_grid_grad_norm_acc;
        }
        if (f_grad != 0)
        {
            f_grad_norm_acc=  std::sqrt( std::accumulate( f_grad_norm.begin(), f_grad_norm.end(), 0.));
            std::cout << "\n\t|| f_grad ||_L2: " << f_grad_norm_acc;
        }
        if (fvd != 0 && f_grad != 0)
        {
            grad_err_acc=  std::sqrt( std::accumulate( grad_err.begin(), grad_err.end(), 0.));
            std::cout << "\n\t|| f_grad - f_grid_grad ||_L2: " << grad_err_acc;
        }
        if (fvd != 0)
        {
            normal_grad_acc=  std::sqrt( std::accumulate( normal_grad.begin(), normal_grad.end(), 0.));
            std::cout << "\n\t|| n^Tf_grid_grad ||_L2(vol): " << normal_grad_acc;
        }
        std::cout << std::endl;
    }

    virtual void visit (const TetraCL& t)
    {
        const InterfaceCommonDataDeformP2CL& cdata= cdata_.get_clone();
        if (cdata.empty())
            return;

        const int tid= omp_get_thread_num();

        std::valarray<double> ones (1., cdata.qdom2d_only_weights.vertex_size ());
        tid0p->area[tid]+= quad_2D( ones, cdata.qdom2d_only_weights);

        std::valarray<double> qfgrid,
            qf;
        if (fvd != 0)
        {
            // XXX: Check, whether incomplete P2-Data exists locally (which is allowed for P2IF_FE, but not handled correctly by this class --> Extend fvd). Likewise for P1IF_FE...
            if (fvd->RowIdx->GetFE() == P2IF_FE)
                resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), t, cdata.qdom2d_only_weights, qfgrid);
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
                resize_and_evaluate_on_vertexes( make_P1Eval( mg, nobnddata, *fvd), t, cdata.qdom2d_only_weights, qfgrid);
//             resize_and_evaluate_on_vertexes( make_P2Eval( mg, nobnddata, *fvd), cdata.qdom_projected, qfgrid);
            tid0p->f_grid_int[tid]+=  quad_2D( qfgrid,        cdata.qdom2d_only_weights);
            tid0p->f_grid_norm[tid]+= quad_2D( qfgrid*qfgrid, cdata.qdom2d_only_weights);
        }
        if (f != 0)
        {
            resize_and_evaluate_on_vertexes( f, t, cdata.qdom2d_full, f_time, qf);
            tid0p->f_int[tid]+= quad_2D( qf, cdata.qdom2d_full);
            tid0p->f_norm[tid]+= quad_2D( qf*qf, cdata.qdom2d_full);
        }
        if (fvd != 0 && f != 0)
        {
            std::valarray<double> qerr= qfgrid - qf;
            tid0p->err[tid]+= quad_2D( qerr*qerr, cdata.qdom2d_full);
        }

        if (fvd != 0)
        {
            if (fvd->RowIdx->GetFE() == P2IF_FE)
            {
                LocalP2CL<> lp2( t, *fvd, nobnddata);
                loc_ngrad.setup( t, cdata);
                for (Uint i= 0; i < 10; ++i)
                    for (Uint j= 0; j < 10; ++j)
                        tid0p->normal_grad[tid]+= lp2[i]*lp2[j]*loc_ngrad.coup[i][j];
            }
        }
//         GridFunctionCL<Point3DCL> qfgradgrid,
//                                   qfgrad;
//         if (fvd != 0) {
//             qfgradgrid.resize( cdata.qdom.vertex_size());
//             if (fvd->RowIdx->GetFE() == P2IF_FE) {
//                 LocalP2CL<> lp2( t, *fvd, nobnddata);
//                 loc_lb.setup( t, cdata);
//                 for (Uint i= 0; i < 10; ++i)
//                     qfgradgrid+= lp2[i]*loc_lb.get_qgradp2( i);
//             }
//             else if (fvd->RowIdx->GetFE() == P1IF_FE) {
// // // XXX Implement this case.
// //                 LocalP1CL<> lp1( t, *fvd, nobnddata);
// //                 loc_lb_p1.setup( t, cdata);
// //                 for (Uint i= 0; i < 4; ++i)
// //                     qfgradgrid+= lp1[i]*loc_lb_p1.get_qgradp1( i);
//             }
//             tid0p->f_grid_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgradgrid, qfgradgrid), cdata.qdom);
//         }
//         if (f_grad != 0) {
//             resize_and_evaluate_on_vertexes( f_grad, cdata.qdom_projected.vertexes(), f_grad_time, qfgrad);
//             for (Uint i= 0; i < cdata.qdom_projected.vertexes().size(); ++i) {
//                 Point3DCL n= cdata.quaqua.local_ls_grad( *cdata.qdom_projected.vertexes()[i].first, cdata.qdom_projected.vertexes()[i].second);
//                 n/= n.norm();
//                 qfgrad[i]-= inner_prod( n, qfgrad[i])*n;
//             }
//             tid0p->f_grad_norm[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qfgrad, qfgrad), cdata.qdom);
//         }
//         if (fvd != 0 && f_grad != 0) {
//             GridFunctionCL<Point3DCL> qerr( qfgradgrid - qfgrad);
//             tid0p->grad_err[tid]+= quad_2D( cdata.qdom_projected.absdets()*dot( qerr, qerr), cdata.qdom);
//         }
    }

    virtual InterfaceL2AccuDeformP2CL* clone (int /*clone_id*/)
    {
        return new InterfaceL2AccuDeformP2CL( *this);
    }
};
#define P2HIGH
void StationaryStrategyP2 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
//std::cout << P << std::endl;
    // Initialize level set and triangulation
    adap.MakeInitialTriang();//init mg
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);//level set numbering
    lset.Phi.SetIdx( &lset.idx);//set discrete lsvel set's index
    // LinearLSInit( mg, lset.Phi, &the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);//initial level set function

    // Setup an interface-P2 numbering
    DROPS::IdxDescCL ifacep2idx( P2IF_FE); //p2 element index decription class
    //std::cout<<"here:"<<P.get<double>("SurfTransp.XFEMReduced")<<std::endl;
    //ifacep2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));//set boundary, only for XFEM
    ifacep2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());//consider boundary
    std::cout << "P2-NumUnknowns: " << ifacep2idx.NumUnknowns() << std::endl;

    // Recover the gradient of the level set function
    IdxDescCL vecp2idx( vecP2_FE);//vector p2 index
    vecp2idx.CreateNumbering( mg.GetLastLevel(), mg);//creat globle index
    VecDescCL lsgradrec( &vecp2idx);//level set gradient recorver
    averaging_P2_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);

    // Compute neighborhoods of the tetras at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);//lattice
    TetraToTetrasT tetra_neighborhoods;//to put tetra neighborhoods
    compute_tetra_neighborhoods( mg, lset.Phi, lset.GetBndData(), lat, tetra_neighborhoods);//compute neighborhood

    QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, &tetra_neighborhoods,
                           P.get<int>( "LevelsetMapper.Iter"),
                           P.get<double>( "LevelsetMapper.Tol"),
                           P.get<std::string>( "LevelsetMapper.Method") == "FixedPointWithLineSearch",
                           P.get<double>( "LevelsetMapper.ArmijoConstant"));

    VecDescCL to_iface( &vecp2idx);
//     {
//         TetraAccumulatorTupleCL accus;
//         InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);
//         accus.push_back( &cdatap2);
//         InterfaceDebugP2CL p2debugaccu( cdatap2);
// //         p2debugaccu.store_offsets( to_iface);
//         p2debugaccu.set_true_area( 4.*M_PI*RadDrop[0]*RadDrop[0]);
//         p2debugaccu.set_ref_dp( &dp_sphere);
//         p2debugaccu.set_ref_abs_det( &abs_det_sphere);
//         accus.push_back( &p2debugaccu);
//         accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     }

    TetraAccumulatorTupleCL accus;//init an accumulator
    InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);//store p2 data
    accus.push_back( &cdatap2);//push_back P2 data

    //set up mass matrix
    DROPS::MatDescCL Mp2( &ifacep2idx, &ifacep2idx);//mass matrix
    //111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
#ifndef P2HIGH
    InterfaceMatrixAccuCL<LocalMassP2CL, InterfaceCommonDataP2CL> accuMp2( &Mp2, LocalMassP2CL(), cdatap2, "Mp2");
#else
    InterfaceMatrixAccuCL<LocalMassP2CLHighQuad, InterfaceCommonDataP2CL> accuMp2( &Mp2, LocalMassP2CLHighQuad(), cdatap2, "Mp2");
#endif
    accus.push_back( &accuMp2);

    //set up stiffness matrix
    DROPS::MatDescCL Ap2( &ifacep2idx, &ifacep2idx);
    //22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
#ifndef P2HIGH
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP2CL, InterfaceCommonDataP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiP2CL( P.get<double>("SurfTransp.Visc")), cdatap2, "Ap2");
#else
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP2CLHighQuad, InterfaceCommonDataP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiP2CLHighQuad( P.get<double>("SurfTransp.Visc")), cdatap2, "Ap2");
#endif
    accus.push_back( &accuAp2);
    //set up right hand side
    DROPS::VecDescCL bp2( &ifacep2idx);
    //33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
#ifndef P2HIGH
    InterfaceVectorAccuCL<LocalVectorP2CL, InterfaceCommonDataP2CL> acculoadp2( &bp2, LocalVectorP2CL( the_rhs_fun, bp2.t), cdatap2);
#else
    InterfaceVectorAccuCL<LocalVectorP2CLHighQuad, InterfaceCommonDataP2CL> acculoadp2( &bp2, LocalVectorP2CLHighQuad( the_rhs_fun, bp2.t), cdatap2);
#endif
    accus.push_back( &acculoadp2);
    accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());//begin tetra loop

//     TetraAccumulatorTupleCL mean_accus;
//     mean_accus.push_back( &cdatap2);
//     InterfaceL2AccuP2CL L2_mean_accu( cdatap2, mg, "P2-mean");
//     L2_mean_accu.set_grid_function( bp2);
//     mean_accus.push_back( &L2_mean_accu);
//     accumulate( mean_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     bp2.Data-= L2_mean_accu.f_grid_int_acc/L2_mean_accu.area_acc;
//     accumulate( mean_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());

//     VectorCL e( 1., bp2.Data.size());
//     VectorCL Ldiag( Ap2.Data.GetDiag());
//     bp2.Data-= dot( VectorCL( e/Ldiag), bp2.Data)/std::sqrt( dot( VectorCL( e/Ldiag), e));

//left hand matrix
    DROPS::MatrixCL Lp2;
    Lp2.LinComb( 1.0, Ap2.Data, 1.0, Mp2.Data);
//   MatrixCL& Lp2= Ap2.Data;
// keep data
#if 1
    DROPS::WriteToFile( Ap2.Data, "ap2_iface.txt", "Ap2");
    DROPS::WriteToFile( Mp2.Data, "mp2_iface.txt", "Mp2");
    DROPS::WriteFEToFile( bp2, mg, "rhsp2_iface.txt", /*binary=*/ false);
#endif

#if 1
//define solver and solve linear equations
    typedef DROPS::SSORPcCL SurfPcT;
//     typedef DROPS::JACPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
    DROPS::VecDescCL xp2( &ifacep2idx);
    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
    DROPS::WriteFEToFile( xp2, mg, "xp2_iface.txt", /*binary=*/ false);
#endif

#if 1

//define direct solver
# if 0
    DROPS::VecDescCL xp2( &ifacep2idx);
    DROPS::DirectSymmSolverCL dsolver(Lp2);
    dsolver.Solve(Lp2,xp2.Data,bp2.Data);
    dsolver.Update(Lp2);
    dsolver.Solve(Lp2,xp2.Data,bp2.Data);
#endif
    //dsolver.Update(A);
    //dsolver.Solve(A,x,b);
    //DROPS::DirectNonSymmSolverCL dnsolver(Lp2);
    //dnsolver.Solve(A,x,b);
# if 0
    DROPS::VecDescCL xp2( &ifacep2idx);
    DROPS::DirectNonSymmSolverCL dnsolver(Lp2);
    dnsolver.Solve(Lp2,xp2.Data,bp2.Data);
#endif

#endif

//define solver and solve linear equations
//    typedef DROPS::SSORPcCL SurfPcT;
//     typedef DROPS::JACPcCL SurfPcT;
//    SurfPcT surfpc;
//    typedef DROPS::GMResSolverCL<SurfPcT> SurfSolverT;
//    SurfSolverT surfsolver( surfpc,10, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
//    DROPS::VecDescCL xp2( &ifacep2idx);
//    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
//    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
//    DROPS::WriteFEToFile( xp2, mg, "xp2_iface.txt", /*binary=*/ false);


//define solver and solve linear equations
//   typedef DROPS::SSORPcCL SurfPcT;
//    typedef DROPS::JACPcCL SurfPcT;
//    SurfPcT surfpc;
//    typedef DROPS::GCRSolverCL<SurfPcT> SurfSolverT;
//    SurfSolverT surfsolver( surfpc,10, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
//    DROPS::VecDescCL xp2( &ifacep2idx);
//    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
//    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
//    DROPS::WriteFEToFile( xp2, mg, "xp2_iface.txt", /*binary=*/ false);



//accumulate errors on every tetrahedron
    TetraAccumulatorTupleCL err_accus;//final tetra error accumulator
    err_accus.push_back( &cdatap2);//push back cdata, include P2 element
    //444444444444444444444444444444444444444444444444444444444444444444444444444444444444444

//       InterfaceL2AccuP2CLHighQuad (const InterfaceCommonDataP2CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
    //      : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), f_grad( 0), f_grad_time( 0.) {}


#ifndef P2HIGH
    InterfaceL2AccuP2CL L2_accu( cdatap2, mg, "P2-solution");//interface L2 accumulator
#else
    InterfaceL2AccuP2CLHighQuad L2_accu( cdatap2, mg, "P2-solution");//error accumulater high quad version
#endif
    L2_accu.set_grid_function( xp2);//set solved solution
    L2_accu.set_function( the_sol_fun, 0.);//set exaction solution
    L2_accu.set_grad_function( the_sol_grad_fun, 0.);//set gradient for exact sol
    err_accus.push_back( &L2_accu);
    accumulate( err_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());//begin accumulating
    DROPS::NoBndDataCL<> nobnd;
    VecDescCL the_sol_vd( &lset.idx);
    //LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    //double H1_err= H1_error_p2( lset.Phi, lset.GetBndData(), make_P2Eval( mg, nobnd, xp2), the_sol_vd,ifacep2idx,Ap2);
    std::cout << "L_2-error: " << L2_accu.err_acc<<std::endl;
//    std::cout << "H_1-error: "<< H1_err <<std::endl;
    std::cout << "H_1-error: " << L2_accu.grad_err_acc<<std::endl;


#if 1
//write out
    {
        std::ofstream os( "quaqua_num_outer_iter.txt");
        for (Uint i= 0; i != quaqua.num_outer_iter.size(); ++i)
            os << i << '\t' << quaqua.num_outer_iter[i] << '\n';
        os << '\n';
        for (Uint i= 0; i != quaqua.num_inner_iter.size(); ++i)
            os << i << '\t' << quaqua.num_inner_iter[i] << '\n';
    }
    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)//write to file
        DROPS::WriteFEToFile( xp2, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path") + "_p2", P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    if (vtkwriter.get() != 0)
    {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, xp2, "InterfaceSolP2"));
        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, lsgradrec), "LSGradRec") );
        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, to_iface), "to_iface") );
        vtkwriter->Write( 0.);
    }
#endif
}


#define P3HIGH0
void StationaryStrategyP3 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
//std::cout << P << std::endl;
    // Initialize level set and triangulation
    adap.MakeInitialTriang();//init mg
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);//level set numbering
    lset.Phi.SetIdx( &lset.idx);//set discrete lsvel set's index
    // LinearLSInit( mg, lset.Phi, &the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);//initial level set function

    // Setup an interface-P3 numbering
    DROPS::IdxDescCL ifacep3idx( P3IF_FE); //p3 element index decription class
    //std::cout<<"here:"<<P.get<double>("SurfTransp.XFEMReduced")<<std::endl;
    //ifacep3idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));//set boundary, only for XFEM
    ifacep3idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());//consider boundary
    std::cout << "P3-NumUnknowns: " << ifacep3idx.NumUnknowns() << std::endl;

    // Recover the gradient of the level set function
    IdxDescCL vecp3idx( vecP3_FE);//vector p3 index
    vecp3idx.CreateNumbering( mg.GetLastLevel(), mg);//creat globle index
    VecDescCL lsgradrec( &vecp3idx);//level set gradient recorver
    averaging_P3_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);//MultiGridCL

    // Compute neighborhoods of the tetras at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);//lattice
    TetraToTetrasT tetra_neighborhoods;//to put tetra neighborhoods
    compute_tetra_neighborhoods( mg, lset.Phi, lset.GetBndData(), lat, tetra_neighborhoods);//compute neighborhood

    QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, &tetra_neighborhoods,
                           P.get<int>( "LevelsetMapper.Iter"),
                           P.get<double>( "LevelsetMapper.Tol"),
                           P.get<std::string>( "LevelsetMapper.Method") == "FixedPointWithLineSearch",
                           P.get<double>( "LevelsetMapper.ArmijoConstant"));

    VecDescCL to_iface( &vecp3idx);
//     {
//         TetraAccumulatorTupleCL accus;
//         InterfaceCommonDataP3CL cdatap3( lset.Phi, lset.GetBndData(), quaqua, lat);
//         accus.push_back( &cdatap3);
//         InterfaceDebugP3CL p3debugaccu( cdatap3);
// //         p3debugaccu.store_offsets( to_iface);
//         p3debugaccu.set_true_area( 4.*M_PI*RadDrop[0]*RadDrop[0]);
//         p3debugaccu.set_ref_dp( &dp_sphere);
//         p3debugaccu.set_ref_abs_det( &abs_det_sphere);
//         accus.push_back( &p3debugaccu);
//         accumulate( accus, mg, ifacep3idx.TriangLevel(), ifacep3idx.GetMatchingFunction(), ifacep3idx.GetBndInfo());
//     }

    TetraAccumulatorTupleCL accus;//init an accumulator
    InterfaceCommonDataP3CL cdatap3( lset.Phi, lset.GetBndData(), quaqua, lat);//store p3 data InterfaceCommonDataP2CL
    accus.push_back( &cdatap3);//push_back P3 data

    //set up mass matrix
    DROPS::MatDescCL Mp3( &ifacep3idx, &ifacep3idx);//mass matrix
    //111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
#ifndef P3HIGH
    InterfaceMatrixAccuCL<LocalMassP3CL, InterfaceCommonDataP3CL> accuMp3( &Mp3, LocalMassP3CL(), cdatap3, "Mp3");//LocalMassP2CL
#else
    InterfaceMatrixAccuCL<LocalMassP3CLHighQuad, InterfaceCommonDataP3CL> accuMp3( &Mp3, LocalMassP3CLHighQuad(), cdatap3, "Mp3");
#endif
    accus.push_back( &accuMp3);

    //set up stiffness matrix
    DROPS::MatDescCL Ap3( &ifacep3idx, &ifacep3idx);
    //22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
#ifndef P3HIGH
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP3CL, InterfaceCommonDataP3CL> accuAp3( &Ap3, LocalLaplaceBeltramiP3CL( P.get<double>("SurfTransp.Visc")), cdatap3, "Ap3");//LocalLaplaceBeltramiP2CL
#else
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiP3CLHighQuad, InterfaceCommonDataP3CL> accuAp3( &Ap3, LocalLaplaceBeltramiP3CLHighQuad( P.get<double>("SurfTransp.Visc")), cdatap3, "Ap3");
    //LocalLaplaceBeltramiP2CLHighQuad
#endif
    accus.push_back( &accuAp3);
    //set up right hand side
    DROPS::VecDescCL bp3( &ifacep3idx);
    //33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
#ifndef P3HIGH
    InterfaceVectorAccuCL<LocalVectorP3CL, InterfaceCommonDataP3CL> acculoadp3( &bp3, LocalVectorP3CL( the_rhs_fun, bp3.t), cdatap3);//LocalVectorP2CL
#else
    InterfaceVectorAccuCL<LocalVectorP3CLHighQuad, InterfaceCommonDataP3CL> acculoadp3( &bp3, LocalVectorP3CLHighQuad( the_rhs_fun, bp3.t), cdatap3);//LocalVectorP2CLHighQuad
#endif
    accus.push_back( &acculoadp3);
    accumulate( accus, mg, ifacep3idx.TriangLevel(), ifacep3idx.GetBndInfo());//begin tetra loop

//     TetraAccumulatorTupleCL mean_accus;
//     mean_accus.push_back( &cdatap3);
//     InterfaceL3AccuP3CL L3_mean_accu( cdatap3, mg, "P3-mean");
//     L3_mean_accu.set_grid_function( bp3);
//     mean_accus.push_back( &L3_mean_accu);
//     accumulate( mean_accus, mg, ifacep3idx.TriangLevel(), ifacep3idx.GetMatchingFunction(), ifacep3idx.GetBndInfo());
//     bp3.Data-= L3_mean_accu.f_grid_int_acc/L3_mean_accu.area_acc;
//     accumulate( mean_accus, mg, ifacep3idx.TriangLevel(), ifacep3idx.GetMatchingFunction(), ifacep3idx.GetBndInfo());

//     VectorCL e( 1., bp3.Data.size());
//     VectorCL Ldiag( Ap3.Data.GetDiag());
//     bp3.Data-= dot( VectorCL( e/Ldiag), bp3.Data)/std::sqrt( dot( VectorCL( e/Ldiag), e));

//left hand matrix
    DROPS::MatrixCL Lp3;
    Lp3.LinComb( 1.0, Ap3.Data, 1.0, Mp3.Data);
//   MatrixCL& Lp3= Ap3.Data;
// keep data
#if 1
    DROPS::WriteToFile( Ap3.Data, "ap3_iface.txt", "Ap3");
    DROPS::WriteToFile( Mp3.Data, "mp3_iface.txt", "Mp3");
    DROPS::WriteFEToFile( bp3, mg, "rhsp3_iface.txt", /*binary=*/ false);
#endif


//define solver and solve linear equations
    typedef DROPS::SSORPcCL SurfPcT;
//     typedef DROPS::JACPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
    DROPS::VecDescCL xp3( &ifacep3idx);
    surfsolver.Solve( Lp3, xp3.Data, bp3.Data, xp3.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
    DROPS::WriteFEToFile( xp3, mg, "xp3_iface.txt", /*binary=*/ false);


//define solver and solve linear equations
//    typedef DROPS::SSORPcCL SurfPcT;
//     typedef DROPS::JACPcCL SurfPcT;
//    SurfPcT surfpc;
//    typedef DROPS::GMResSolverCL<SurfPcT> SurfSolverT;
//    SurfSolverT surfsolver( surfpc,10, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
//    DROPS::VecDescCL xp3( &ifacep3idx);
//    surfsolver.Solve( Lp3, xp3.Data, bp3.Data, xp3.RowIdx->GetEx());
//    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
//    DROPS::WriteFEToFile( xp3, mg, "xp3_iface.txt", /*binary=*/ false);


//define solver and solve linear equations
//   typedef DROPS::SSORPcCL SurfPcT;
//    typedef DROPS::JACPcCL SurfPcT;
//    SurfPcT surfpc;
//    typedef DROPS::GCRSolverCL<SurfPcT> SurfSolverT;
//    SurfSolverT surfsolver( surfpc,10, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);
//    DROPS::VecDescCL xp3( &ifacep3idx);
//    surfsolver.Solve( Lp3, xp3.Data, bp3.Data, xp3.RowIdx->GetEx());
//    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';
//    DROPS::WriteFEToFile( xp3, mg, "xp3_iface.txt", /*binary=*/ false);



//accumulate errors on every tetrahedron
    TetraAccumulatorTupleCL err_accus;//final tetra error accumulator
    err_accus.push_back( &cdatap3);//push back cdata, include P3 element
    //444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
#ifndef P3HIGH
    InterfaceL2AccuP3CL L3_accu( cdatap3, mg, "P3-solution");//interface L3 accumulator InterfaceL2AccuP2CL
#else
    InterfaceL3AccuP3CLHighQuad L3_accu( cdatap3, mg, "P3-solution");//error accumulater high quad version
#endif
    L3_accu.set_grid_function( xp3);//set solved solution
    L3_accu.set_function( the_sol_fun, 0.);//set exaction solution
    L3_accu.set_grad_function( the_sol_grad_fun, 0.);//set gradient for exact sol
    err_accus.push_back( &L3_accu);
    accumulate( err_accus, mg, ifacep3idx.TriangLevel(), ifacep3idx.GetBndInfo());//begin accumulating
#if 1
//write out
    {
        std::ofstream os( "quaqua_num_outer_iter.txt");
        for (Uint i= 0; i != quaqua.num_outer_iter.size(); ++i)
            os << i << '\t' << quaqua.num_outer_iter[i] << '\n';
        os << '\n';
        for (Uint i= 0; i != quaqua.num_inner_iter.size(); ++i)
            os << i << '\t' << quaqua.num_inner_iter[i] << '\n';
    }
    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)//write to file
        DROPS::WriteFEToFile( xp3, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path") + "_p3", P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<> nobnd;
    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    if (vtkwriter.get() != 0)
    {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, xp3, "InterfaceSolP3"));
        vtkwriter->Register( make_VTKScalar(      make_P3Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Register( make_VTKVector( make_P3Eval( mg, nobnd_vec, lsgradrec), "LSGradRec") );
        vtkwriter->Register( make_VTKVector( make_P3Eval( mg, nobnd_vec, to_iface), "to_iface") );
        vtkwriter->Write( 0.);
    }
#endif
}



void StationaryStrategyDeformationP2 (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    // Initialize level set and triangulation
    adap.MakeInitialTriang();
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    // LinearLSInit( mg, lset.Phi, &the_lset_fun);
    LSInit( mg, lset.Phi, the_lset_fun, 0.);

    // Setup an interface-P2 numbering
    DROPS::IdxDescCL ifacep2idx( P2IF_FE);
    ifacep2idx.GetXidx().SetBound( P.get<double>("SurfTransp.XFEMReduced"));
    ifacep2idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "P2-NumUnknowns: " << ifacep2idx.NumUnknowns() << std::endl;

    // Compute the mesh deformation
    IdxDescCL vecp2idx( vecP2_FE);
    vecp2idx.CreateNumbering( mg.GetLastLevel(), mg);
    IdxDescCL p2idx( P2_FE);
    p2idx.CreateNumbering( mg.GetLastLevel(), mg);
    VecDescCL deformation( &vecp2idx);
    LocalQuaMapperCL locqua (mg, lset.Phi,
                             /*maxiter*/ P.get<int>( "LevelsetMapper.Iter"),
                             /*tol*/ P.get<double>( "LevelsetMapper.Tol"),
                             /*armijo_c*/ P.get<double>( "LevelsetMapper.ArmijoConstant"),
                             /*max_damping_steps*/ P.get<Uint>( "LevelsetMapper.MaxDampingSteps"));
    locqua.set_trust_region (P.get<double>( "LevelsetMapper.TrustRegion"))
    .set_deformation_method (P.get<std::string>( "LevelsetMapper.DeformationMethod") == "map_local_level_sets" ? LocalQuaMapperCL::MAP_LOCAL_LEVEL_SETS : LocalQuaMapperCL::MAP_ZERO_LEVEL_SETS);
    LocalQuaMapperDistanceP2CL locquap2(locqua); // Provides the interface for the Oswald-projection class.
    VecDescCL locdist_vd ( &p2idx);
    OswaldProjectionP2AccuCL<LocalQuaMapperDistanceP2CL> loc_dist_accu(locquap2, locdist_vd);
    loc_dist_accu.set_level_set_function (&lset.Phi, &lset.GetBndData(), &PrincipalLatticeCL::instance (1))
    .set_check_averaging (true);
    TetraAccumulatorTupleCL accus2;
    accus2.push_back( &loc_dist_accu);
    LocalQuaMapperDeformationP2CL locquadefp2(locqua); // Provides the interface for the Oswald-projection class.
    OswaldProjectionP2AccuCL<LocalQuaMapperDeformationP2CL> loc_def_accu(locquadefp2, deformation);
    loc_def_accu.set_level_set_function (&lset.Phi, &lset.GetBndData(), &PrincipalLatticeCL::instance (1))
    .set_check_averaging (true);
    accus2.push_back( &loc_def_accu);
    accumulate( accus2, mg, p2idx.TriangLevel(), p2idx.GetBndInfo());

    // Compute neighborhoods of the tetras at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);

    VecDescCL to_iface( &vecp2idx);
//     {
//         TetraAccumulatorTupleCL accus;
//         InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);
//         accus.push_back( &cdatap2);
//         InterfaceDebugP2CL p2debugaccu( cdatap2);
// //         p2debugaccu.store_offsets( to_iface);
//         p2debugaccu.set_true_area( 4.*M_PI*RadDrop[0]*RadDrop[0]);
//         p2debugaccu.set_ref_dp( &dp_sphere);
//         p2debugaccu.set_ref_abs_det( &abs_det_sphere);
//         accus.push_back( &p2debugaccu);
//         accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     }

    TetraAccumulatorTupleCL accus;
    InterfaceCommonDataDeformP2CL cdatap2( lset.Phi, lset.GetBndData(), deformation, lat);
    accus.push_back( &cdatap2);
    // Setup a P1 numbering
    DROPS::IdxDescCL p1idx( P1_FE);
    p1idx.CreateNumbering( mg.GetLastLevel(), mg);
    std::cout << "P1-NumUnknowns: " << p1idx.NumUnknowns() << std::endl;
    VecDescCL d_iface_vd (&p1idx);
    InterfaceApproxErrorDeformAccuCL ifaceerroraccu (&cdatap2, /*yg*/ 0, &d_iface_vd);
    ifaceerroraccu.set_d  (&sphere_dist)
    .set_Dd (&d_sphere_dist);
    accus.push_back( &ifaceerroraccu);

//     accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetMatchingFunction(), ifacep2idx.GetBndInfo());
//     exit (0);

    DROPS::MatDescCL Mp2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalMassDeformP2CL, InterfaceCommonDataDeformP2CL> accuMp2( &Mp2, LocalMassDeformP2CL(), cdatap2, "Mp2");
    accus.push_back( &accuMp2);
    DROPS::MatDescCL Ap2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalLaplaceBeltramiDeformP2CL, InterfaceCommonDataDeformP2CL> accuAp2( &Ap2, LocalLaplaceBeltramiDeformP2CL( P.get<double>("SurfTransp.Visc")), cdatap2, "Ap2");
    accus.push_back( &accuAp2);
    DROPS::MatDescCL Anp2( &ifacep2idx, &ifacep2idx);
    InterfaceMatrixAccuCL<LocalNormalLaplaceDeformP2CL, InterfaceCommonDataDeformP2CL> accuAnp2( &Anp2, LocalNormalLaplaceDeformP2CL (P.get<double>("SurfTransp.NormalLaplaceCoefficent")), cdatap2, "Anp2");
    accus.push_back( &accuAnp2);
    DROPS::VecDescCL bp2( &ifacep2idx);
    InterfaceVectorAccuCL<LocalVectorDeformP2CL, InterfaceCommonDataDeformP2CL> acculoadp2( &bp2, LocalVectorDeformP2CL( the_rhs_fun, bp2.t), cdatap2);
    accus.push_back( &acculoadp2);

    accumulate( accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

    double h = P.get<DROPS::Point3DCL>("Mesh.E1")[0]/P.get<double>("Mesh.N1")*std::pow(2., -P.get<double>("Mesh.AdaptRef.FinestLevel"));

    std::cout << "h: " << h  << '\n';

    DROPS::MatrixCL Lp2;
    Lp2.LinComb (1.0, Ap2.Data, 1.0, Mp2.Data, 1.0, Anp2.Data);
//     MatrixCL& Lp2= Mp2.Data;
//
//     DROPS::WriteToFile( Ap2.Data, "ap2_iface.txt", "Ap2");
//     DROPS::WriteToFile( Mp2.Data, "mp2_iface.txt", "Mp2");
//     DROPS::WriteToFile( Anp2.Data, "anp2_vol.txt", "Anp2");
//     DROPS::WriteFEToFile( bp2, mg, "rhsp2_iface.txt", /*binary=*/ false);

    typedef DROPS::SSORPcCL SurfPcT;
// //     typedef DROPS::JACPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Solver.Iter"), P.get<double>("SurfTransp.Solver.Tol"), true);

    DROPS::VecDescCL xp2( &ifacep2idx);
    surfsolver.Solve( Lp2, xp2.Data, bp2.Data, xp2.RowIdx->GetEx());
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    TetraAccumulatorTupleCL err_accus;
    err_accus.push_back( &cdatap2);
    InterfaceL2AccuDeformP2CL L2_accu( cdatap2, mg, "deformed P2-solution");
    L2_accu.set_grid_function( xp2);
    L2_accu.set_function( the_sol_fun, 0.);
    L2_accu.set_grad_function( the_sol_grad_fun, 0.);
    err_accus.push_back( &L2_accu);
    accumulate( err_accus, mg, ifacep2idx.TriangLevel(), ifacep2idx.GetBndInfo());

    if (P.get<int>( "SurfTransp.SolutionOutput.Freq") > 0)
        DROPS::WriteFEToFile( xp2, mg, P.get<std::string>( "SurfTransp.SolutionOutput.Path") + "_p2", P.get<bool>( "SurfTransp.SolutionOutput.Binary"));

    DROPS::NoBndDataCL<> nobnd;
    DROPS::NoBndDataCL<Point3DCL> nobnd_vec;
    VecDescCL the_sol_vd( &lset.idx);
    LSInit( mg, the_sol_vd, the_sol_fun, /*t*/ 0.);
    if (vtkwriter.get() != 0)
    {
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
        vtkwriter->Register( make_VTKIfaceScalar( mg, xp2, "InterfaceSolP2"));
        vtkwriter->Register( make_VTKScalar(      make_P2Eval( mg, nobnd, the_sol_vd),  "TrueSol"));
        vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, deformation), "deformation") );
//         vtkwriter->Register( make_VTKVector( make_P2Eval( mg, nobnd_vec, to_iface), "to_iface") );
        vtkwriter->Register( make_VTKScalar( make_P1Eval( mg, nobnd, d_iface_vd), "d_iface") );
        vtkwriter->Register( make_VTKScalar( make_P2Eval( mg, nobnd, locdist_vd), "locdist") );
        vtkwriter->Write( 0.);
    }
}


//void lsFun(double x, double y, double z, double *value)
//{
//    *value = x * x + y * y + z * z - 1.0;
//}
//
//
//void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
//{
//    grad[0] = x + x;
//    grad[1] = y + y;
//    grad[2] = z + z;
//}

//static void
//rhsIntFunP1(double x, double y, double z, double *f)
//{
//    *f = (x + y + z) << 1+1;
//}


#define _TEST
#ifdef TEST
void
ls(double x, double y, double z, double *value)
/* the level set function */
{
    *value = x * x + y * y + z * z - .5 * .5;
}


void
ls_grad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    grad[0] = x + x;
    grad[1] = y + y;
    grad[2] = z + z;
}

void
u(double x, double y, double z, double *f)
/* the integrand */
{
    *f = (x + y + z) * (x + y + z);
}
#endif

int main (int argc, char* argv[])
{
    //template for the use of high order quad module in phg

#ifdef TEST

    double ref;
    double res;
    int n;
    int order = 7,type = 0;
    double tet[4][3] =  {{1., 0., 0.},	/* vertex 0 */
        {0., 1., 0.},	/* vertex 1 */
        {0., 0., 1.},	/* vertex 2 */
        {0., 0., 0.}
    };	/* vertex 3 */
    if (type == 0)
    {
        ref = (0.2231747704246810387019576057274844651312);
    }
    else
    {
        ref = (0.02231747704246810387019576057274844651312);
        if (type > 0)
            ref = (0.1) - ref;
    }

//typedef enum {DOF_PROJ_NONE, DOF_PROJ_DOT, DOF_PROJ_CROSS} DOF_PROJ;


    n = phgQuadInterface2(
            ls,		/* the level set function */
            2,		/* polynomial order of the level set function */
            ls_grad,		/* the gradient of the level set function */
            tet,		/* coordinates of the vertices of the tetra */
            u,		/* the integrand */
            1,		/* dimension of the integrand */
            DOF_PROJ_NONE,	/* projection type for surface integral */
            1,		/* integration type (-1, 0, 1) */
            order,		/* order of the 1D Gaussian quadrature */
            &res,		/* the computed integral */
            NULL		/* pointer returning the computed rule */
        );
    std::cout<<"the domain quad result is: "<<res<<std::endl;

    n = phgQuadInterface2(
            ls,		/* the level set function */
            2,		/* polynomial order of the level set function */
            ls_grad,	/* the gradient of the level set function */
            tet,		/* coordinates of the vertices of the tetra */
            u,		/* the integrand */
            1,		/* dimension of the integrand */
            DOF_PROJ_NONE,	/* projection type for surface integral */
            type,		/* integration type (-1, 0, 1) */
            order,		/* order of the 1D Gaussian quadrature */
            &res,		/* the computed integral */
            NULL		/* pointer returning the computed rule */
        );
    std::cout<<"the surface quad result is: "<<res<<std::endl;

    std::cout<< "the relative error is: "<<fabs(res-ref)/ref<<std::endl;
    getchar();

#endif // TEST


    try
    {
        ScopeTimerCL timer( "main");

        DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/surfactant/surfactant/surfactantxyz.json");
        //P.read_json("../../param/surfactant/surfactant/surfactantxyz.json");
        std::cout << P << std::endl;

        DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

        std::cout << "Setting up interface-PDE.\n";
        WindVelocity= P.get<DROPS::Point3DCL>("SurfTransp.Exp.Velocity");
        RadDrop=      P.get<DROPS::Point3DCL>("SurfTransp.Exp.RadDrop");
        PosDrop=      P.get<DROPS::Point3DCL>("SurfTransp.Exp.PosDrop");
        RadTorus=     P.get<DROPS::Point2DCL>("SurfTransp.Exp.RadTorus");
        the_wind_fun= invecmap[P.get<std::string>("SurfTransp.Exp.Wind")];
        the_lset_fun= inscamap[P.get<std::string>("SurfTransp.Exp.Levelset")];
        the_normal_fun= invecmap[P.get<std::string>("SurfTransp.Exp.Normal")];
        the_rhs_fun=  inscamap[P.get<std::string>("SurfTransp.Exp.Rhs")];
        the_sol_fun=  inscamap[P.get<std::string>("SurfTransp.Exp.Solution")];
        the_sol_grad_fun = invecmap[P.get<std::string>("SurfTransp.Exp.SurfGradSol")];

        if (P.get<std::string>("SurfTransp.Exp.Solution") == "LaplaceBeltramixyzSol")
            the_sol_grad_fun=  &laplace_beltrami_xyz_sol_grad;
        for (Uint i= 0; i < 6; ++i)
            bf_wind[i]= the_wind_fun;

        std::cout << "Setting up domain:\n";
        std::unique_ptr<MGBuilderCL> builder( make_MGBuilder( P));
        DROPS::MultiGridCL mg( *builder);
        typedef DistMarkingStrategyCL MarkerT;
        MarkerT marker( the_lset_fun, P.get<double>( "Mesh.AdaptRef.Width"),
                        P.get<int>( "Mesh.AdaptRef.CoarsestLevel"), P.get<int>( "Mesh.AdaptRef.FinestLevel"));
        //mark band

        DROPS::AdapTriangCL adap( mg, &marker);

        // DROPS::LevelsetP2CL lset( mg, lsbnd, sf);
        DROPS::LevelsetP2CL& lset( *LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")) );

        if (P.get<int>("VTK.Freq",0))
            vtkwriter= std::unique_ptr<VTKOutCL>( new VTKOutCL(
                    adap.GetMG(),
                    "DROPS data",
                    P.get<int>("Time.NumSteps")/P.get<int>("VTK.Freq") + 1,
                    P.get<std::string>("VTK.VTKDir"),
                    P.get<std::string>("VTK.VTKName"),
                    P.get<std::string>("VTK.TimeFileName"),
                    P.get<int>("VTK.Binary"),
                    P.get<bool>("VTK.UseOnlyP1"),
                    false, /* <- P2DG */
                    -1,    /* <- level */
                    P.get<bool>("VTK.ReUseTimeFile")));
        if(true)//(P.get<bool>( "SurfTransp.Exp.StationaryPDE"))
        {
            if (P.get<std::string>("LevelsetMapper.DeformationMethod") != "")
                StationaryStrategyDeformationP2( mg, adap, lset);
            else
            {
                if (P.get<int>( "SurfTransp.FEDegree") == 1)
                    StationaryStrategyP1( mg, adap, lset);
                else
                    StationaryStrategyP2( mg, adap, lset);//p2 fem
            }
        }
        else
            Strategy( mg, adap, lset);

        delete &lset;
        rusage usage;
        getrusage( RUSAGE_SELF, &usage);
        std::cout << "ru_maxrss: " << usage.ru_maxrss << " kB.\n";
        return 0;
    }
    catch (DROPS::DROPSErrCL err)
    {
        err.handle();
    }
}
