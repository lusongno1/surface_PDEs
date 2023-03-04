/// \file
/// \brief Some settings about PDE paras, and some global variables, functions, and classes
/// \author LSEC: Song Lu

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
#include "sfpde.h"

using namespace DROPS;
//**********************************************set global variables*******************************************************/
double tet[4][3];
int iG;
int jG;
int orderG = 11;//10;
double gradTri[4][3];
int GLOBAL_TMP_COUNT = 0;

//***************************************************define test case*******************************************************/
static DROPS::RegisterScalarFunction regsca_level_set_function_drops( "LevelSetFunDrops", level_set_function_drops);
static DROPS::RegisterScalarFunction regsca_xyz_rhs( "xyzRhs", xyz_rhs);
static DROPS::RegisterScalarFunction regsca_laplace_beltrami_xyz_sol( "LaplaceBeltramixyzSol", laplace_beltrami_xyz_sol);

double zero_fun (const DROPS::Point3DCL& p, double)
{
    return 0.0;
}


void oneFun(double x, double y, double z, double *value)
{
    //*value = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
     *value = 1;
}


// test case 1
//define right hand side and true solution
//my test case，f = 3*(x+y+z) for problem -\Delta u + u = f
//then u = f/3
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return 3*(p[0]+p[1]+p[2]);//p.norm();
    //return 3*(p[0]+p[1]+p[2])/p.norm();
}
//my test case u=x+y+z
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p[0]+p[1]+p[2]);//p.norm();
    //return (p[0]+p[1]+p[2])/p.norm();
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{1,1,1};
    return tmp;
}
#endif


// test case 2, constant
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    return 1;
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return 1;
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{0,0,0};
    return tmp;
}
#endif


// test case 3, a general test case
//define right hand side and true solution
//u = a*|x|^2/(12+|x|^2)*(3x1^2x2-x2^3)
//f = a*(3x1^2x2-x2^3)
#if 1
double level_set_function_drops(const DROPS::Point3DCL& p, double)//directly modified in routine
{

    double x = p[0],y=p[1],z=p[2];
    return x * x + y * y + z * z - 1.0;
    //return p.norm()-1.0;
}


void lsFun(double x, double y, double z, double *value)
{
    *value = x * x + y * y + z * z - 1.0;
}


void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
{
    grad[0] = x + x;
    grad[1] = y + y;
    grad[2] = z + z;
}
double a(1.0);
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    return 3.0*p[0]*p[0]*p[1]-p[1]*p[1]*p[1];
    //return a/std::pow( p.norm(), 3.)*(3.*p[0]*p[0]*p[1]-p[1]*p[1]*p[1]);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p.norm_sq()/(12.+p.norm_sq()))*xyz_rhs(p,0.);
}
//DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
//{
////   DROPS::Point3DCL tmp{6*a/13*p[0]*p[1],-3*a/13*p[1]*p[1],0};
//    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)
//                          *( DROPS::MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) -
//                             (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
//    return tmp;// This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
//}


DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
     DROPS::Point3DCL tmp{6*a/13*p[0]*p[1],3*a/13*(p[0]*p[0]-p[1]*p[1]),0};
    //   DROPS::Point3DCL tmp{6*a/13*p[0]*p[1],-3*a/13*(p[0]*p[0]-p[1]),0};
//    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)
//                          *( DROPS::MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) -
//                             (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
    return tmp;// This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
}




#endif


// test case 3b, case 3 with a surface shift
#if 0
double shift = sqrt(2)-1+0.3;
double level_set_function_drops(const DROPS::Point3DCL& p, double)//directly modified in routine
{

    double x = p[0],y=p[1],z=p[2];

    x = x+shift;
    y = y+shift;
    z = z+shift;
    return x * x + y * y + z * z - 1.0;
    //return p.norm()-1.0;
}


void lsFun(double x, double y, double z, double *value)
{
    x = x+shift;
    y = y+shift;
    z = z+shift;
    *value = x * x + y * y + z * z - 1.0;
}


void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
{
    x = x+shift;
    y = y+shift;
    z = z+shift;
    grad[0] = x + x;
    grad[1] = y + y;
    grad[2] = z + z;
}
double a(1.0);
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    return 3.0*p[0]*p[0]*p[1]-p[1]*p[1]*p[1];
    //return a/std::pow( p.norm(), 3.)*(3.*p[0]*p[0]*p[1]-p[1]*p[1]*p[1]);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p.norm_sq()/(12.+p.norm_sq()))*xyz_rhs(p,0.);
}
DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
//   DROPS::Point3DCL tmp{6*a/13*p[0]*p[1],-3*a/13*p[1]*p[1],0};
    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)
                          *( DROPS::MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) -
                             (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
    return tmp;// This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
}

#endif


//test case 4
//define right hand side and true solution
//u = x*y*z
//f = x*y*z - 12*x*y*z*(x^2 + y^2 + z^2 - 2)
//#if 1
//double xyz_rhs (const DROPS::Point3DCL& p, double)
//{
//
//    return p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
//}
//double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
//{
//    return p[0]*p[1]*p[2];
//}
//
//
//DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
//{
//    DROPS::Point3DCL tmp{p[1]*p[2],p[0]*p[2],p[0]*p[1]};
//    return tmp;
//}
//
//double level_set_function_drops (const DROPS::Point3DCL& p, double)
//{
//    DROPS::Point3DCL RadDrop(1,1,1);
//    DROPS::Point3DCL PosDrop(0,0,0);
//    DROPS::Point3DCL x( p - PosDrop);
//    //double value=0;
//    //lsFun(x[0], x[1], x[2], &value);
//    return x.norm() - RadDrop[0];
//    //return value;
//}
//static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);
//
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
//#endif


//test case 5
//define level set function:atom
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    double x = p[0];
    double y = p[1];
    double z = p[2];

    //return p[0]+p[1]+p[2];//p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
    return pow(y,2)*sin(x)+exp(z);//p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[1]*p[2];
}


DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{p[1]*p[2],p[0]*p[2],p[0]*p[1]};
    return tmp;
}




double level_set_function_drops (const DROPS::Point3DCL& p, double)
{

    double x = p[0];
    double y = p[1];
    double z = p[2];
    //DROPS::Point3DCL RadDrop(1,1,1);
    //DROPS::Point3DCL PosDrop(0,0,0);
    //DROPS::Point3DCL x( p - PosDrop);
    //double value=0;
    //lsFun(x[0], x[1], x[2], &value);
   // double result = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
    double result = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
    return result;
    //return value;
}
static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);

void lsFun(double x, double y, double z, double *value)
{
    //*value = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
     *value = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
}


void lsGrad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    //grad[0] = x*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1);
    grad[0] = x*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1);
    //grad[1] = (y*((x*x)*2.75625E+5+(y*y)*7.77924E+5+(z*z)*2.75625E+5-1.00775E+6))/3.90625E+5;
    grad[1] = y*(-1.28E+2/2.5E+1)+y*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1)*2.8224;
    //grad[2] = z*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1);
    grad[2] = z*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1);
}
#endif


//test case 6
//define level set, Donut
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return p[0]+p[1]+p[2];//p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[1]*p[2];
}


DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{p[1]*p[2],p[0]*p[2],p[0]*p[1]};
    return tmp;
}




double level_set_function_drops (const DROPS::Point3DCL& p, double)
{

    double x = p[0];
    double y = p[1];
    double z = p[2];
    //DROPS::Point3DCL RadDrop(1,1,1);
    //DROPS::Point3DCL PosDrop(0,0,0);
    //DROPS::Point3DCL x( p - PosDrop);
    //double value=0;
    //lsFun(x[0], x[1], x[2], &value);
   // double result = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
    double result = sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)-3.0/5.0;;
    return result;
    //return value;
}
static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);

void lsFun(double x, double y, double z, double *value)
{
    //*value = pow((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1,2.0)-(y*y)*(6.4E+1/2.5E+1)-1.3E+1/1.0E+1;
     *value = sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)-3.0/5.0;
}


void lsGrad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    //grad[0] = x*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1);
    grad[0] = x*1.0/sqrt(x*x+y*y)*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)*(sqrt(x*x+y*y)-1.0);
    //grad[1] = (y*((x*x)*2.75625E+5+(y*y)*7.77924E+5+(z*z)*2.75625E+5-1.00775E+6))/3.90625E+5;
    grad[1] = y*1.0/sqrt(x*x+y*y)*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)*(sqrt(x*x+y*y)-1.0);
    //grad[2] = z*((x*x)/4.0+(y*y)*(4.41E+2/6.25E+2)+(z*z)/4.0+9.0/1.0E+1);
    grad[2] = z*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z);
}
#endif




//test case 6
//define level set function:tooth
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    return p[0]+p[1]+p[2];

    //return p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[1]*p[2];
}


DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{p[1]*p[2],p[0]*p[2],p[0]*p[1]};
    return tmp;
}




double level_set_function_drops (const DROPS::Point3DCL& p, double)
{
    double x = p[0];
    double y = p[1];
    double z = p[2];
    //DROPS::Point3DCL RadDrop(1,1,1);
    //DROPS::Point3DCL PosDrop(0,0,0);
    //DROPS::Point3DCL x( p - PosDrop);
    //double value=0;
    //lsFun(x[0], x[1], x[2], &value);
    double result = (x*x)*(-1.6E+1/2.5E+1)+(x*x*x*x)*(2.56E+2/6.25E+2)-(y*y)*(1.6E+1/2.5E+1)+(y*y*y*y)*(2.56E+2/6.25E+2)-(z*z)*(1.6E+1/2.5E+1)+(z*z*z*z)*(2.56E+2/6.25E+2);
    //return value;
}
static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);

void lsFun(double x, double y, double z, double *value)
{
    *value = (x*x)*(-1.6E+1/2.5E+1)+(x*x*x*x)*(2.56E+2/6.25E+2)-(y*y)*(1.6E+1/2.5E+1)+(y*y*y*y)*(2.56E+2/6.25E+2)-(z*z)*(1.6E+1/2.5E+1)+(z*z*z*z)*(2.56E+2/6.25E+2);
}

void lsGrad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    grad[0] = x*((x*x)*3.2E+1-2.5E+1)*(3.2E+1/6.25E+2);
    grad[1] = y*((y*y)*3.2E+1-2.5E+1)*(3.2E+1/6.25E+2);
    grad[2] = z*((z*z)*3.2E+1-2.5E+1)*(3.2E+1/6.25E+2);
}
#endif

//test for the gyroid, heart-shape and torus, justest have a thy
#if 0
double level_set_function_drops(const DROPS::Point3DCL& p, double)//directly modified in routine
{
    //std::cout<<M_PI<<std::endl;
    //double phi = cos(M_PI*p[0])*sin(M_PI*p[1])+
    //cos(M_PI*p[1])*sin(M_PI*p[2])+cos(M_PI*p[2])*sin(M_PI*p[0]);
    double x = p[0],y=p[1],z=p[2];
    //double phi = 2*pow(x-0.5,2)-8*pow(y-0.5,3)-16*pow(z-0.5,4)-1/50;
    //double phi = pow((pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1),3)-pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3);
    double phi = sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)-3.0/5.0;;
    return phi;
//    return std::sqrt( std::pow( RadTorus[0] - std::sqrt(p[0]*p[0] + p[1]*p[1]), 2) + std::pow( p[2], 2)) - RadTorus[1];
}
static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);

//void lsFun(double x, double y, double z, double *value)
//{
//    *value = x * x + y * y + z * z - 1.0;
//}


//void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
//{
//    grad[0] = x + x;
//    grad[1] = y + y;
//    grad[2] = z + z;
//}

void lsFun(double x, double y, double z, double *value)
{
    double R = 1;
    double r = 0.6;
    //double phi = pow((pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1),3)-pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3);
    double phi = sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)-3.0/5.0;;
    *value = phi;//2*pow(x-0.5,2)-8*pow(y-0.5,3)-16*pow(z-0.5,4)-1/50;
    //cos(M_PI*x)*sin(M_PI*y)+cos(M_PI*y)*sin(M_PI*z)+cos(M_PI*z)*sin(M_PI*x);
}
//
void lsGrad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    grad[0] = x*1.0/sqrt(x*x+y*y)*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)*(sqrt(x*x+y*y)-1.0);
    grad[1] = y*1.0/sqrt(x*x+y*y)*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)*(sqrt(x*x+y*y)-1.0);
    grad[2] = z*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z);

    //grad[0] = 6*x*pow(pow(x,2) + (9*pow(y,z))/4 + pow(z,2) - 1,2) - 2*x*pow(z,3);
    //grad[1] =      (27*y*pow(pow(x,2) + (9*pow(y,z))/4 + pow(z,2) - 1,2))/2 - (9*y*pow(z,3))/40;
    //grad[2] = 6*z*pow(pow(x,2) + (9*pow(y,z))/4 + pow(z,2) - 1,2) - (27*pow(y,z)*pow(z,2))/80 - 3*pow(x,2)*pow(z,2);
    //grad[0] =   4*x - 2;
    //grad[1]  = -24*pow(y - 1/2,2);
    //grad[2] = -64*pow(z - 1/2,3);
    // grad[0] = M_PI*cos(M_PI*x)*cos(M_PI*z) - M_PI*sin(M_PI*x)*sin(M_PI*y);
    // grad[1] = M_PI*cos(M_PI*x)*cos(M_PI*y) - M_PI*sin(M_PI*y)*sin(M_PI*z);
    // grad[2] = M_PI*cos(M_PI*y)*cos(M_PI*z) - M_PI*sin(M_PI*x)*sin(M_PI*z);
}

double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    double x = p[0],y=p[1],z=p[2];
    //return pow(x,2)*sin(y)*exp(z);

    //return p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
    //double tmp = x+y+z;
    //return std::sqrt( p[2]*p[2] + std::pow( std::sqrt( p[0]*p[0] + p[1]*p[1]) - 1.0, 2));
    return x*sin(y)*exp(z);


    //return tmp;
}

double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return 1.0;
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{1.0,1.0,1.0};
    return tmp;
}
#endif



//****************************************************some useful funtions ***************************************************************/

void vecMinus(double a[3],double b[3],double (&result)[3])
{
    for(int i=0; i<3; i++)
    {
        result[i] = a[i]-b[i];
    }
}

void crossMul(double a[3],double b[3],double (&p)[3])
{
    p[0] = a[1]*b[2] - a[2]*b[1];
    p[1] = a[2]*b[0] - a[0]*b[2];
    p[2] = a[0]*b[1] - a[1]*b[0];
}
//template <typename T>
//double dotP3(T a,T b)
//{
//
//    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
//}


double getBaryCoord(double tetra[4][3],int i,double x,double y,double z)
{
    double pValue = 0;
    double v0[3] = {tetra[i][0],tetra[i][1],tetra[i][2]};
    int idx = 0;
    double vGround[3][3];
    for(int j=0; j<4; j++)
    {
        if(j==i)
            continue;
        for(int k=0; k<3; k++)
            vGround[idx][k] = tetra[j][k];
        idx++;
    }
    double vec1[3];
    double vec2[3];
    double n[3];
    vecMinus(vGround[1],vGround[0],vec1);
    vecMinus(vGround[2],vGround[0],vec2);
    crossMul(vec1,vec2,n);
    double n_norm = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    for(int j=0; j<3; j++)
        n[j] /=n_norm;
    double vecV[3] = {v0[0] - vGround[0][0],v0[1] - vGround[0][1],v0[2] - vGround[0][2]};
    double vecX[3] = {x - vGround[0][0],y - vGround[0][1],z - vGround[0][2]};
    double valXYZ = dotP3(vecX,n);
    double valV = dotP3(vecV,n);
    //assert((valV>=0&&valXYZ>=0)||(valV<=0&&valXYZ<=0));
    pValue = valXYZ/valV;
    //assert()
    return pValue;
}

DROPS::BaryCoordCL getBaryCoords(double tetra[4][3],double x,double y,double z)
{
    double BaryCoordArr[4];
    for(int i=0; i<4; i++)
    {
        BaryCoordArr[i] = getBaryCoord(tet,i, x,y,z);
    }
    const DROPS::BaryCoordCL& BaryCoord{BaryCoordArr[0],BaryCoordArr[1],BaryCoordArr[2],BaryCoordArr[3]};
    return BaryCoord;
}


void GetTet2DArr(const DROPS::TetraCL& t,double tet[4][3])
{
    for (int i= 0; i < 4; ++i)
    {
        auto vtx = t.GetVertex(i);
        auto coord = vtx->GetCoord();
        for(int j=0; j<3; j++)
        {

            tet[i][j] = coord[j];

        }
    }

}


void getSfNormalVec(double x,double y,double z,double (&n)[3])
{
    double ls_grad[3];
    lsGrad(x,y,z,ls_grad);
    double ls_grad_norm = std::sqrt(ls_grad[0]*ls_grad[0]+ls_grad[1]*ls_grad[1]+ls_grad[2]*ls_grad[2]);
    for(int i=0; i<3; i++)
    {
        n[i] = ls_grad[i]/ls_grad_norm;
    }

}

//template <typename T>
void getSurfaceGradient(DROPS::Point3DCL v,double n[3],double (&sf_grad)[3])
{
    double proj_norm = 0;
    for(int i=0; i<3; i++)
        proj_norm += v[i]*n[i];
    for(int i=0; i<3; i++)
    {
        sf_grad[i] = v[i] - proj_norm*n[i];
    }
}

void getSurfaceGradient(DROPS::Point3DCL v,double n[3],DROPS::SVectorCL<3> &sf_grad)
{
    double proj_norm = 0;
    for(int i=0; i<3; i++)
        proj_norm += v[i]*n[i];
    for(int i=0; i<3; i++)
    {
        sf_grad[i] = v[i] - proj_norm*n[i];
    }
}


void ouput_valarray(std::valarray<double> v)
{
    std::cout<<"begin output valarray:"<<std::endl;
    for(int i=0; i<v.size(); i++)
    {
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}
void cout2txt(double a)
{
    std::ofstream mycout("./debug.txt",std::ios_base::app);
    mycout<<a<<std::endl;
    mycout.close();
}

int nc = 1;
void coutTet(const DROPS::TetraCL& t)
{
    std::cout<<std::endl<<nc++<<":"<<std::endl;;
    for (int i= 0; i < 4; ++i)
    {
        auto vtx = t.GetVertex(i);
        auto coord = vtx->GetCoord();
        for(int j=0; j<3; j++)
        {

            auto tmp = coord[j];
            std::cout<<tmp<<" ";

        }
        std::cout<<std::endl;
    }

}

void getEigVec(DROPS::VecDescCL &eigVec,int &row)
{
    // 读文件
    std::ifstream inFile("eig_vec.txt", std::ios::in);
    std::string lineStr;
    std::vector<std::vector<double>> strMatrix;
    while (getline(inFile, lineStr))
    {
        // 打印整行字符串
        //cout << lineStr << endl;
        // 存成二维表结构
        std::stringstream ss(lineStr);
        std::string str;
        std::vector<double> lineArray;
        // 按照逗号分隔
        while (getline(ss, str, ','))
            lineArray.push_back(stof(str));
        strMatrix.push_back(lineArray);
    }

    //int row = 1;
    for(int i=0;i<strMatrix[row].size();i++)
    {
        eigVec.Data[i] = strMatrix[row][i];

    }

}





