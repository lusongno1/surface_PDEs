/// \file
/// \brief Including file: Some setting about PDE para, and some global function
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
//#pragma once
#include "misc/container.h"
#include "geom/simplex.h"
#include "num/discretize.h"
#include <fstream>
#include <math.h>
#include "geom/multigrid.h"
#include "num/lattice-eval.h"
#include <cmath>
#include "misc/funcmap.h"
#include <fstream>
#include <sstream>
#include <iostream>
//#include <vector.h>
#define BOOST_NO_EXCEPTIONS
#define BOOST_EXCEPTION_DISABLE

#ifndef DROPS_SFPDE_H
#define DROPS_SFPDE_H

using namespace DROPS;
//namespace DROPS
//{
double xyz_rhs (const DROPS::Point3DCL& p, double);
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double);
void lsFun(double x, double y, double z, double *value);
void lsGrad(double x, double y, double z, double *grad);
void oneFun(double x, double y, double z, double *value);
//}

void vecMinus(double a[3],double b[3],double (&result)[3]);
void crossMul(double a[3],double b[3],double (&p)[3]);
template <typename T>
double dotP3(T a,T b)
{

    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
double getBaryCoord(double tetra[4][3],int i,double x,double y,double z);
DROPS::BaryCoordCL getBaryCoords(double tetra[4][3],double x,double y,double z);
void GetTet2DArr(const DROPS::TetraCL& t,double tet[4][3]);
void getSfNormalVec(double x,double y,double z,double (&n)[3]);
//template <typename T>
void getSurfaceGradient(DROPS::Point3DCL v,double n[3],double (&sf_grad)[3]);
void getSurfaceGradient(DROPS::Point3DCL v,double n[3],DROPS::SVectorCL<3> &sf_grad);
DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double);
void ouput_valarray(std::valarray<double> v);
void cout2txt(double a);
void coutTet(const DROPS::TetraCL& t);
double level_set_function_drops (const DROPS::Point3DCL& p, double);
void getEigVec(DROPS::VecDescCL &eigVec,int &row);
double zero_fun (const DROPS::Point3DCL& p, double);

extern double tet[4][3];
extern int iG;
extern int jG;
extern int orderG;
extern double gradTri[4][3];//store gradients of shape functions
extern int GLOBAL_TMP_COUNT;
//extern DROPS::LocalP2CL<> localP2Set[10];



#endif
