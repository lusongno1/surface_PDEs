#ifndef DROPS_FEMP3_H
#define DROPS_FEMP3_H


#include "num/accumulator.h"
#include "geom/principallattice.h"
//#include "num/lattice-eval.h"
#include "geom/subtriangulation.h"
#include "geom/multigrid.h"
#include "num/quadrature.h"
#include "surfactant/ifacetransp.h"




//#include "surfactant/ifacetransp.h"
//#include "misc/params.h"
//#include "geom/builder.h"
//#include "levelset/levelset.h"
//#include "levelset/adaptriang.h"
//#include "levelset/surfacetension.h"
//#include "levelset/levelsetmapper.h"
//#include "out/output.h"
//#include "out/vtkOut.h"
//#include "misc/dynamicload.h"
//#include "misc/funcmap.h"
//#include "misc/omp_variable.h"
//#include "geom/subtriangulation.h"
//#include "num/gradient_recovery.h"
//#include "surfactant/sfpde.h"

using namespace DROPS;

template<class Data, class _BndData, class _VD>
class P3EvalCL;


void averaging_P3_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad);

typedef std::pair<const TetraCL*, BaryCoordCL> TetraBaryPairT;
typedef std::vector<TetraBaryPairT>            TetraBaryPairVectorT;

//template <class T, class ResultIterT>
//inline ResultIterT
//evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultIterT result_iterator);
//
//template <class T, class ResultContT>
//inline ResultContT&
//resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultContT& result_container);
//
//template <class PEvalT, class ResultIterT>
//inline ResultIterT
//evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultIterT result_iterator);
//
//template <class PEvalT, class ResultContT>
//inline ResultContT&
//resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultContT& result_container);

//template<class BndData_, class VD_>
//  P3EvalCL<typename BndData_::bnd_type, BndData_, VD_>
//    make_P3Eval (const MultiGridCL& mg, BndData_& bnd, VD_& vd);

    // Create a P3EvalCL without the agonizing template-pain.
template<class BndData_, class VD_>
  P3EvalCL<typename BndData_::bnd_type, BndData_, VD_>
    make_P3Eval (const MultiGridCL& mg, BndData_& bnd, VD_& vd);


template <class P3FuncT, class Cont>
void RestrictP3(const TetraCL& s, const P3FuncT& f, Cont& c);

template <class VecDescT, class BndDataT, class Cont>
void RestrictP3(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c);

//    template <class T, class ResultContT>
//  inline ResultContT&
//  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultContT& result_container);
//template <class PEvalT, class ResultIterT>
//  inline ResultIterT
//  evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultIterT result_iterator);
//template <class PEvalT, class ResultContT>
//  inline ResultContT&
//  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultContT& result_container);
//****************************************************some useful classed ***************************************************************/
//DROPS::LocalP2CL<> localP2Set[10];
//std::vector<double,4> v4type;
using v3 = DROPS::SVectorCL<3>;
using v4 = DROPS::SVectorCL<4>;
using v43 = DROPS::SMatrixCL<4,3>;
class TetrahedronFECL
{
private:
    v43 coordinates;//(4,std::vector<double>(3));
    v4 baryCoordTmp;
public:

    TetrahedronFECL(v43 coordinates):coordinates(coordinates) {};
    void getBaryCoord(double x,double y,double z)
    {
        for(DROPS::Uint i=0; i<4; i++)
        {
            double pValue = 0;
            double v0[3] = {coordinates(i,0),coordinates(i,1),coordinates(i,2)};
            int idx = 0;
            double vGround[3][3];
            for(int j=0; j<4; j++)
            {
                if(j==i)
                    continue;
                for(int k=0; k<3; k++)
                    vGround[idx][k] = coordinates(j,k);
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
            baryCoordTmp[i] = valXYZ/valV;
            //assert()
        }
    };
};

class TetrahedronP3FECL:public TetrahedronFECL
{
};


//**************************************************************************
// Class:   FE_P3CL                                                        *
// Purpose: Shape functions and their gradients for piecewise quadratic,   *
//          continuous finite elements on the reference tetrahedron        *
//          The number of the H-functions refers to the number of the      *
//          (mid-) vertex in the tetrahedron as defined in topo.h, where   *
//          the degree of freedom is located.                              *
//**************************************************************************
class FE_P3CL
{
  private:
    static const double _D2H[10][3][3];

  public:
    // default ctor, copy-ctor, assignment-op, dtor

    static const Uint NumDoFC= 10;

    // restriction of the shape functions to reference edge
    static double H0(double v1) { return 1. +v1*(2.*v1 -3.); }
    static double H1(double v1) { return v1*(2.*v1 -1.); }
    static double H2(double v1) { return 4.*v1*(1. -v1); }

    // restriction of the shape functions to reference face
    static double H0(double v1, double v2) { const double sum= v1 + v2; return 1. +sum*(2.*sum -3.); }
    static double H1(double v1, double)    { return v1*(2.*v1 -1.); }
    static double H2(double, double v2)    { return v2*(2.*v2 -1.); }
    static double H3(double v1, double v2) { return 4.*v1*( 1. -(v1 + v2) ); }
    static double H4(double v1, double v2) { return 4.*v2*( 1. -(v1 + v2) ); }
    static double H5(double v1, double v2) { return 4.*v1*v2; }

    // restriction of the shape functions to reference tetrahedron
    static double H0(double v1, double v2, double v3) { const double sum= v1 + v2 + v3; return 1. +sum*(2.*sum -3.); }
    static double H1(double v1, double, double)       { return v1*(2.*v1 -1.); }
    static double H2(double, double v2, double)       { return v2*(2.*v2 -1.); }
    static double H3(double, double, double v3)       { return v3*(2.*v3 -1.); }
    static double H4(double v1, double v2, double v3) { return 4.*v1*( 1. -(v1 + v2 + v3) ); }
    static double H5(double v1, double v2, double v3) { return 4.*v2*( 1. -(v1 + v2 + v3) ); }
    static double H6(double v1, double v2, double)    { return 4.*v1*v2; }
    static double H7(double v1, double v2, double v3) { return 4.*v3*( 1. -(v1 + v2 + v3) ); }
    static double H8(double v1, double, double v3)    { return 4.*v1*v3; }
    static double H9(double, double v2, double v3)    { return 4.*v2*v3; }
    static inline double H (Uint dof, double v1, double v2, double v3);

    // restriction of the shape functions to reference tetrahedron, barycentric coordinates
    static inline double H0(const BaryCoordCL& p) { return p[0]*(2.*p[0] - 1.); }
    static inline double H1(const BaryCoordCL& p) { return p[1]*(2.*p[1] - 1.); }
    static inline double H2(const BaryCoordCL& p) { return p[2]*(2.*p[2] - 1.); }
    static inline double H3(const BaryCoordCL& p) { return p[3]*(2.*p[3] - 1.); }
    static inline double H4(const BaryCoordCL& p) { return 4.*p[0]*p[1]; }
    static inline double H5(const BaryCoordCL& p) { return 4.*p[0]*p[2]; }
    static inline double H6(const BaryCoordCL& p) { return 4.*p[1]*p[2]; }
    static inline double H7(const BaryCoordCL& p) { return 4.*p[0]*p[3]; }
    static inline double H8(const BaryCoordCL& p)    { return 4.*p[1]*p[3]; }
    static inline double H9(const BaryCoordCL& p)    { return 4.*p[2]*p[3]; }
    static inline double H (Uint dof, const BaryCoordCL& p);

    template <class Cont>
      static inline typename ValueHelperCL<Cont>::value_type
      val(const Cont& c, const BaryCoordCL& p) {
          return LinearCombinationCL<Cont, typename ValueHelperCL<Cont>::value_type>::do_it( c, H0( p), H1( p), H2( p), H3( p), H4( p),  H5( p), H6( p), H7( p), H8( p), H9( p));
      }
    template <class Cont>
      static inline typename ValueHelperCL<Cont>::value_type
      val(const Cont& c, double v1, double v2, double v3) {
          return c[0] * H0( v1, v2, v3) + c[1] * H1( v1, v2, v3) + c[2] * H2( v1, v2, v3) + c[3] * H3( v1, v2, v3)
               + c[4] * H4( v1, v2, v3) + c[5] * H5( v1, v2, v3) + c[6] * H6( v1, v2, v3) + c[7] * H7( v1, v2, v3)
               + c[8] * H8( v1, v2, v3) + c[9] * H9( v1, v2, v3);
      }

    // pt[0]...pt[numpt-1] are coordinates where the shape-functions are evaluated.
    // v is an array of 10 valarrays. They are resized to have numpt components.
    // v[i] contains H_i( pt[0])...H_i( pt[numpt-1])
    static void ApplyAll(Uint numpt, const BaryCoordCL* const pt, std::valarray<double>* v);

    // gradients of the shape functions on the reference tetrahedron.
    // To obtain the gradient on tetra T: See comments in FE_P1CL.
    static inline SVectorCL<3> DH0Ref(double, double, double);
    static inline SVectorCL<3> DH1Ref(double, double, double);
    static inline SVectorCL<3> DH2Ref(double, double, double);
    static inline SVectorCL<3> DH3Ref(double, double, double);
    static inline SVectorCL<3> DH4Ref(double, double, double);
    static inline SVectorCL<3> DH5Ref(double, double, double);
    static inline SVectorCL<3> DH6Ref(double, double, double);
    static inline SVectorCL<3> DH7Ref(double, double, double);
    static inline SVectorCL<3> DH8Ref(double, double, double);
    static inline SVectorCL<3> DH9Ref(double, double, double);


    // DHREef(i, ...) == DHiRef(...)
    static inline SVectorCL<3> DHRef(Uint dof, double v1, double v2, double v3);

    // D2HRef(i) == second derivative of H_i; this is constant
    // diff( H_d(x), k,i) = sum( M_ij*M_kl*D2HRef(d, j, l), j,l=0..2)
    static inline double D2HRef(Uint dof, Uint r, Uint s)
        { return _D2H[dof][r][s]; }

    // Laplace(H_d) on a tetrahedron T; M:= transpose(inverse(A)), where A is the matrix
    // of the affine transformation that maps the reference tetrahedron onto T.
    static inline double Laplace(Uint dof, const SMatrixCL<3,3>& M);

    // The barycentric coordinates of the dofs.
    static const BaryCoordCL bary_coord[10];
};

template<class T= double>
class LocalP3CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;
    typedef FE_P3CL FETYPE;//FE_P2CL

  protected:
    typedef LocalP3CL<T> self_;

  public:
    LocalP3CL() : base_type( value_type(), FE_P3CL::NumDoFC) {}
    LocalP3CL(const value_type& t): base_type( t, FE_P3CL::NumDoFC) {}
    // Initialize from a given function
    LocalP3CL(const TetraCL&, instat_fun_ptr , double= 0.0);
    // Initialize from VecDescCL and boundary-data
    template<class BndDataT>
      LocalP3CL(const TetraCL&, const VecDescCL&, const BndDataT&);
    // Initialize from PiEvalCL
    template <class P3FunT>
      LocalP3CL(const TetraCL&, const P3FunT&);
    // Initialize from LocalP1CL
    LocalP3CL(const LocalP1CL<T>&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(LocalP3CL, T, base_type)

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const TetraCL&, instat_fun_ptr, double= 0.0);
    template<class BndDataT>
      inline self_&
      assign(const TetraCL&, const VecDescCL&, const BndDataT&);
    template<class BndDataT>
      inline self_&
      assign_on_tetra(const TetraCL&, const VecDescCL&, const BndDataT&);
    template <class P3FunT>
      inline self_&
      assign(const TetraCL&, const P3FunT&);
    inline self_&
    assign(const LocalP1CL<T>&);

    // pointwise evaluation in barycentric coordinates
    inline value_type operator()(const BaryCoordCL&) const;
};
//**************************************************************************
// Class:   LocalP3CL                                                      *
//**************************************************************************
template<class T>
  inline LocalP3CL<T>&
  LocalP3CL<T>::assign(const TetraCL& s, instat_fun_ptr f, double t)
{
    for (Uint i= 0; i< NumVertsC; ++i)
        (*this)[i]= f( s.GetVertex( i)->GetCoord(), t);
    for (Uint i= 0; i< NumEdgesC; ++i)
        (*this)[i+NumVertsC]= f( GetBaryCenter( *s.GetEdge( i)), t);
    return *this;
}

template<class T>
  template<class BndDataT>
    inline LocalP3CL<T>&
    LocalP3CL<T>::assign(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecT& v= vd.Data;
    if (vd.RowIdx->IsDG())
    { // This is just for P3
        Uint idx_num = vd.RowIdx->GetIdx();
        Uint first = s.Unknowns(idx_num);
        for (int i = 0; i < 10; ++i)
        {
            (*this)[i] = DoFT::get( vd.Data, first++);
        }
        return *this;
    }
    const Uint tlvl= s.GetLevel();
    const Uint vlvl= vd.GetLevel();
    const Uint idx= vd.RowIdx->GetIdx();
    if (tlvl == vlvl) {
        for (Uint i= 0; i< NumVertsC; ++i)
            (*this)[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
        for (Uint i= 0; i< NumEdgesC; ++i)
            (*this)[i+NumVertsC]= !bnd.IsOnDirBnd( *s.GetEdge( i))
                ? DoFT::get( v, s.GetEdge( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetEdge( i), vd.t);
    }
    else {
        if (tlvl < vlvl) RestrictP3( s, vd, bnd, *this);
        else throw DROPSErrCL( "LocalP3CL::Assign: Prolongation not implemented.\n");
    }
    return *this;
}

template<class T>
  template<class BndDataT>
    inline LocalP3CL<T>&
    LocalP3CL<T>::assign_on_tetra(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();

    //const Uint tlvl= s.GetLevel();
    //const Uint flvl= vd.GetLevel();

    if (vd.RowIdx->IsDG())
    { // This is just for P3
        Uint first = s.Unknowns(idx);
        for (int i = 0; i < 10; ++i)
        {
            (*this)[i] = DoFT::get( v, first++);
        }
    }
    else
    {
        for (Uint i= 0; i< NumVertsC; ++i)
            (*this)[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
        for (Uint i= 0; i< NumEdgesC; ++i)
            (*this)[i+NumVertsC]= !bnd.IsOnDirBnd( *s.GetEdge( i))
                ? DoFT::get( v, s.GetEdge( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetEdge( i), vd.t);
    }
    return *this;
}

template<class T>
  template<class P3FunT>
    inline LocalP3CL<T>&
    LocalP3CL<T>::assign(const TetraCL& s, const P3FunT& f)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecDescCL& vd = *(f.GetSolution());
    const Uint tlvl= s.GetLevel();
    const Uint flvl= vd.GetLevel();
    //const Uint idx= vd.RowIdx->GetIdx();

    if (vd.RowIdx->IsDG())
    { // This is just for P3
        if (tlvl != flvl)
            throw DROPSErrCL( "LocalP3CL::Assign: Prolongation not implemented.\n");
        Uint idx_num = vd.RowIdx->GetIdx();
        Uint first = s.Unknowns(idx_num);
        for (int i = 0; i < 10; ++i)
        {
            (*this)[i] = DoFT::get( vd.Data, first++);
        }
        return *this;
    }
    else

    {
        const Uint tlvl= s.GetLevel();
        const Uint flvl= f.GetLevel();
        if (tlvl == flvl)
            f.GetDoF( s, *this);
        else
            if (tlvl < flvl) RestrictP3( s, f, *this);
            else throw DROPSErrCL( "LocalP3CL::Assign: Prolongation not implemented.\n");
    }
    return *this;
}

template<class T>
  inline LocalP3CL<T>&
  LocalP3CL<T>::assign(const LocalP1CL<T>& p1)
{
    for (size_t i= 0; i < 4; ++i)
        (*this)[i]= p1[i];
    for (size_t i= 0; i < 6; ++i)
        (*this)[i + 4]= 0.5*(p1[VertOfEdge( i, 0)] + p1[VertOfEdge( i, 1)]);
    return *this;
}


template<class T>
  LocalP3CL<T>::LocalP3CL(const TetraCL& s, instat_fun_ptr f , double t)
  : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( s, f, t);
}

template<class T>
  template <class P3FunT>
    LocalP3CL<T>::LocalP3CL(const TetraCL& s, const P3FunT& f)
    : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( s, f);
}

template<class T>
  template<class BndDataT>
    LocalP3CL<T>::LocalP3CL(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd)
    : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( s, vd, bnd);
}

template<class T>
  LocalP3CL<T>::LocalP3CL (const LocalP1CL<T>& p1)
    : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( p1);
}

template<class T>
  inline typename LocalP3CL<T>::value_type
  LocalP3CL<T>::operator() (const BaryCoordCL& p) const
{
    return FE_P3CL::val( *this, p);
}


class P3DiscCL
{
  public:
    // gradients on reference tetra
    static void GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[10]);
    static void GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10]);
    static void GetGradientsOnRef( Quad5CL<Point3DCL> GRef[10]);
    // The 2nd arg points to 3 vertices of the triangle
    static void GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[10], const BaryCoordCL* const);
    // p3[i] contains a LocalP3CL-object that is initialized with FE_P3CL::Hi
    static void GetP3Basis( LocalP3CL<> p3[10]);
    // p3[i] contains a Quad5_2DCL-object that is initialized with FE_P3CL::Hi
    static void GetP3Basis( Quad5_2DCL<> p3[10], const BaryCoordCL* const p);
    // compute gradients
    static void GetGradients( LocalP1CL<Point3DCL> G[10], const LocalP1CL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<4; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad2CL<Point3DCL> G[10], Quad2CL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<5; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad5CL<Point3DCL> G[10], Quad5CL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<Quad5DataCL::NumNodesC; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad5_2DCL<Point3DCL> G[10], Quad5_2DCL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<Quad5_2DDataCL::NumNodesC; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradient( Quad2CL<Point3DCL> &G, Quad2CL<Point3DCL> &GRef, const SMatrixCL<3,3> &T)
    { for (int j=0; j<5; ++j) G[j]= T*GRef[j]; }
    static void GetGradient( Quad5CL<Point3DCL> &G, Quad5CL<Point3DCL> &GRef, const SMatrixCL<3,3> &T)
    { for (int j=0; j<Quad5DataCL::NumNodesC; ++j) G[j]= T*GRef[j]; }
    /// compute gradient of a function provided as LocalP3CL<double> object
    template<class GradT>
    static void GetFuncGradient( GradT& gradF, const LocalP3CL<>& F, const GradT G[10])
    { gradF= F[0]*G[0]; for (int i=1; i<10; ++i) gradF+= F[i]*G[i]; }
    // Compute the Hessians H[d]= M*Href[d]*M^T
    static void GetHessians (SMatrixCL<3,3> H[10], const SMatrixCL<3,3>& M) {
        for (Uint d= 0; d< 10; ++d) {
            std::memset( &H[d], 0, 3*3*sizeof( double));
            for (Uint i= 0; i < 3; ++i)
                for (Uint j= 0; j < 3; ++j)
                    for (Uint k= 0; k < 3; ++k)
                        for (Uint l= 0; l < 3; ++l)
                            H[d](i,j)+= M(i, k)*M(j, l)*FE_P3CL::D2HRef( d, k, l);
        }
    }
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 1
    static inline SVectorCL<3> Quad( const TetraCL& tetra, instat_vector_fun_ptr, Uint, double= 0.0);
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 2
    template<class valT>
    static inline valT Quad( valT f[10], int i);
    // returns int phi_i phi_j dx
    static inline double GetMass( int i, int j);
    // returns int phi_i dx
    static inline double GetLumpedMass( int i) { return i<4 ? -1./110. : 1./30.; }
};


class LocalP3GradientCL
{
  private:
    const VecDescCL& f_;
    const BndDataCL<>& fbnd_;

    SMatrixCL<3,3> M;
    LocalP1CL<Point3DCL> GradRefLP1[10],
                         GradLP1[10],
                         p1grad;
    LocalP3CL<> p3;//LocalP2CL
    LocalP3CL<Point3DCL> p3grad;

  public:
    typedef Point3DCL value_type;
    static const int num_components= 3;

    LocalP3GradientCL (const VecDescCL& f, const BndDataCL<>& fbnd)
        : f_( f), fbnd_( fbnd)  { P3DiscCL::GetGradientsOnRef( GradRefLP1); }//P2DiscCL

    void set_tetra (const TetraCL* t);
    value_type&       operator[] (size_t i)       { return p3grad[i]; }
    const value_type& operator[] (size_t i) const { return p3grad[i]; }
    bool invalid_p (size_t /*i*/) const { return false; }
    void finalize_accumulation () const {}
};

/// \brief Collect indices of unknowns, boundary-segments and boundary
///     conditions on a tetrahedron.
///
/// This is convenient for discretisation of operators in the Setup-routines.
class LocalNumbP3CL
{
  public:
    /// \brief Field of unknown-indices; NoIdx, iff the degree of freedom lies
    /// on a boundary without unknowns. (Formerly called Numb.)
    IdxT     num   [10];
    /// \brief On boundaries, the number of the relevant BndSegDataCL-object
    /// in the corresponding BndDataCL-object, else NoBndC.
    BndIdxT  bndnum[10];
    /// \brief The relevant BndCondT, NoBC in the interior dofs.
    BndCondT bc    [10];

    /// \brief The default constructor leaves everything uninitialized.
    LocalNumbP3CL() {}
    /// \brief Read indices, boundary-segment numbers and boundary conditions
    /// from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      LocalNumbP3CL(const TetraCL&, const IdxDescCL&, const BndDataT&);

    /// \brief Read indices only
    /// from a tetrahedron.
    LocalNumbP3CL(const TetraCL&, const IdxDescCL&);

    /// \brief Read indices, boundary-segment numbers and boundary conditions
    ///     from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      void
      assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd);

    /// \brief Compute the indices only.
    /// Only num is set up.
    void assign_indices_only (const TetraCL& s, const IdxDescCL& idx);

    /// \brief True, iff index i has a dof associated with it.
    bool WithUnknowns(IdxT i) const { return num[i] != NoIdx; }
};

#if 1
template <typename LocalP3T>
class OswaldProjectionP3AccuCL : public TetraAccumulatorCL
{
  private:
    LocalP3T loc_;
    std::valarray<double>* n_,
                         * n_invalid_;
    bool check_averaging_;
    VecDescCL& avg_;

    LocalNumbP3CL numg;//LocalNumbP2CL

    const VecDescCL*   ls;      // a P3-level-set function
    const BndDataCL<>* lsetbnd; // boundary data for the level set function
    std::valarray<double> ls_loc;
    const PrincipalLatticeCL* lat;

    OswaldProjectionP3AccuCL& set_n (std::valarray<double>* n) { // The clones must refer to the n_ of thread 0.
        n_= n;
        return *this;
    }
    OswaldProjectionP3AccuCL& set_n_invalid (std::valarray<double>* n) { // The clones must refer to the n_invalid of thread 0.
        n_invalid_= n;
        return *this;
    }

  public:
    OswaldProjectionP3AccuCL (LocalP3T loc, VecDescCL& avg)
        : loc_( loc), n_( 0), n_invalid_( 0), check_averaging_( false), avg_( avg), ls( 0), lsetbnd( 0), lat( 0) {}

    OswaldProjectionP3AccuCL& set_check_averaging (bool b= true) {
        check_averaging_= b;
        return *this;
    }

    OswaldProjectionP3AccuCL& set_level_set_function (const VecDescCL* lsarg, const BndDataCL<>* lsetbndarg, const PrincipalLatticeCL* latarg) {
        ls= lsarg;
        lsetbnd= lsetbndarg;
        lat= latarg;
        ls_loc.resize ( lat ? lat->vertex_size () : 0);
        return *this;
    }

    virtual void begin_accumulation   () {
        n_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
        if (check_averaging_)
            n_invalid_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
    }
    virtual void finalize_accumulation() {
        loc_.finalize_accumulation ();
        if (check_averaging_)
            for (size_t i= 0; i < n_->size (); ++i)
                if (n_[0][i] == 0 && n_invalid_[0][i] > 0)
                    std::cerr << "OswaldProjectionP3AccuCL::finalize_accumulation: No local value for " << i << "; invalid_p: " << n_invalid_[0][i] << ".\n";
        delete n_;
        delete n_invalid_;
    }

    virtual void visit (const TetraCL& t) {
        if (ls != 0) {
            LocalP3CL<> locp3_ls( t, *ls, *lsetbnd);
            evaluate_on_vertexes( locp3_ls, *lat, Addr( ls_loc));
            if (equal_signs( ls_loc))
                return;
        }
        loc_.set_tetra( &t);
        numg.assign_indices_only( t, *avg_.RowIdx);
        for (Uint i= 0; i < 10; ++i) {
            if (!numg.WithUnknowns( i))
                continue;
            const IdxT dof= numg.num[i];
        if (loc_.invalid_p (i)) {
            if (check_averaging_)
                ++n_invalid_[0][dof/loc_.num_components];
            continue;
        }
            double& n= n_[0][dof/loc_.num_components]; // This assumes that the local gradient is in the components dof/3..dof/3 + 2.
            n+= 1.;
            typedef typename LocalP3T::value_type value_type;
            const value_type& oldavg= DoFHelperCL<value_type, VectorCL>::get( avg_.Data, dof);
            DoFHelperCL<value_type, VectorCL>::set( avg_.Data, dof, ((n - 1.)/n)*oldavg + (1./n)*loc_[i]);
        }
    }

    virtual TetraAccumulatorCL* clone (int /*clone_id*/) {
        OswaldProjectionP3AccuCL* p= new OswaldProjectionP3AccuCL( loc_, avg_);
        p->set_n( n_)
          .set_n_invalid (n_invalid_)
          .set_check_averaging (check_averaging_)
          .set_level_set_function (ls, lsetbnd, lat);
       return p;
    }
};
#endif

void averaging_P3_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad);


class InterfaceCommonDataP3CL : public TetraAccumulatorCL
{
private:
    InterfaceCommonDataP3CL** the_clones;

    const VecDescCL*   ls;      // P3-level-set
    const BndDataCL<>* lsetbnd; // boundary data for the level set function

    const PrincipalLatticeCL* lat;

    bool compute_quaddomains_;

public:
    /// common data @{
    LocalP3CL<> locp3_ls;

    LocalP3CL<>          p3[10];
    LocalP1CL<Point3DCL> gradrefp3[10];

    std::valarray<double>     ls_loc;
    SurfacePatchCL            surf;
    QuadDomain2DCL            qdom;
    ProjectedQuadDomain2DCL   qdom_projected;
    QuaQuaMapperCL            quaqua;

    const PrincipalLatticeCL& get_lattice () const
    {
        return *lat;
    }
    /// @}

    LocalP3CL<> get_local_p3_ls (const TetraCL& t) const
    {
        return LocalP3CL<>( t, *ls, *lsetbnd);
    };

    const InterfaceCommonDataP3CL& get_clone () const
    {
        const int tid= omp_get_thread_num();
        return tid == 0 ? *this : the_clones[tid][0];
    }

    bool empty () const
    {
        return surf.empty();
    }

    void compute_absdet (bool b)
    {
        qdom_projected.compute_absdets( b);
    }
    void compute_quaddomains (bool b)
    {
        compute_quaddomains_= b;
    }

    void set_lattice (const PrincipalLatticeCL& newlat)
    {
        lat= &newlat;
        ls_loc.resize( lat->vertex_size());
    }

    InterfaceCommonDataP3CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
                             const QuaQuaMapperCL& quaquaarg, const PrincipalLatticeCL& lat_arg);//InterfaceCommonDataP2CL
    virtual ~InterfaceCommonDataP3CL () {}

    virtual void begin_accumulation ()
    {
        the_clones= new InterfaceCommonDataP3CL*[omp_get_max_threads()];
        the_clones[0]= this;
    }
    virtual void finalize_accumulation()
    {
        delete[] the_clones;
    }

    virtual void visit (const TetraCL& t)
    {
        surf.clear();
        locp3_ls.assign( t, *ls, *lsetbnd);
        evaluate_on_vertexes( locp3_ls, *lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            return;
        surf.make_patch<MergeCutPolicyCL>( *lat, ls_loc);
        if (surf.empty())
            return;

        if (compute_quaddomains_)
        {
            make_CompositeQuad5Domain2D ( qdom, surf, t);
            quaqua.set_point( &t, BaryCoordCL()); // set the current tetra.
            qdom_projected.assign( surf, qdom, quaqua);
        }
    }

    virtual InterfaceCommonDataP3CL* clone (int clone_id)
    {
        return the_clones[clone_id]= new InterfaceCommonDataP3CL( *this);
    }
};

class LocalMassP3CL
{
private:
    std::valarray<double> qp3[10];

public:
    static const FiniteElementT row_fe_type= P3IF_FE,
                                col_fe_type= P3IF_FE;

    double coup[10][10];

    void setup (const TetraCL&, const InterfaceCommonDataP3CL& cdata)
    {
#if 0
        SVectorCL<4> tmp(1.0,2.0,3.0,4.0);
        std::cout<< cdata.p3[7](tmp)<<std::endl;
#endif
        for (int i= 0; i < 10; ++i)
            resize_and_evaluate_on_vertexes ( cdata.p3[i], cdata.qdom, qp3[i]);
        //caculate discrete shape function values on triangles

        for (int i= 0; i < 10; ++i)
        {
            coup[i][i]= quad_2D( cdata.qdom_projected.absdets()*qp3[i]*qp3[i], cdata.qdom);
            //cal interations between differerent shape functions
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D( cdata.qdom_projected.absdets()*qp3[j]*qp3[i], cdata.qdom);
        }
    }

    LocalMassP3CL () {}
};


/// \brief Trafo of the interfacial gradient on the linear interface to the quadratic iface under a QuaQuaMapperCL.
/// Computes W from La. 5.1 of the high order paper. The transformation of the gradient in the discretization requires W^{-1}.
//void gradient_trafo (const TetraCL& tet, const BaryCoordCL& xb, const QuaQuaMapperCL& quaqua, const SurfacePatchCL& p, SMatrixCL<3,3>& W);

class LocalLaplaceBeltramiP3CL
{
private:
    double D_; // diffusion coefficient

    LocalP1CL<Point3DCL> gradp3[10];
    GridFunctionCL<Point3DCL> qgradp3[10];

    GridFunctionCL<Point3DCL> nl;
    GridFunctionCL<SMatrixCL<3,3> > Winv;

public:
    static const FiniteElementT row_fe_type= P3IF_FE,
                                col_fe_type= P3IF_FE;

    double coup[10][10];

    const GridFunctionCL<Point3DCL>& get_qgradp3 (size_t i)
    {
        return qgradp3[i];
    }

    void setup (const TetraCL& t, const InterfaceCommonDataP3CL& cdata)
    {
        if (cdata.surf.normal_empty())
            cdata.surf.compute_normals( t);
        resize_and_scatter_piecewise_normal( cdata.surf, cdata.qdom, nl);

        Winv.resize( cdata.qdom.vertex_size());
        QRDecompCL<3,3> qr;
        SVectorCL<3> tmp;
        for (Uint i= 0; i < cdata.qdom.vertex_size(); ++i)
        {
            gradient_trafo( t, cdata.qdom.vertex_begin()[i], cdata.quaqua, cdata.surf, qr.GetMatrix());
            qr.prepare_solve();
            for (Uint j= 0; j < 3; ++j)
            {
                tmp= std_basis<3>( j + 1);
                qr.Solve( tmp);
                Winv[i].col( j, tmp);
            }
        }

        double dummy;
        SMatrixCL<3,3> T;
        GetTrafoTr( T, dummy, t);
        P3DiscCL::GetGradients( gradp3, cdata.gradrefp3, T);
        for (int i= 0; i < 10; ++i)
        {
            resize_and_evaluate_on_vertexes ( gradp3[i], cdata.qdom, qgradp3[i]);
            for (Uint j= 0; j < qgradp3[i].size(); ++j)
            {
                tmp=  qgradp3[i][j] - inner_prod( nl[j], qgradp3[i][j])*nl[j];
                qgradp3[i][j]= Winv[j]*tmp;
            }
        }

        for (int i= 0; i < 10; ++i)
        {
            coup[i][i]= quad_2D( cdata.qdom_projected.absdets()*dot( qgradp3[i], qgradp3[i]), cdata.qdom);
            for(int j= 0; j < i; ++j)
                coup[i][j]= coup[j][i]= quad_2D( cdata.qdom_projected.absdets()*dot( qgradp3[j], qgradp3[i]), cdata.qdom);
        }
    }

    LocalLaplaceBeltramiP3CL (double D)
        :D_( D) {}
};



/// \brief Compute the P3 load vector corresponding to the function f on a single tetra.
class LocalVectorP3CL
{
private:
    instat_scalar_fun_ptr f_;
    double time_;

    std::valarray<double> qp3,
        qf;

public:
    static const FiniteElementT row_fe_type= P3IF_FE;

    double vec[10];

    LocalVectorP3CL (instat_scalar_fun_ptr f, double time) : f_( f), time_( time) {}

    void setup (const TetraCL& t, const InterfaceCommonDataP3CL& cdata, const IdxT numr[10])
    {
        //   coutTet(t);
        resize_and_evaluate_on_vertexes( f_, cdata.qdom_projected.vertexes(), time_, qf);
        qp3.resize( cdata.qdom.vertex_size());
        int IdxPos = -1;
        for (Uint i= 0; i < 10; i++)
        {
            if (numr[i] == NoIdx)
                continue;
            IdxPos = i;
            evaluate_on_vertexes( cdata.p3[i], cdata.qdom, Addr( qp3));
            //  std::valarray<double> qp30{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
            vec[i]= quad_2D( cdata.qdom_projected.absdets()*qf*qp3, cdata.qdom);
            //    cout2txt(vec[i]);
            //   ouput_valarray(qf);
            //     ouput_valarray(qp30);
            //
            //  ouput_valarray( cdata.qdom_projected.absdets());

        }
        //if(IdxPos!=-1)
         //   cout2txt(vec[IdxPos]);
        //   std::cout<<cdata.qdom_projected.absdets()<<std::endl;
        // getchar();
    }
};


//**************************************************************************
// Class:   P3EvalCL                                                       *
// Template Parameter:                                                     *
//          Data     - The result-type of this finite-element-function on  *
//                     the multigrid                                       *
//          _BndData - Class-type that contains functions that describe,   *
//                     how to handle values in boundary-simplices:         *
//                     bool IsOnDirBnd(VertexCL) - iff true, we use        *
//                     Data GetDirBndValue(T) to obtain the function value *
//          _VD      - (const-) VecDescBaseCL<> like type - type of the    *
//                     container that stores numerical data of this finite-*
//                     element function. The class must contain the typedef*
//                     DataType representing the type used for storing     *
//                     numerical data.                                     *
// Purpose: Abstraction that represents boundary-data and VecDescCL-objects*
//          as a function on the multigrid, that can be evaluated on       *
//          vertices, edges, faces and tetrahedrons via val()-functions and*
//          coordinates in the reference tetrahedron.                      *
//          Degrees of freedom can be set, via SetDoF.                     *
//          We provide GetDoF(S) for simplices S to store all relevant     *
//          numerical data via push_back() in a container. This container  *
//          can be passed to a special val() function and allows for faster*
//          evaluation of the FE-function, if several evaluations on the   *
//          same simplex are necessary.                                    *
//          Generally, evaluations on lower-dimensional simplices are      *
//          faster as only a smaller amount of shape-functions has to be   *
//          evaluated.                                                     *
//**************************************************************************
template<class Data, class _BndData, class _VD>
class P3EvalCL
{
public:
    typedef Data     DataT;
    typedef _BndData BndDataCL;
    typedef _VD      VecDescT;

    typedef P3EvalCL<Data, _BndData, _VD> _self;
    typedef P3EvalCL<Data, _BndData, typename ConstHelperCL<_VD>::stripped_type> modifiable_type;
    typedef P3EvalCL<Data, _BndData, typename ConstHelperCL<_VD>::const_type>    const_type;

    typedef LocalP3CL<DataT> LocalFET;

protected:
    // numerical data
    VecDescT*          _sol;
    // boundary-data
    BndDataCL*         _bnd;
    // the multigrid
    const MultiGridCL* _MG;

    inline DataT // helper-function to evaluate on a vertex; use val() instead
    GetDoF(const VertexCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s, _sol->t)
            : DoFHelperCL<DataT,typename VecDescT::DataType>::get(
            _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()));
    }
    inline DataT // helper-function to evaluate on an edge; use val() instead
    GetDoF(const EdgeCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s, _sol->t)
            : DoFHelperCL<DataT,typename VecDescT::DataType>::get(
            _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()));
    }

public:
    P3EvalCL() :_sol( 0), _bnd( 0), _MG( 0){}
    P3EvalCL(VecDescT* sol, BndDataCL* bnd, const MultiGridCL* MG)
        :_sol( sol), _bnd( bnd), _MG( MG) {}
    //default copy-ctor, dtor, assignment-op
    // copying P3EvalCL-objects is safe - it is a flat copy, which is fine,
    // as P3EvalCL does not take possession of the pointed to _sol, _bnd and _MG.

    void // set / get the container of numerical data
    SetSolution(VecDescT* sol) { _sol= sol; }
    VecDescT*
    GetSolution() const { return _sol; }
    void // set / get the container of boundary-data
    SetBndData(BndDataCL* bnd) { _bnd= bnd; }
    BndDataCL*
    GetBndData() const { return _bnd; }
    const MultiGridCL& // the multigrid we refer to
    GetMG() const { return *_MG; }
    Uint // Triangulation level of this function
    GetLevel() const { return _sol->GetLevel(); }
    // The time at which boundary data is evaluated.
    double GetTime() const { return _sol->t; }

    inline bool UnknownsMissing(const TetraCL& t) const;
    // True, iff the function can be evaluated on the given simplex.
    inline bool IsDefinedOn(const VertexCL&) const;
    inline bool IsDefinedOn(const EdgeCL&) const;
    inline bool IsDefinedOn(const TetraCL&, Uint) const;
    inline bool IsDefinedOn(const TetraCL&) const;

    // evaluation on vertices
    inline void // set the degree of freedom in the vertex; fails if, we are on a Dirichlet-boundary
    SetDoF(const VertexCL&, const DataT&);
    template<class _Cont>
      inline void
      GetDoF(const VertexCL& s, _Cont& c) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&) const;
    inline DataT
    val(const VertexCL& s) const;

    // evaluation on edges
    inline void // set the degree of freedom on the edge; fails if, we are on a Dirichlet-boundary
    SetDoF(const EdgeCL&, const DataT&);
    template<class _Cont>
      inline void
      GetDoF(const EdgeCL& s, _Cont& c) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double) const;
    inline DataT
    val(const EdgeCL& s, double v1) const;
    inline DataT // for your convenience: val(edge) == val(edge, 0.5)
    val(const EdgeCL& s) const;

    // evaluation on faces
    template<class _Cont>
      inline void
      GetDoF(const TetraCL& s, Uint i, _Cont& c) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double, double) const;
    inline DataT
    val(const TetraCL& s, Uint i, double v1, double v2) const;

    // evaluation on a tetrahedron
    template<class _Cont>
      inline void
      GetDoF(const TetraCL& s, _Cont& c) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double, double, double) const;
    inline DataT
    val(const TetraCL& s, double v1, double v2, double v3) const;
    inline DataT
    val(const TetraCL& s, const BaryCoordCL&) const;
    template<class _Cont>
      inline static Data
      val(const _Cont&, const BaryCoordCL&);
};



/// \brief Accumulate L2-norms and errors on the higher order zero level.
/// Works for P1IF_FE, P3IF_FE, and C-functions. All functions are evaluated on the P3-levelset.
class InterfaceL2AccuP3CL : public TetraAccumulatorCL
{
private:
    const InterfaceCommonDataP3CL& cdata_;
    const MultiGridCL& mg;
    std::string name_;

    NoBndDataCL<> nobnddata;
    const VecDescCL* fvd;

    instat_scalar_fun_ptr f;
    double f_time;

    LocalLaplaceBeltramiP3CL loc_lb;

    instat_vector_fun_ptr f_grad;
    double f_grad_time;

    InterfaceL2AccuP3CL* tid0p; // The object in OpenMP-thread 0, in which the following variables are updated.
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

    InterfaceL2AccuP3CL (const InterfaceCommonDataP3CL& cdata, const MultiGridCL& mg_arg, std::string name= std::string())
        : cdata_( cdata), mg( mg_arg), name_( name), fvd( 0), f( 0), f_time( 0.),  loc_lb( 1.), f_grad( 0), f_grad_time( 0.) {}
    virtual ~InterfaceL2AccuP3CL () {}

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
        std::cout << "InterfaceL2AccuP3CL::begin_accumulation";
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
        std::cout << "InterfaceL2AccuP3CL::finalize_accumulation";
        if (name_ != std::string())
            std::cout << " for \"" << name_ << "\":";//"P3-solution"
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
        const InterfaceCommonDataP3CL& cdata= cdata_.get_clone();//close cdata
        if (cdata.empty())//empty cdata, finish
            return;

        const int tid= omp_get_thread_num();//multithread

        tid0p->area[tid]+= quad_2D( cdata.qdom_projected.absdets(), cdata.qdom);//for parallelization

        std::valarray<double> qfgrid,//to store solution's dicrete value on triangle
            qf;
        if (fvd != 0)
        {
            // XXX: Check, whether incomplete P3-Data exists locally (which is allowed for P3IF_FE, but not handled correctly by this class --> Extend fvd). Likewise for P1IF_FE...
            if (fvd->RowIdx->GetFE() == P3IF_FE)
                resize_and_evaluate_on_vertexes( make_P3Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);//make_P2Eval
            else if (fvd->RowIdx->GetFE() == P1IF_FE)
                resize_and_evaluate_on_vertexes( make_P1Eval( mg, nobnddata, *fvd), t, cdata.qdom, qfgrid);
//             resize_and_evaluate_on_vertexes( make_P3Eval( mg, nobnddata, *fvd), cdata.qdom_projected, qfgrid);
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
            cout2txt(tmp);
        }

        GridFunctionCL<Point3DCL> qfgradgrid,
                       qfgrad;
        if (fvd != 0)
        {
            qfgradgrid.resize( cdata.qdom.vertex_size());
            if (fvd->RowIdx->GetFE() == P3IF_FE)
            {
                LocalP3CL<> lp3( t, *fvd, nobnddata);
                loc_lb.setup( t, cdata);
                for (Uint i= 0; i < 10; ++i)
                    qfgradgrid+= lp3[i]*loc_lb.get_qgradp3( i);
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

    virtual InterfaceL2AccuP3CL* clone (int /*clone_id*/)
    {
        return new InterfaceL2AccuP3CL( *this);
    }
};


//template <class T, class ResultContT>
//  inline ResultContT&
//  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraBaryPairVectorT& pos, double t, ResultContT& result_container)
//{
//    result_container.resize( pos.size());
//    evaluate_on_vertexes( f, pos, t, sequence_begin( result_container));
//    return result_container;
//}
//
//template <class PEvalT, class ResultIterT>
//  inline ResultIterT
//  evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultIterT result_iterator)
//{
//    typename PEvalT::LocalFET loc_f;
//    const TetraCL* prev_tetra= 0;
//    for (Uint i= 0; i < pos.size(); ++i) {
//        if (prev_tetra != pos[i].first) {
//            prev_tetra= pos[i].first;
//            loc_f.assign( *pos[i].first, f);
//        }
//        *result_iterator++= loc_f( pos[i].second);
//    }
//    return result_iterator;
//}
//
//template <class PEvalT, class ResultContT>
//  inline ResultContT&
//  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraBaryPairVectorT& pos, ResultContT& result_container)
//{
//    result_container.resize( pos.size());
//    evaluate_on_vertexes( f, pos, sequence_begin( result_container));
//    return result_container;
//}
//
//
//
template<class BndData_, class VD_>
  P3EvalCL<typename BndData_::bnd_type, BndData_, VD_>
    make_P3Eval (const MultiGridCL& mg, BndData_& bnd, VD_& vd)
{
    return P3EvalCL<typename BndData_::bnd_type, BndData_, VD_>( &vd, &bnd, &mg);
}

//**************************************************************************
// RestrictP3: Stores the DoF-values of a P3-function corresponding to vd  *
//     and bnd for tetrahedron s in the container c.                       *
// Precondition: vd is a VecDescCL for a P3-function on level l, bnd is a  *
//     BndDataCL and s a tetrahedron on a level <= l. c is a container     *
//     (component access with []) that can hold at least 10 values of f's  *
//     return type.                                                        *
// Postcondition: c contains the value of f in the 10 DoF in the order used*
//     by FE_P3CL.                                                         *
//**************************************************************************
//template <class VecDescT, class BndDataT, class Cont>
//void RestrictP3(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c);//RestrictP2
template <class VecDescT, class BndDataT, class Cont>
void RestrictP3(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c)
{
    const Uint slvl= s.GetLevel();
    const Uint flvl= vd.GetLevel();
    Assert( slvl<=flvl, DROPSErrCL("RestrictP3: Tetra is on a finer level"
            "than the function."), ~0);

    typedef typename VecDescT::DataType VecT;
    typedef DoFHelperCL< typename BndDataT::bnd_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();
    for (Uint i= 0; i<NumVertsC; ++i)
        c[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
    for (Uint i= 0; i<NumEdgesC; ++i) {
        const EdgeCL& e= *s.GetEdge( i);
        c[i+NumVertsC]= (slvl < flvl && e.IsRefined())
            ? ( !bnd.IsOnDirBnd( *e.GetMidVertex())
                ? DoFT::get( v, e.GetMidVertex()->Unknowns( idx))
                : bnd.GetDirBndValue( *e.GetMidVertex(), vd.t))
            : ( !bnd.IsOnDirBnd( e)
                ? DoFT::get( v, e.Unknowns( idx))
                : bnd.GetDirBndValue( e, vd.t));
    }
}

template <class P3FuncT, class Cont>
void RestrictP3(const TetraCL& s, const P3FuncT& f, Cont& c)
{
    RestrictP3( s, *f.GetSolution(), *f.GetBndData(), c);
}


//**************************************************************************
// Class:   P3EvalCL                                                       *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  inline void
  P3EvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
            DROPSErrCL( "P3EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."),
            DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns( _sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P3EvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P3EvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return c[0];
}

template<class Data, class _BndData, class _VD>
  inline Data
  P3EvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF( s);
}


template<class Data, class _BndData, class _VD>
  inline void
  P3EvalCL<Data, _BndData, _VD>::SetDoF(const EdgeCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
            DROPSErrCL( "P3EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-edge."),
            DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns( _sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P3EvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF( *s.GetVertex( 0));
    c[1]= GetDoF( *s.GetVertex( 1));
    c[2]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P3EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0] * FE_P3CL::H0( v1) + c[1] * FE_P3CL::H1( v1) + c[2] * FE_P3CL::H2( v1);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P3EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF( *s.GetVertex(0))*FE_P3CL::H0( v1)
          + GetDoF( *s.GetVertex(1))*FE_P3CL::H1( v1)
          + GetDoF( s)*FE_P3CL::H2( v1);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P3EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s) const
{
    return GetDoF( s);
}


template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P3EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, Uint face, _Cont& c) const
{
    for(Uint i= 0; i < 3; ++i) {
        c[i]= GetDoF( *s.GetVertex( VertOfFace( face, i)));
        c[i+3]= GetDoF( *s.GetEdge( EdgeOfFace( face, i)));
    }
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P3EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2) const
{
    return c[0] * FE_P3CL::H0(v1, v2) + c[1] * FE_P3CL::H1(v1, v2)
         + c[2] * FE_P3CL::H2(v1, v2) + c[3] * FE_P3CL::H3(v1, v2)
         + c[4] * FE_P3CL::H4(v1, v2) + c[5] * FE_P3CL::H5(v1, v2);
}


template<class Data, class _BndData, class _VD>
  inline Data
  P3EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, Uint face, double v1, double v2) const
{
    return  GetDoF( *s.GetVertex( VertOfFace( face, 0)))*FE_P3CL::H0( v1, v2)
           +GetDoF( *s.GetVertex( VertOfFace( face, 1)))*FE_P3CL::H1( v1, v2)
           +GetDoF( *s.GetVertex( VertOfFace( face, 2)))*FE_P3CL::H2( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 0)))*FE_P3CL::H3( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 1)))*FE_P3CL::H4( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 2)))*FE_P3CL::H5( v1, v2);
}


template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P3EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF( *s.GetVertex( i));
    for (Uint i= 0; i < NumEdgesC; ++i)
        c[i+NumVertsC]= GetDoF( *s.GetEdge( i));
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P3EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return  c[0] * FE_P3CL::H0( v1, v2, v3) + c[1] * FE_P3CL::H1( v1, v2, v3)
          + c[2] * FE_P3CL::H2( v1, v2, v3) + c[3] * FE_P3CL::H3 (v1, v2, v3)
          + c[4] * FE_P3CL::H4( v1, v2, v3) + c[5] * FE_P3CL::H5( v1, v2, v3)
          + c[6] * FE_P3CL::H6( v1, v2, v3) + c[7] * FE_P3CL::H7( v1, v2, v3)
          + c[8] * FE_P3CL::H8( v1, v2, v3) + c[9] * FE_P3CL::H9( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P3EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P3CL::H0( v1, v2, v3)
           +GetDoF( *s.GetVertex( 1))*FE_P3CL::H1( v1, v2, v3)
           +GetDoF( *s.GetVertex( 2))*FE_P3CL::H2( v1, v2, v3)
           +GetDoF( *s.GetVertex( 3))*FE_P3CL::H3( v1, v2, v3)
           +GetDoF( *s.GetEdge( 0))*FE_P3CL::H4( v1, v2, v3)
           +GetDoF( *s.GetEdge( 1))*FE_P3CL::H5( v1, v2, v3)
           +GetDoF( *s.GetEdge( 2))*FE_P3CL::H6( v1, v2, v3)
           +GetDoF( *s.GetEdge( 3))*FE_P3CL::H7( v1, v2, v3)
           +GetDoF( *s.GetEdge( 4))*FE_P3CL::H8( v1, v2, v3)
           +GetDoF( *s.GetEdge( 5))*FE_P3CL::H9( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P3EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, const BaryCoordCL& p) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P3CL::H0( p)
           +GetDoF( *s.GetVertex( 1))*FE_P3CL::H1( p)
           +GetDoF( *s.GetVertex( 2))*FE_P3CL::H2( p)
           +GetDoF( *s.GetVertex( 3))*FE_P3CL::H3( p)
           +GetDoF( *s.GetEdge( 0))*FE_P3CL::H4( p)
           +GetDoF( *s.GetEdge( 1))*FE_P3CL::H5( p)
           +GetDoF( *s.GetEdge( 2))*FE_P3CL::H6( p)
           +GetDoF( *s.GetEdge( 3))*FE_P3CL::H7( p)
           +GetDoF( *s.GetEdge( 4))*FE_P3CL::H8( p)
           +GetDoF( *s.GetEdge( 5))*FE_P3CL::H9( p);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P3EvalCL<Data, _BndData, _VD>::val(const _Cont& c, const BaryCoordCL& p)
{
    return  c[0] * FE_P3CL::H0( p) + c[1] * FE_P3CL::H1( p)
          + c[2] * FE_P3CL::H2( p) + c[3] * FE_P3CL::H3( p)
          + c[4] * FE_P3CL::H4( p) + c[5] * FE_P3CL::H5( p)
          + c[6] * FE_P3CL::H6( p) + c[7] * FE_P3CL::H7( p)
          + c[8] * FE_P3CL::H8( p) + c[9] * FE_P3CL::H9( p);
}

template<class Data, class _BndData, class _VD>
inline bool P3EvalCL<Data, _BndData, _VD>::UnknownsMissing(const TetraCL& t) const
{
    const Uint idx= _sol->RowIdx->GetIdx();
    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd();
         it!=end; ++it)
        if ( !IsDefinedOn( **it)) return true;
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd();
         it!=end; ++it)
        if ( !(_bnd->IsOnDirBnd( **it) || (*it)->Unknowns.Exist(idx) )) return true;
    return false;
}


template<class Data, class _BndData, class _VD>
inline bool P3EvalCL<Data, _BndData, _VD>::IsDefinedOn(const VertexCL& v) const
{
    return _bnd->IsOnDirBnd( v)
           || (v.Unknowns.Exist() && v.Unknowns.Exist( _sol->RowIdx->GetIdx()));
}

template<class Data, class _BndData, class _VD>
inline bool P3EvalCL<Data, _BndData, _VD>::IsDefinedOn(const EdgeCL& e) const
{
    return IsDefinedOn( *e.GetVertex( 0)) && IsDefinedOn( *e.GetVertex( 1))
           && (_bnd->IsOnDirBnd( e)
               || (e.Unknowns.Exist() && e.Unknowns.Exist( _sol->RowIdx->GetIdx())));
}

template<class Data, class _BndData, class _VD>
inline bool P3EvalCL<Data, _BndData, _VD>::IsDefinedOn(
    const TetraCL& t, const Uint face) const
{
    const Uint idx= _sol->RowIdx->GetIdx();
    for (Uint i=0; i<3; ++i) {
        const VertexCL& v= *t.GetVertex( VertOfFace( face, i));
        if (!IsDefinedOn( v)) return false;
        const EdgeCL* const ep= t.GetEdge( EdgeOfFace( face, i));
        if (!(_bnd->IsOnDirBnd( *ep)
              || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx))))
            return false;
    }
    return true;
}


template<class Data, class _BndData, class _VD>
inline bool P3EvalCL<Data, _BndData, _VD>::IsDefinedOn(const TetraCL& t) const
{
    for (Uint i=0; i<NumVertsC; ++i)
        if (!IsDefinedOn( *t.GetVertex( i))) return false;
    const Uint idx= _sol->RowIdx->GetIdx();
    for (Uint i=0; i<NumEdgesC; ++i) {
        const EdgeCL* const ep= t.GetEdge( i);
        if (!(_bnd->IsOnDirBnd( *ep)
              || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx))))
            return false;
    }
    return true;
}

//#ifndef CreateNumbOnInterfaceP3
//#define CreateNumbOnInterfaceP3

#if 0
/// \brief Routine to number P3-unknowns on the vertices and edges surrounding an
/// interface.
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are cut by the zero level of lset.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// A more user friendly interface is provided by IdxDescCL::CreateNumbOnInterface.

void CreateNumbOnInterfaceP3 (const Uint idx, IdxT& counter, Uint stride,
        const MultiGridCL::TriangVertexIteratorCL& vbegin,
        const MultiGridCL::TriangVertexIteratorCL& vend,
        const MultiGridCL::TriangEdgeIteratorCL& ebegin,
        const MultiGridCL::TriangEdgeIteratorCL& eend,
        const MultiGridCL::TriangTetraIteratorCL& begin,
        const MultiGridCL::TriangTetraIteratorCL& end,
        const VecDescCL& ls, const BndDataCL<>& lsetbnd, double omit_bound= -1./*default to using all dof*/)
{
    //const size_t stride= 1;

    LocalP3CL<> p3[10];
    for (int i= 0; i < 10; ++i)
        p3[i][i]= 1.;
    // first set NoIdx in all vertices
    for (MultiGridCL::TriangVertexIteratorCL vit= vbegin; vit != vend; ++vit) {
        vit->Unknowns.Prepare(idx);
        vit->Unknowns.Invalidate(idx);
    }
    for (MultiGridCL::TriangEdgeIteratorCL   eit= ebegin; eit != eend; ++eit) {
        eit->Unknowns.Prepare(idx);
        eit->Unknowns.Invalidate(idx);
    }
    // then create numbering of vertices at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);
    std::valarray<double> ls_loc( lat.vertex_size());
    LocalP3CL<> locp3_ls;
    SPatchCL<3> patch;
    QuadDomainCodim1CL<3> qdom;
    std::valarray<double> qp3; // P3-shape-function as integrand
    for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
        locp3_ls.assign( *it, ls, lsetbnd);
        evaluate_on_vertexes( locp3_ls, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        const double limit= omit_bound < 0. ? 0. : omit_bound*std::pow( it->GetVolume()*6, 4./3.);
        patch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        if (patch.empty())
            continue;

        make_CompositeQuad5Domain2D( qdom, patch, *it);
        qp3.resize( qdom.vertex_size());
        for (Uint i= 0; i < 10; ++i) {
            UnknownHandleCL& unknowns= i < 4 ? const_cast<VertexCL*>( it->GetVertex( i))->Unknowns
                                             : const_cast<EdgeCL*>( it->GetEdge( i - 4))->Unknowns;
            if (unknowns.Exist( idx))
                continue;

            evaluate_on_vertexes( p3[i], qdom, Addr( qp3));
            if (quad_codim1( qp3*qp3, qdom) > limit) {
                unknowns( idx)= counter;
                counter+= stride;
            }
        }
    }
}
#endif
#endif
