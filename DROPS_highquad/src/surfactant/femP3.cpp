#include <surfactant/femP3.h>


using namespace DROPS;
//**************************************************************************
// Class:   FE_P3CL                                                        *
//**************************************************************************
inline double
FE_P3CL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      case 4: return H4( v1, v2, v3);
      case 5: return H5( v1, v2, v3);
      case 6: return H6( v1, v2, v3);
      case 7: return H7( v1, v2, v3);
      case 8: return H8( v1, v2, v3);
      case 9: return H9( v1, v2, v3);
      default: throw DROPSErrCL("FE_P3CL::H: Invalid shape function.");
    };
}

inline double
FE_P3CL::H(Uint dof, const BaryCoordCL& p)
{
    switch(dof) {
      case 0: return H0( p);
      case 1: return H1( p);
      case 2: return H2( p);
      case 3: return H3( p);
      case 4: return H4( p);
      case 5: return H5( p);
      case 6: return H6( p);
      case 7: return H7( p);
      case 8: return H8( p);
      case 9: return H9( p);
      default: throw DROPSErrCL("FE_P3CL::H: Invalid shape function.");
    };
}

inline SVectorCL<3>
FE_P3CL::DH0Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret(4.*(v1 + v2 + v3) -3.);
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH1Ref(double v1, double, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v1 -1.; ret[1]= ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH2Ref(double, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= 0.; ret[1]= 4*v2 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH3Ref(double, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 0.; ret[2]= 4.*v3 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH4Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*( 1. - (2.*v1 + v2 + v3) ); ret[1]= ret[2]= -4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH5Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= -4.*v2; ret[1]= 4.*( 1. -(2.*v2 + v1 + v3) );
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH6Ref(double v1, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v2; ret[1]= 4.*v1; ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH7Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[1]= -4.*v3; ret[2]= 4.*( 1. -(2.*v3 + v1 + v2) );
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH8Ref(double v1, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v3; ret[1]= 0.; ret[2]= 4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH9Ref(double, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 4.*v3; ret[2]= 4.*v2;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DHRef(Uint dof, double v1, double v2, double v3)
{
    switch (dof)
    {
      case 0: return DH0Ref(v1, v2, v3);
      case 1: return DH1Ref(v1, v2, v3);
      case 2: return DH2Ref(v1, v2, v3);
      case 3: return DH3Ref(v1, v2, v3);
      case 4: return DH4Ref(v1, v2, v3);
      case 5: return DH5Ref(v1, v2, v3);
      case 6: return DH6Ref(v1, v2, v3);
      case 7: return DH7Ref(v1, v2, v3);
      case 8: return DH8Ref(v1, v2, v3);
      case 9: return DH9Ref(v1, v2, v3);
      default: throw DROPSErrCL("FE_P3CL::DHRef: Invalid shape function.");
    };
}

inline double
FE_P3CL::Laplace(Uint dof, const SMatrixCL<3,3>& M)
{
    double ret= 0.;
    for (Uint i=0; i<3; ++i)
      for(Uint j=0; j<3; ++j)
        for (Uint k=0; k<3; ++k)
        {
            ret+= M(i,j)*M(i,k)*D2HRef(dof, j, k);
        }
    return ret;
}




//**************************************************************************
// Class:   LocalP3CL                                                      *
//**************************************************************************

template<class T>
void ExtendP1onChild( const LocalP3CL<T>& isoP3, int child, LocalP3CL<T>& P1onParent)
{
    static const int childOrder[8]= { 0, 1, 7, 6, 2, 4, 5, 3};
    const Uint ch= childOrder[child];
    // children ordered such that
    // A) children 0,...,3 located at corners of parent, parent vertex ch is also a child vertex
    // B) children 4,...,7 located inside parent forming an octahedron, one face is on parent face ch-4
    // Note: data can be looked up in topo.h/cpp

    // first step: compute values on parent's vertices
    if (ch<4) // case A
    {
        const Uint pv= ch; // parent vertex which is also vertex of child

        P1onParent[pv]= isoP3[pv]; // copy value for vertex pv
        // 1D interpolation along the 3 edges starting at parent vertex pv
        for (Uint v=0; v<4; ++v)
            if (v != pv)
                P1onParent[v]= 2*isoP3[EdgeByVert(v,pv)+4] - isoP3[pv];
    }
    else // case B
    {
        const Uint pf= ch - 4; // parent face containing one of the child's faces
        // parent face contains 3 vertices of child, remaining fourth vertex is located on parent's edge pe (= edge 1, 4, 1, 4).
        const Uint ev= (pf + 2) % 4; // parent vertex in face pf which is located on parent edge pe (= vertex 2, 3, 0, 1)

        // 2D interpolation on parent face
        for (Uint i=0; i<3; ++i) {
            const Uint fv= VertOfFace(pf,i);
            P1onParent[fv]= 0.;
            for (Uint j=0; j<3; ++j) {
                const Uint fe= EdgeOfFace(pf,j);
                const T val= isoP3[fe+4];
                if (fv == VertOfEdge(fe,0))
                    P1onParent[fv]+= val;
                else if (fv == VertOfEdge(fe,1))
                    P1onParent[fv]+= val;
                else // vert fv is opposite to edge fe in face pf
                    P1onParent[fv]-= val;
             }
        }

        // 1D interpolation along edge pe whose midvertex is remaining 4th vertex of child
        P1onParent[OppVert(pf)]= 2*isoP3[EdgeByVert(OppVert(pf),ev)+4] - P1onParent[ev];
    }

    // second step: linear interpolation of edge values
    for (Uint e=0; e<6; ++e)
        P1onParent[e+4]= 0.5*(P1onParent[VertOfEdge(e,0)] + P1onParent[VertOfEdge(e,1)]);
}

void P3DiscCL::GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[10])
{
    for (int i= 0; i < 10; ++i)
    {
        GRef[i][0]= FE_P3CL::DHRef( i, 0,0,0);//FE_P2CL
        GRef[i][1]= FE_P3CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P3CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P3CL::DHRef( i, 0,0,1);
    }
}

void P3DiscCL::GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
    {
        GRef[i][0]= FE_P3CL::DHRef( i, 0,0,0);
        GRef[i][1]= FE_P3CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P3CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P3CL::DHRef( i, 0,0,1);
        GRef[i][4]= FE_P3CL::DHRef( i, 0.25,0.25,0.25);
    }
}

void P3DiscCL::GetGradientsOnRef( Quad5CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
        for (int j=0; j<Quad5DataCL::NumNodesC; ++j)
        {
            const BaryCoordCL& Node= Quad5DataCL::Node[j];
            GRef[i][j]= FE_P3CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P3DiscCL::GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[10],
    const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<Point3DCL>::SetInterface( p, NodeInTetra);
    for (int i= 0; i < 10; ++i)
        for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
            const BaryCoordCL& Node= NodeInTetra[j];
            GRef[i][j]= FE_P3CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P3DiscCL::GetP3Basis( LocalP3CL<> p3[10])
{
    for (int i= 0; i < 10; ++i) {
        for (int j= 0; j < 10; ++j)
            p3[i][j]= 0;
        p3[i][i]= 1;
    }
}

void P3DiscCL::GetP3Basis( Quad5_2DCL<> p3[10], const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<>::SetInterface( p, NodeInTetra);
    for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
        const BaryCoordCL& Node= NodeInTetra[j];
        p3[0][j]= FE_P3CL::H0( Node);
        p3[1][j]= FE_P3CL::H1( Node);
        p3[2][j]= FE_P3CL::H2( Node);
        p3[3][j]= FE_P3CL::H3( Node);
        p3[4][j]= FE_P3CL::H4( Node);
        p3[5][j]= FE_P3CL::H5( Node);
        p3[6][j]= FE_P3CL::H6( Node);
        p3[7][j]= FE_P3CL::H7( Node);
        p3[8][j]= FE_P3CL::H8( Node);
        p3[9][j]= FE_P3CL::H9( Node);
    }
}

void LocalP3GradientCL::set_tetra (const TetraCL* t)
{
    double det; // dummy
    GetTrafoTr( M, det, *t);
    P3DiscCL::GetGradients( GradLP1, GradRefLP1, M);
    p3.assign( *t, f_, fbnd_);
    p1grad= Point3DCL();
    for (Uint i= 0; i < 10; ++i)
        p1grad+= p3[i]*GradLP1[i];
    for (Uint i= 0; i < 4; ++i)
        p3grad[i]= p1grad[i];
    for (Uint i= 0; i < 6; ++i)
        p3grad[i+4]= 0.5*(p1grad[VertOfEdge( i, 0)] + p1grad[VertOfEdge( i, 1)]);

}


void
LocalNumbP3CL::assign_indices_only (const TetraCL& s, const IdxDescCL& idx)
{
    const Uint sys= idx.GetIdx();
    if (!idx.IsDG())
    {
        for (Uint i= 0; i < 4; ++i)
            num[i]= s.GetVertex( i)->Unknowns.Exist( sys) ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
        for(Uint i= 0; i < 6; ++i)
            num[i+4]= s.GetEdge( i)->Unknowns.Exist( sys) ? s.GetEdge( i)->Unknowns( sys)   : NoIdx;
    }
    else
    {
        Uint first = s.Unknowns(sys);
        for (int i = 0; i < 10; ++i)
            num[i] = first++;
    }
}

LocalNumbP3CL::LocalNumbP3CL(const TetraCL& s, const IdxDescCL& idx)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL-object to be used.
{
    this->assign_indices_only( s, idx);
}


void averaging_P3_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad)
{
    LocalP3GradientCL loc( f, fbnd);
    OswaldProjectionP3AccuCL<LocalP3GradientCL> accu( loc, grad);//OswaldProjectionP2AccuCL
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accumulate( accus, mg, f.RowIdx->TriangLevel(), f.RowIdx->GetBndInfo());
};


//// Create a P3EvalCL without the agonizing template-pain.
//template<class BndData_, class VD_>
//  P3EvalCL<typename BndData_::bnd_type, BndData_, VD_>
//    make_P3Eval (const MultiGridCL& mg, BndData_& bnd, VD_& vd)
//{
//    return P3EvalCL<typename BndData_::bnd_type, BndData_, VD_>( &vd, &bnd, &mg);
//}




InterfaceCommonDataP3CL::InterfaceCommonDataP3CL (const VecDescCL& ls_arg, const BndDataCL<>& lsetbnd_arg,
        const QuaQuaMapperCL& quaquaarg, const PrincipalLatticeCL& lat_arg)
    : ls( &ls_arg), lsetbnd( &lsetbnd_arg), lat( &lat_arg),  compute_quaddomains_( true), ls_loc( lat->vertex_size()),
      quaqua( quaquaarg)
{
    qdom_projected.compute_absdets( true);
    P3DiscCL::GetGradientsOnRef( gradrefp3);
    for (Uint i= 0; i < 10 ; ++i)
        p3[i][i]= 1.; // P3-Basis-Functions
}


