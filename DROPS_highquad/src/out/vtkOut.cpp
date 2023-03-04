/// \file vtkOut.cpp
/// \brief solution output in VTK XML format
/// \author LNM RWTH Aachen: Jens Berger, Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "out/vtkOut.h"
#include "geom/simplex.h"
#include "geom/deformation.h"
#include "misc/base64.h"

namespace DROPS
{

VTKOutCL::VTKOutCL(const MultiGridCL& mg, const std::string& dataname, Uint numsteps,
                   const std::string& dirname, const std::string& filename,
                   const std::string& pvdfilename, bool binary, bool onlyP1, bool P2DG,
                   Uint lvl, bool reusepvd, bool usedeformed)
/** Beside constructing the VTKOutCL, this function computes the number of
    digits, that are used to decode the time steps in the filename.
\param mg        MultiGridCL that contains the geometry
\param dataname  name of the data
\param numsteps  number of time steps
\param filename  prefix of all files (e.g. vtk/output)
\param binary    Write out files in binary format.
\param onlyP1    output special for P1 FE (less output)
\param P2DG      output for DG FE (discontinuous data)
\param lvl       Multigrid level
*/
    : mg_(mg), timestep_(0), numsteps_(numsteps), descstr_(dataname),
      dirname_(dirname), filename_(filename), pvdfilename_(pvdfilename),
      binary_(binary), onlyP1_(onlyP1), P2DG_(P2DG), geomwritten_(false),
      vAddrMap_(), eAddrMap_(), coords_(), tetras_(), lvl_(lvl),
      numPoints_(0), numTetras_(0), reusepvd_(reusepvd), usedeformed_(usedeformed)
{
    if (!dirname.empty() && *dirname.rbegin()!='/' )
        dirname_+= '/';
    decDigits_= 1;
    while( numsteps>9){ ++decDigits_; numsteps/=10; }
    // reference coordinates for local P2 function
    // todo: clean up
    refcoords[1][0] = 1.0;
    refcoords[2][1] = 1.0;
    refcoords[3][2] = 1.0;
    refcoords[4][0] = 0.5;
    refcoords[5][0] = refcoords[5][1] = 0.5;
    refcoords[6][1] = 0.5;
    refcoords[7][2] = 0.5;
    refcoords[8][0] = refcoords[8][2] = 0.5;
    refcoords[9][1] = refcoords[9][2] = 0.5;
    for (int i = 0; i < 10; ++i)
    {
        double tmp=1.0;
        for (int j = 0; j < 3; ++j)
        {
            refcoords_b[i][j] = refcoords[i][j];
            tmp -= refcoords[i][j];
        }
        refcoords_b[i][3] = tmp;
    }

}

VTKOutCL::~VTKOutCL ()
{
    for (std::map<std::string,VTKVariableCL*>::iterator it= vars_.begin(); it != vars_.end(); ++it)
        delete it->second;
}

void VTKOutCL::Register (VTKVariableCL& var)
{
    if( vars_.find(var.varName()) != vars_.end())
        std::cout << "Error! Variable name is used twice! No registration of the variable is carried out!" << std::endl;
    else
        vars_[var.varName()]= &var;
}

void VTKOutCL::Write ( double time, bool writeDistribution)
{
    PutGeom( time, writeDistribution);
    for( std::map<std::string, VTKVariableCL*>::iterator it= vars_.begin(); it != vars_.end(); ++it) {
        it->second->put( *this);
    }
    Commit();
}

void VTKOutCL::AppendTimecode( std::string& str) const
/** Appends a time-code to the filename*/
{
    char format[]= "%0Xi",
         postfix[8];
    format[2]= '0' + char(decDigits_);
    std::sprintf( postfix, format, timestep_);
    str+= postfix;
}

void VTKOutCL::CheckFile( const std::ofstream& os) const
/** Checks if a file is open*/
{
    if (!os)
        throw DROPSErrCL( "VTKOutCL: error while opening file!");
}

void VTKOutCL::NewFile(double time, __UNUSED__ bool writeDistribution)
/** Each process opens a new file and writes header into it*/
{
    std::string filename(filename_);
#ifdef _PAR
   ProcCL::AppendProcNum(filename);
   filename+="_";
#endif
    AppendTimecode(filename);
    filename+= ".vtu";
    file_.open((dirname_+filename).c_str());
    if ( !file_){
        CreateDirectory( dirname_);
        file_.open((dirname_+filename).c_str());
    }
    CheckFile( file_);
    PutHeader();
    wrotePointDataLine_= false;
// The file that links the data from all the separate (but nevertheless valid) XML VTK files is exclusively generated by the master-processor
#ifdef _PAR
    IF_MASTER {
        std::string masterfilename( filename_);
        AppendTimecode( masterfilename);
        masterfilename += ".pvtu";
        std::ofstream masterfile( (dirname_+masterfilename).c_str());
        masterfile<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                  <<"<PUnstructuredGrid GhostLevel=\"0\">\n"
                  <<"\t<PPoints>\n"
                  <<"\t\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\" "<<(binary_? "format=\"binary\"":"format=\"ascii\"")<<"/>\n"
                  <<"\t</PPoints>";
        WriteVarNames( masterfile, true);
        const VectorBaseCL<float> x;
        for( VTKvarMapT::iterator it=vars_.begin(); it!=vars_.end(); it++)
        {
            WriteValues(x, it->first, it->second->GetDim(), &masterfile);
        }
        masterfile<<"\n\t</PPointData>"
                  <<"\n\t<PCells>\n"
                  <<"\t\t<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n"
                  <<"\t\t<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n"
                  <<"\t\t<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>\n"
                  <<"\t</PCells>\n";
        if (writeDistribution)
            masterfile << "\t<PCellData Scalars=\"processor\">\n"
                       << "\t\t<PDataArray type=\"Int32\" Name=\"processor\" format=\"ascii\"/>\n"
                       << "\t</PCellData>\n";

        for( int p = 0; p < ProcCL::Size(); p++)
        {
            std::string parfilename;
            std::stringstream helper;
            helper << filename_ << "." << std::setfill('0') << std::setw( int( log10( (float)ProcCL::Size()))+1) << p << "_";

            parfilename=helper.str();
            AppendTimecode( parfilename);
            parfilename+= ".vtu";
            masterfile << "<Piece Source=\"" << parfilename << "\"/>\n";
        }
        masterfile << "</PUnstructuredGrid>\n"
                   << "</VTKFile>";
        GenerateTimeFile( time, masterfilename);
    }
#else
        GenerateTimeFile( time, filename);
#endif
}

void VTKOutCL::GenerateTimeFile( double time, const std::string & name) const
{
    std::string timefilename(pvdfilename_);
    timefilename+=".pvd";
    timefilename=dirname_+timefilename;
    if(timestep_==0 && !reusepvd_)
    {
        std::ofstream timefile(timefilename.c_str());
        timefile << "<?xml version=\"1.0\"?>\n"
                    "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                    "<Collection>\n"
                    "\t<DataSet timestep=\""<< time <<"\" group=\"\" part=\"0\" file=\"" << name <<"\"/>\n"
                    "</Collection>\n"
                    "</VTKFile>";
        timefile.close();
    }
    else
    {
        std::ofstream timefile;
        timefile.open(timefilename.c_str(),std::ios_base::in);
        if(timefile.is_open())
        {
            timefile.seekp(-24,std::ios_base::end);
            timefile << "\t<DataSet timestep=\""<< time <<"\" group=\"\" part=\"0\" file=\"" << name <<"\"/>\n"
                        "</Collection>\n"
                        "</VTKFile>";
            timefile.close();
        }
        else std::cerr << "Could not open VTK Timefile!" << std::endl;
    }
}

void VTKOutCL::PutHeader()
/** Writes the header into the VTK file*/
{
    file_ << "<?xml version=\"1.0\"?>\n"     // this is just the XML declaration, it's unnecessary for the actual VTK file
             "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
             "<UnstructuredGrid>\n";
}

void VTKOutCL::PutFooter()
/** Closes the file XML conform*/
{
    file_ <<"\n\t</PointData>"
            "\n</Piece>"
            "\n</UnstructuredGrid>"
            "\n</VTKFile>";
}

void VTKOutCL::GatherCoord()
/** Iterates over all vertices and edges and assigns a consecutive number and their coordinates to them */
{
    // Get number of vertices and edges in last triangulation level
    Uint numVertices=0, numEdges=0;

    for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(lvl_); it!=mg_.GetTriangVertexEnd(lvl_); ++it)
        numVertices++;
    for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(lvl_); it!=mg_.GetTriangEdgeEnd(lvl_); ++it)
        numEdges++;

    numTetras_=std::distance(mg_.GetTriangTetraBegin(lvl_),mg_.GetTriangTetraEnd(lvl_));  // calculates the number of tetrahedra

    if (!P2DG_)
        numPoints_= numLocPoints_= numVertices + numEdges;
    else
        numPoints_= numLocPoints_ = 10 * numTetras_;

    coords_.resize(3*numLocPoints_);

    if (P2DG_ && usedeformed_)
      throw DROPSErrCL("P2DG and usedeformed in VTK!");

    ///\todo instead of v/eAddrMap_, one could use a P2 index instead (cf. EnsightOutCL)
    if (!P2DG_)
    {
        if (usedeformed_)
        {
            MeshDeformationCL & md = MeshDeformationCL::getInstance();
            ///\todo instead of v/eAddrMap_, one could use a P2 index instead (cf. EnsightOutCL)
            Uint counter=0;
            for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(lvl_); it!=mg_.GetTriangVertexEnd(lvl_); ++it){
                // store a consecutive number for the vertex
                vAddrMap_[&*it]= counter;

                // Put coordinate of the vertex into the field of coordinates
                for (int i=0; i<3; ++i)
                    coords_[3*counter+i]= (float)md.GetTransformedVertexCoord(*it)[i]; // (float)it->GetCoord()[i];
                ++counter;
            }

            for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(lvl_); it!=mg_.GetTriangEdgeEnd(lvl_); ++it){
                // store a consecutive number for the edge
                eAddrMap_[&*it]= counter;

                // Put coordinate of the barycenter of the edge into the field of coordinates
                // const Point3DCL baryCenter= GetBaryCenter(*it);
                for (int i=0; i<3; ++i)
                    coords_[3*counter+i]= (float)md.GetTransformedEdgeBaryCenter(*it)[i]; // (float)baryCenter[i];
                ++counter;
            }
        }
        else
        {
            ///\todo instead of v/eAddrMap_, one could use a P2 index instead (cf. EnsightOutCL)
            Uint counter=0;
            for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(lvl_); it!=mg_.GetTriangVertexEnd(lvl_); ++it){
                // store a consecutive number for the vertex
                vAddrMap_[&*it]= counter;
                // Put coordinate of the vertex into the field of coordinates
                for (int i=0; i<3; ++i)
                    coords_[3*counter+i]= (float)it->GetCoord()[i];
                ++counter;
            }

            for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(lvl_); it!=mg_.GetTriangEdgeEnd(lvl_); ++it){
                // store a consecutive number for the edge
                eAddrMap_[&*it]= counter;

                // Put coordinate of the barycenter of the edge into the field of coordinates
                const Point3DCL baryCenter= GetBaryCenter(*it);
                for (int i=0; i<3; ++i)
                    coords_[3*counter+i]= (float)baryCenter[i];
                ++counter;
            }
        }
    }
    if (P2DG_)
    {
        ///\todo instead of v/eAddrMap_, one could use a P2 index instead (cf. EnsightOutCL)
        Uint counter=0;
        Point3DCL mid(1.0/3.0);
        const double inshift = 1e-6;
        for (MultiGridCL::const_TriangTetraIteratorCL it= mg_.GetTriangTetraBegin(lvl_); it!=mg_.GetTriangTetraEnd(lvl_); ++it){
            // store a consecutive number for the edge
            tAddrMap_[&*it]= counter;
            for (int i = 0; i < 10; ++i)
            {
                Point3DCL tmp = (1-inshift)*refcoords[i] + inshift*mid;
                Point3DCL p = GetWorldCoord(*it,tmp);//refcoords[i]);
                for (int j=0; j<3; ++j)
                    coords_[3*counter+j]= (float)p[j];
                counter++;
            }
            // Put coordinate of the barycenter of the edge into the field of
        }
    }

}

template <class T>
void WriteBase64( const VectorBaseCL<T>& x, std::ostream& os)
/** Writes out a base64 encoded data-stream */
{
    const size_t s= x.size()*sizeof( T);
    const unsigned char* p= reinterpret_cast<const unsigned char*>( &s);
    // If sizeof( size_t) > 4, this only works in little endian... which we hard-code in the vtu-file in any case.
    Base64Encoding::encode( p, p + 4, os, /*wrap_lines*/ false);
    p= reinterpret_cast<const unsigned char*>( Addr( x));
    Base64Encoding::encode( p, p + s, os, /*wrap_lines*/ false);
}

void VTKOutCL::WriteCoords()
/** Each process writes out its coordinates. */
{
    file_<< "<Piece NumberOfPoints=\""<<numPoints_<<"\" NumberOfCells=\""<<numTetras_<<"\">"
            "\n\t<Points>"
            "\n\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << ( binary_ ? "binary\">\n\t\t" : "ascii\">\n\t\t");

    if (binary_)
        WriteBase64(coords_, file_);
    else
        for (Uint i=0; i<numPoints_; ++i)
            file_<< coords_[3*i+0] << ' ' << coords_[3*i+1] << ' ' << coords_[3*i+2]<< ' ';

    file_<< "\n\t\t</DataArray> \n"
            "\t</Points>\n";
}

void VTKOutCL::GatherTetra()
/** Gathers tetrahedra in an array*/
{
    numTetras_=std::distance(mg_.GetTriangTetraBegin(lvl_),mg_.GetTriangTetraEnd(lvl_));  // calculates the number of tetrahedra
    if (P2DG_)
        tetras_.resize((10)*numTetras_);      //
    else
        tetras_.resize((onlyP1_?4:10)*numTetras_);      // (four vertices) * Number of Tetrahedra respectively (four vertices + six edges) * Number of Tetrahedra

    // Gathers connectivities
    Uint counter=0;
    for (MultiGridCL::const_TriangTetraIteratorCL it= mg_.GetTriangTetraBegin(lvl_); it!=mg_.GetTriangTetraEnd(lvl_); ++it){ //loop over all tetrahedra
        if (P2DG_){
            Uint offset = tAddrMap_[&*it];
            for (int i= 0; i<10; ++i)
                tetras_[counter++] = offset + i;
        }
        else
        {
            for (int vert= 0; vert<4; ++vert)
                tetras_[counter++] = vAddrMap_[it->GetVertex(vert)];
            if(!onlyP1_)
            {
                for (int eddy=0; eddy<6; ++eddy)
                    tetras_[counter++] = eAddrMap_[it->GetEdge(eddy)];
            std::swap(tetras_[counter-4],tetras_[counter-5]);    // Permutation needed to make DROPS and VTK compatible (different numeration)
            }
        }
    }
    Assert(counter==(onlyP1_? 4:10)*numTetras_, DROPSErrCL("VTKOutCL::GatherTetra: Mismatching number of tetrahedra"), ~0);
}

void VTKOutCL::WriteTetra( bool writeDistribution)
/** Writes the tetrahedra into the VTK file*/
{
    file_   << "\t<Cells>\n"
               "\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
// Binary output for the connectivity data seems useless (using >5 byte per integer), because it only blows up the amount of needed storage space, but it's implemented anyway,
// for the sake of completeness.
//  if(binary_)
//  {
//      file_<<"binary\">"/*\n\t\t*/;
//      // Write out connectivities
//      WriteBase64(tetras_, file_);
//  }
//  else
//  {
        file_   <<"ascii\">\n\t\t";
        // Write out connectivities
        if(onlyP1_)
            for (Uint i=0; i<numTetras_; ++i)
            {
                file_ << tetras_[4*i+0] << ' '<< tetras_[4*i+1] << ' '<< tetras_[4*i+2] << ' '<< tetras_[4*i+3] << " ";
            }
        else
            for (Uint i=0; i<numTetras_; ++i)
            {
                file_ << tetras_[10*i+0] << ' '<< tetras_[10*i+1] << ' '<< tetras_[10*i+2] << ' '<< tetras_[10*i+3] << ' '
                      << tetras_[10*i+4] << ' '<< tetras_[10*i+5] << ' '<< tetras_[10*i+6] << ' '<< tetras_[10*i+7] << ' '
                      << tetras_[10*i+8] << ' '<< tetras_[10*i+9] << " ";
            }
//  }
    file_ << "\n\t\t</DataArray>\n"
             "\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t";
    if(onlyP1_)
        for(Uint i=1; i<=numTetras_; ++i) file_ << i*4<<" ";
    else
        for(Uint i=1; i<=numTetras_; ++i) file_ << i*10<<" ";
    file_ << "\n\t\t</DataArray>"
             "\n\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t";
    const char* tetraType= (onlyP1_? "10 ":"24 ");
    for(Uint i=1; i<=numTetras_; ++i)
        file_ << tetraType;
    file_ << "\n\t\t</DataArray>"
             "\n\t</Cells>";

    if ( writeDistribution)
        WriteDistribution();
}

void VTKOutCL::WriteDistribution()
/** Writes the distribution-data into the file (as CellData)*/
{
#ifdef _PAR
    file_ << '\n'
          << "\t<CellData>\n"
          << "\t\t<DataArray type=\"Int32\" Name=\"processor\" format=\"ascii\">\n"
          << "\t\t";
    int c=ProcCL::MyRank();
    for( Uint i=0; i < numTetras_; ++i)
            file_<< c << " ";
    file_   << "\n\t\t</DataArray>\n"
            << "\t</CellData>\n";
#endif
}

void VTKOutCL::WriteVarNames(std::ofstream& file, bool masterfile)
{
    std::vector<std::string> scalarvalued;
    std::vector<std::string> vectorvalued;

    for (VTKvarMapT::iterator it= vars_.begin(); it != vars_.end(); ++it)
    {
        if (it->second->GetDim()==1) scalarvalued.push_back(it->first);
        if (it->second->GetDim()==3) vectorvalued.push_back(it->first);
    }
    file << "\n\t<" << (masterfile? "P":"") << "PointData ";

    if(!scalarvalued.empty())
    {
        file << "Scalars=\"" << scalarvalued[0];
        for (size_t i= 1; i < scalarvalued.size(); ++i)
            file << "," << scalarvalued[i];
        file << "\"";
    }
    if (!vectorvalued.empty())
    {
        if (!scalarvalued.empty())
            file << " ";
        file << "Vectors=\"" << vectorvalued[0];
        for (size_t i= 1; i<vectorvalued.size();++i)
            file << "," << vectorvalued[i];
        file << "\"";
    }

    file << ">";
}

void VTKOutCL::WriteValues( const VectorBaseCL<float>& allData, const std::string& name, int numData, std::ofstream* filePtr)
/** Writes out the calculated numerical data*/
{
    std::ofstream& file= filePtr ? *filePtr : file_;

    file << "\n\t\t<" << ( !filePtr ? "" : "P") << "DataArray type=\"Float32\" Name=\"" << name << "\""
            " NumberOfComponents=\"" << numData << "\" format=\""
         << ( binary_ ? "binary\"" : "ascii\"") << ( !filePtr ? "" : "/") << ">";
    if( !filePtr) // omit for master file
    {
        file << "\n\t\t";
        if ( binary_)
            WriteBase64(allData, file_);
        else {
            for ( Uint i=0; i<allData.size(); ++i)
                file << allData[i] << ' ';
        }
        file << "\n\t\t</DataArray>";
    }
}

void VTKOutCL::PutGeom(double time, bool writeDistribution)
/** At first the geometry is put into the VTK file. Therefore this procedure
    opens the file and writes description into the file.
    \param writeDistribution Flag indicator whether distribution-data should be written in the file (as CellData)
*/
{
    NewFile(time,writeDistribution);
    Clear();
    GatherCoord();
    GatherTetra();
    WriteCoords();
    WriteTetra( writeDistribution);
    WriteVarNames(file_,/*masterfile=*/0);
}

void VTKOutCL::Clear()
{
    if (numPoints_>0){

        vAddrMap_.clear();
        eAddrMap_.clear();
        tAddrMap_.clear();
        coords_.resize(0);
        tetras_.resize(0);
    }
}

} // end of namespace DROPS
