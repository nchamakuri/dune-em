/////////////////////////////////////////////////////////////////
//       Program to simulate single heartbeat (no S2)          //
//       Used for computing the L2 error for different grid    //
/////////////////////////////////////////////////////////////////

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<sys/stat.h>
#include <unistd.h>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/common/timer.hh>
#include <dune/common/float_cmp.hh> 	

#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/grid/geometrygrid/grid.hh>
#include <dune/grid/common/gridinfo.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/instationaryfilenamehelper.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/finiteelement/l2orthonormal.hh>
#include<dune/pdelab/finiteelement/qkdglagrange.hh>
#include<dune/pdelab/solver/newton.hh>
#include<dune/pdelab/function/discretegridviewfunction.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
// dune-heart includes
#include<dune/electromechanics/common/cardiac_util.hh>
#include<dune/electromechanics/common/ms_regularized.hh>
#include<dune/electromechanics/localoperators/electromechanical_intra_coupled_lop.hh>

#include "schemes.hh"
#include "driver_em.hh"

#define STRUCTURED
#define CUBE

template <int dim>
class GridTransformation
: public Dune :: AnalyticalCoordFunction< double, dim, dim, GridTransformation <dim> >{
    typedef GridTransformation This;
    typedef Dune :: AnalyticalCoordFunction< double, dim, dim, This > Base;

    public:
    typedef typename Base :: DomainVector DomainVector;
    typedef typename Base :: RangeVector RangeVector;

    GridTransformation(int my_rank_):
        giveOutput(1), my_rank(my_rank_){}

    void evaluate(const DomainVector &x, RangeVector &y) const{
        y = x;
    }
  private:
    mutable int giveOutput;
    int my_rank;
};


//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
    {
    if(helper.rank()==0)
      std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }

    // open ini file and parse it in
    ConfigParser::parse((argc >= 2 ? argv[1] : "electromechanical2D.ini"));

    const int dim = 2;
    // initialize the cells and size of the grid
    std::vector<int> nel(2);   // number of grid pints in x and y directions
        std::vector<double> L(2);  // domain size considering [0,L0]x[0,L1]
        std::array<int, 2> Par;
        nel[0] =   ConfigParser::get("global.xcells",80);
        nel[1] =   ConfigParser::get("global.ycells",80);
        L[0] =  ConfigParser::get<double>("global.xmax",1.0);
        L[1] =  ConfigParser::get<double>("global.ymax",1.0);
        int Partitioningx = ConfigParser::get("global.xprocs",1);
        int Partitioningy = ConfigParser::get("global.yprocs",1);
        int refinements = ConfigParser::get("global.refinements",0);
        Par[0] = Partitioningx; Par[1] = Partitioningy;
    
    if(helper.rank() == 0)
       std::cout<<" helper size = "<<helper.size()<<" x and y partition = "<<Par[0]<<" x "<< Par[1]<<std::endl;

	// ========= Setup YaspGrid ======== //
    std::array<std::vector<double>,dim>  coords;
        for (int i = 0; i < dim; i++){
              double h = L[i] / nel[i];
              coords[i].resize(nel[i] + 1);
              coords[i][0] = 0.0;
              for (int j = 1; j <= nel[i]; j++){
                    coords[i][j] = coords[i][j-1] + h;
              }
        }
        std::bitset<dim> periodic(false);
        int overlap = 2;
        typedef Dune::YaspGrid<dim,Dune::TensorProductCoordinates<double,dim>> YGRID;
        YaspPartition<dim,std::array<int,dim>> yp(Par);
        YGRID yaspgrid(coords,periodic,overlap,helper.getCommunicator(),(Dune::YLoadBalance<2>*)&yp);
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
        if(refinements > 0) yaspgrid.globalRefine(refinements);
        typedef YGRID::LeafGridView YGV;
        const YGV ygv = yaspgrid.leafGridView();

        if(helper.rank() == 0){
              Dune::gridinfo(yaspgrid);
              std::cout << "Number of elements per processor: " << ygv.size(0) << std::endl;
              std::cout << "Number of nodes per processor: "    << ygv.size(2) << std::endl;
        }

    // ========= Grid Transformation ======== //
    typedef GridTransformation<dim> GRID_TRAFO;
        GRID_TRAFO gTrafo(helper.rank());

        typedef typename Dune::GeometryGrid<YGRID,GRID_TRAFO> GRID;
        GRID grid(yaspgrid,gTrafo);
        if(helper.rank() == 0)
            std::cout << "Grid transformation complete" << std::endl;

        //Define Grid view
    typedef typename GRID::LeafGridView GV;
    const GV gv = grid.leafGridView();
    if(helper.rank() == 0)
         std::cout << "Grid view set up complete" << std::endl;
        

#ifdef CUBE 
    driver_em(gv,TaylorHood_21_Quadrilateral(gv),helper);
#else
   driver_em(gv,TaylorHood_21_Triangle(gv),helper);
#endif
   }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
  }
}
