
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-ofset: 2 -*-

#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
//#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/grid/geometrygrid/grid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/common/mcmgmapper.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#endif
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

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelement/qkdglagrange.hh>
#include<dune/pdelab/finiteelement/l2orthonormal.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/common/shared_ptr.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/localoperator/l2.hh>
#include <dune/grid/common/gridinfo.hh> // definition of gridinfo
#include<dune/pdelab/solver/newton.hh>

// dune-em includes
#include<dune/electromechanics/common/cardiac_util.hh>
#include<dune/electromechanics/common/helper.hh>
#include<dune/electromechanics/common/ms_regularized.hh>
#include<dune/electromechanics/common/output.hh>
#include<dune/electromechanics/localoperators/disc.hh>
#include<dune/electromechanics/common/odesolvers.hh>
#include<dune/electromechanics/localoperators/taylorhoodmechanicallop.hh>
#include<dune/electromechanics/elasticity_setup/dirichletBC.hh>

#define STRUCTURED
#define CUBE
#define USE_YASPGRID 1
//#define UGGRID
#define USE_ALUGRID 0

#include "finiteelasticity.hh"
#include "electromechanical.hh"
 

template<class GV>
class BiDHelper
{
public:
    enum{dim = GV::dimension,porder =1, qorder=2, isParallel=!Dune::MPIHelper::isFake};
    typedef typename GV::Grid::ctype CT;
    typedef double Real; 
#if USE_YASPGRID
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,CT,Real,porder> FEM;
#else
    typedef Dune::PDELab::PkLocalFiniteElementMap<GV,CT,Real,porder> FEM;
#endif
    //type of ionic model
   
    typedef MitchellShaefferRegularizedModel::IonicModel<GV,Real> IonicModel;
    //typedef LuoRudy91Model::IonicModel<GV,Real> IonicModel;
   
    enum {ncomp = IonicModel::ncomp};
    // choose ODE solver to solve the parabolic and gating variables
    // Others are ROS2, ROS3p, ROS3PL (Rosenbrock type methods)
    typedef ROS3PL<GV> ODESolver;
    using ES = Dune::PDELab::OverlappingEntitySet<GV>;
    using CON = Dune::PDELab::OverlappingConformingDirichletConstraints;
    
    // now setup the Grid Function Space for coupled parabolic and Gating variables
    typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed> VBE;
  //typedef Dune::LocalSparseMatrix<Real,typename IonicModel::SpeciesList, typename IonicModel::SpeciesInteractionList> MatrixBlock;
    typedef Dune::FieldMatrix<double,ncomp,ncomp> MatrixBlock;
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    typedef Dune::PDELab::GridFunctionSpace<ES,FEM,CON,Dune::PDELab::ISTL::VectorBackend<>> GFS0;
    typedef typename Dune::PDELab::EntityBlockedOrderingTag OrderTag;
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS0,ncomp,VBE,OrderTag> VGFS; 
   
    // setup GFS for extra cellular space
    typedef Dune::PDELab::ISTL::VectorBackend<> VBE1; 
    typedef Dune::PDELab::GridFunctionSpace<ES,FEM,CON,Dune::PDELab::ISTL::VectorBackend<>> UGFS;
   
    // Now setup the GFS for conductivity tensors
    typedef Dune::PDELab::P0LocalFiniteElementMap<CT,Real,dim> P0FEM;
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,Dune::PDELab::NoConstraints,VBE1> P0GFS;
    typedef Dune::PDELab::PowerGridFunctionSpace<P0GFS,dim,VBE1,OrderTag> TGFS;
    using U0 = Dune::PDELab::Backend::Vector<TGFS,Real>;

    // initialization of solution vectors, for coupled Vm and gating variables
    using V = Dune::PDELab::Backend::Vector<VGFS,Real>;

    // extracellular potential 
    using U = Dune::PDELab::Backend::Vector<UGFS,Real>;

    // Grid operator
    typedef typename UGFS::template ConstraintsContainer<Real>::Type UCC;
    typedef typename VGFS::template ConstraintsContainer<Real>::Type VCC;
       
    // choose the linear solver
    //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_Jacobi<VGFS> LinearSolver;
    typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<VGFS,VCC> LinearSolver;	
    
    using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
    typedef Dune::PDELab::GridFunctionSubSpace<VGFS,Path0> GFSV;
    typedef Dune::PDELab::DiscreteGridFunction<GFSV,V> DGFSV;

    // // create objects
    IonicModel ionicmodel;
    FEM fem;
    ES es;
    CON vcon;
    CON ucon; 
    GFS0 gfs0; 
    VGFS vgfs;
    UGFS ugfs;
    V vinit;
    U uinit;
    
    P0FEM p0fem;
    P0GFS p0gfs;
    TGFS tgfs;
    U0 sigma;
 
    BiDHelper(const GV& gv_)
        : gv(gv_), fem(gv_), ionicmodel(), es(gv_), vcon(), ucon(), gfs0(es,fem,vcon), vgfs(gfs0), ugfs(es,fem,ucon), vinit(vgfs,0.0), uinit(ugfs,0.0), p0fem(gv_.template begin<0>()->type()), p0gfs(gv_,p0fem),tgfs(p0gfs),sigma(tgfs,0.0)
    
    {
      if(gv.comm().rank()==0)
        std::cout << "BidHelper: Grid size: elements =" << gv.size(0)<<" on rank "<<gv.comm().rank()<<std::endl;
        // initialize the solution according to Ionic model called
          ionicmodel.InitialSolution(gv,vgfs,vinit);
   /*      using Dune::PDELab::Backend::native;
         std::string sfile = ConfigParser::get<std::string>("Data.finsoultion");
         std::stringstream strfilebuff;
         strfilebuff<<sfile<<"_"<<gv.comm().rank()<<".txt";
         std::string fname = strfilebuff.str();
         readVWDatafromFile(gv,native(vinit),fname);
*/
    }

    void preprocess()
    {
        // initial and end simulation time
        t0 = 0.0;
        tn_dns =   ConfigParser::get("global.tenddns",6.0);
        tau = ConfigParser::get("global.tau",4.0E-2);
            
	// Liinear solver components
        maxlineariters =  ConfigParser::get("LinearSolver.maxiters",500);; 
        reduction       = ConfigParser::get("LinearSolver.reduction",2.0E-6);

        // verbosity for solver
        verbose_primal = ConfigParser::get("LinearSolver.verbose_primal",0);
               
        // flags for different vtkoutput writers
        vtkoutprimal = ConfigParser::get("vtkout.primal",false);
        vtkoutinterface = ConfigParser::get("vtkout.interface",false);
        modulo = ConfigParser::get("vtkout.modulo",3);
        vtkcomps = ConfigParser::get("vtkout.vtkcomps",1);
   
        // filename for vtkout
        vtkoutpath = ConfigParser::get<std::string>("vtkout.vtkoutpath");
        fnprimal = ConfigParser::get<std::string>("vtkout.fnprimal");
            
        
        // flags to check the first and second derivatives
        cputimings = ConfigParser::get("debug.cputimings",false);

        // compute the total GFS size;
        gfssize = compute_internalnodes(gv);

             
    }
 
    void postprocess(){}
	

    ~BiDHelper(){}

public:
    GV gv;
    Real t0, tn_dns, tau,reduction;
    bool vtkoutprimal;
    bool vtkoutinterface;
    std::string fnprimal, sfname, vtkoutpath;
    bool cputimings;
    unsigned modulo, vtkcomps, verbose_primal, maxlineariters;
    int gfssize;
};



template<class GV, class P0GFS, class TensorData>
class ElectroMechanicalSimulation
{

    //define the HELPER class
    typedef BiDHelper<GV> BidHelper;
    typedef typename BidHelper::Real Real;
    BidHelper bidhelper;
    enum{dim=GV::dimension, ncomp = BidHelper::ncomp, qorder=BidHelper::qorder,isParallel=!Dune::MPIHelper::isFake};
    
    typedef typename BidHelper::FEM FEM;
    typedef typename BidHelper::CON   CON;
    typedef typename BidHelper::IonicModel IonicModel;
    typedef typename BidHelper::VGFS   VGFS;
    typedef typename BidHelper::UGFS   UGFS;
    // typedef typename BidHelper::P0GFS   P0GFS;
    typedef typename BidHelper::TGFS   TGFS; 
    typedef typename BidHelper::V   V;
    typedef typename BidHelper::U   U;
    typedef typename BidHelper::U0   U0;
    typedef typename BidHelper::GFSV GFSV;
    typedef typename BidHelper::DGFSV DGFSV;
    typedef typename BidHelper::VCC VCC;
    typedef typename BidHelper::MBE MBE;

    // Define the finite elasticity model
    typedef FiniteElasticity<GV,IonicModel> FiEl;
    FiEl fiel;
    typedef typename FiEl::GFSUP GFSUP;
    typedef typename FiEl::UP UP;
    typedef typename FiEl::MatrixBackend MatrixBackend;

    // Make grid operator space; Matrix Type is same for all 
    typedef Dune::PDELab::BiDomain::MassMatrixLocalOperator<IonicModel,ncomp> MLOP; 
    typedef Dune::PDELab::GridOperator<VGFS,VGFS,MLOP,MBE,Real,Real,Real,VCC,VCC> MGOS;
    typedef typename MGOS::Jacobian MatrixType; // MatrixType is same for all

    using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
    using Path1 = Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>;
    typedef Dune::PDELab::GridFunctionSubSpace<GFSUP,Path0> DisplacementSubGFS;
    typedef Dune::PDELab::VectorDiscreteGridFunctionGradient<DisplacementSubGFS,UP> VDGFGradU;
    typedef typename GFSUP::template ConstraintsContainer<Real>::Type CCUP;   
    typedef Dune::Heart::Elasticity::Scalar_BC<GV,FiEl,Real> BC;
  
public:
    ElectroMechanicalSimulation(const GV& gv_): gv(gv_), bidhelper(gv_), fiel(gv_,bidhelper.ionicmodel)
    {
      if(gv.comm().rank()==0)
        std::cout << "BidSim: Con: Grid size: elements =" << gv.size(0)<<" on rank " << gv.comm().rank()<<" size "<<gv.comm().size()<<std::endl;
    }

    // instantiate the preprocess if needed
    void preprocess()
    {
        // now call the helper class proprocess
        bidhelper.preprocess();
    }

    // now run the main optimization algorithm
    void run(const P0GFS& p0gfs, const TensorData& sigma, Dune::GeometryType::BasicType elemtype)
    {
      if(gv.comm().rank()==0)
         std::cout << " BidSim:run: Grid size: elements =" << gv.size(0)<<" on rank " << gv.comm().rank()<<" size "<<gv.comm().size()<<std::endl;
        using Dune::PDELab::Backend::native;
        VCC vcc;
        vcc.clear();
        NeumannBC<GV,ncomp> vnbc(gv);     
        Dune::PDELab::constraints(vnbc,bidhelper.vgfs,vcc);
        
        if(gv.comm().rank()==0)
            std::cout << "constrained dofs for Dirichlet BC =" << vcc.size() 
                      << " of " << bidhelper.vgfs.globalSize() << std::endl; 

	// making the cpnstraints map and initilizing it from the function in the finite elasticity model
        CCUP ccup;
        ccup.clear();
  	BC bux(gv,fiel,0), buy(gv,fiel,1); // velocity in x and y directions
	BC bp(gv,fiel);                      // pressure
	auto bu = Dune::PDELab::CompositeConstraintsParameters<BC,BC>(bux,buy);	
	auto bup = Dune::PDELab::CompositeConstraintsParameters<decltype(bu),BC>(bu,bp);
        Dune::PDELab::constraints(bup,fiel.gfsup,ccup); // assemble constraints
/*	
        // Define boundary constraints for u and p
        auto buxlambda = [](const auto& x){ // Dirichlet for x-component of velocity at y = 0
          if (x[1]<1e-5) return true;//if (x[1]<1e-5) return true;
          return false;
        };
        auto bux = Dune::PDELab::makeBoundaryConditionFromCallable(gv,buxlambda);
        auto buylambda = [](const auto& x){ // Dirichlet for y-component of velocity at y = 0
          if (x[1]<1e-5) return true;//if (x[1]<1e-5) return true;
          return false;
        };

        auto buy = Dune::PDELab::makeBoundaryConditionFromCallable(gv,buylambda);
        auto bu = Dune::PDELab::CompositeConstraintsParameters<decltype(bux),decltype(buy)>(bux,buy);
        auto bplambda = [](const auto& x){return false;};
        auto bp = Dune::PDELab::makeBoundaryConditionFromCallable(gv,bplambda);
        auto bup = Dune::PDELab::CompositeConstraintsParameters<decltype(bu),decltype(bp)>(bu,bp);
        Dune::PDELab::constraints(bup,fiel.gfsup,ccup); // assemble constraints
*/
        if(gv.comm().rank()==0)
            std::cout << "elasticity constrained dofs=" << ccup.size() 
                 << " from " << fiel.gfsup.globalSize() << std::endl;

        UP up(fiel.up);
        DisplacementSubGFS displacementsubgfs(fiel.gfsup); // subspace
        VDGFGradU vdgfGradU(displacementsubgfs,up); // current displacement gradient as a grid function
        // evaluate Deformation gradient F
        auto getF = [&](const auto& e, const auto& x ){
            Dune::FieldMatrix<Real,dim,dim> s;
            vdgfGradU.evaluate(e,x,s);
            for(std::size_t i=0; i<dim; ++i)
                s[i][i] += 1.0;
            return s;
        };

        V v(bidhelper.vinit);
        GFSV gfsv(bidhelper.vgfs);
        DGFSV dgfsv(gfsv,v);
        // evaluate transmembrane potentional 
        auto getv = [&](const auto& e, const auto& x ){
            Dune::FieldVector<Real,1> s;
            dgfsv.evaluate(e,x,s);
            return s;
        };

        Dune::Timer watch;   
        Dune::PDELab::set_constrained_dofs(vcc,0.0,bidhelper.vinit); 
        MBE mbe(dim==2 ? 9:27);// Maximal number of nonzeroes per row can be cross-checked by patternStatistics().   
        // Make mass matrix operator
        MLOP mlop(bidhelper.ionicmodel,qorder); 
        MGOS mgos(bidhelper.vgfs,vcc,bidhelper.vgfs,vcc,mlop,mbe); 
        MatrixType m(mgos); m = 0.0;   
        mgos.jacobian(bidhelper.vinit,m);
        //if(gv.comm().rank()==0)
          std::cout << "Mass matrix assemble time: =" << watch.elapsed()<<std::endl;     
        
        // Make  stiffness operator for intracellular space
        typedef Dune::PDELab::BiDomain::StiffnessMatrixLocalOperator_mechano_intra<TGFS,U0,IonicModel,decltype(getF),ncomp> ALOP; 
        ALOP alop(bidhelper.tgfs,bidhelper.sigma,bidhelper.ionicmodel,getF,qorder);  //true 
        typedef Dune::PDELab::GridOperator<VGFS,VGFS,ALOP,MBE,Real,Real,Real,VCC,VCC> AGOS;
        AGOS agos(bidhelper.vgfs,vcc,bidhelper.vgfs,vcc,alop,mbe);
        MatrixType ai(agos); ai = 0.0;
        agos.jacobian(bidhelper.vinit,ai); 

        if(gv.comm().rank()==0)
        std::cout<<" cpu times for assembling mass and stiffmatrices "<<watch.elapsed()<<std::endl;

        if(gv.comm().rank()==0)
          std::cout<<"frob norm "<< native(m).frobenius_norm()<<"  "<<native(ai).frobenius_norm2()<<" "<<bidhelper.vinit.two_norm()<<std::endl;
        if(gv.comm().rank()==0)
            std::cout<<" Run the bidomain program "<<std::endl;
        // now define the classes
    	ElectroMechanicalSolver<GV,P0GFS,TensorData,VCC,BidHelper,CCUP,FiEl,decltype(getv)> emsolve(gv,p0gfs,sigma,
                vcc,bidhelper,ccup,fiel,getv,qorder);       
   
        emsolve.preprocess();
        emsolve.run(v,up,m,ai,agos,mgos,bidhelper.cputimings,bidhelper.vtkoutprimal);
    }

    // post process
    void postprocess()
    {
    }
    
    // destructor
    ~ElectroMechanicalSimulation()
    {
    }
public:
    GV gv;
};


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
        //initialize Mpi
        Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
        if(Dune::MPIHelper::isFake)
            std::cout<< "This is a sequential program." << std::endl;
        else
            {
                if(helper.rank()==0)
                    std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
            }
  
        using namespace Dune;
    
        // the GridSelector :: GridType is defined in gridtype.hh and is 
        // set during compilation
        
       
        ConfigParser::parse((argc >= 2 ? argv[1] : "parameters.ini"));
        BidConfigSet::parse((argc >= 2 ? argv[1] : "parameters.ini"));//just for writing the dataset in that folder
                                                                     //should be avoided double instantiation!!!!!

#if USE_YASPGRID
        // create the grid 
        const int dim = 2;
       
        // initialize the cells and size of the grid
        std::vector<int> nel(2);   // number of grid pints in x and y directions
        std::vector<double> L(2);  // domain size considering [0,L0]x[0,L1]
        std::array<int, 2> Par;
        nel[0] =   ConfigParser::get("global.xcells",69);
        nel[1] =   ConfigParser::get("global.ycells",69);
        L[0] =  ConfigParser::get("global.xmax",1.0);
        L[1] =  ConfigParser::get("global.ymax",1.0);
        int Partitioningx = ConfigParser::get("global.xprocs",1);
        int Partitioningy = ConfigParser::get("global.yprocs",1);
        int refinements = ConfigParser::get("global.refinements",0);
        Par[0] = Partitioningx; Par[1] = Partitioningy;
	
        Dune::Timer watch;
        std::vector<double> times(3);

        // ========= Setup YaspGrid ======== //
        watch.reset();
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
        typedef Dune::YaspGrid<dim,Dune::TensorProductCoordinates<double,dim>> HOSTGRID;
        YaspPartition<dim,std::array<int,dim>> yp(Par);
        HOSTGRID hgrid(coords,periodic,overlap,helper.getCommunicator(),(Dune::YLoadBalance<2>*)&yp);
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
        if(refinements > 0) hgrid.globalRefine(refinements);
        typedef HOSTGRID::LeafGridView YGV;
        const YGV ygv = hgrid.leafGridView();
        //typedef int U0; U0 sigma = 0; // dummy sigma parameter for YASP grid construction
        //typedef int P0GFS; P0GFS p0gfs = 0; // dummy sigma parameter for YASP grid construction
        
        if(helper.rank() == 0){
              Dune::gridinfo(hgrid);
              std::cout << "Number of elements per processor: " << ygv.size(0) << std::endl;
              std::cout << "Number of nodes per processor: "    << ygv.size(2) << std::endl;
        }
#else
        // create the grid 
        const int dim = 2;
        using HOSTGRID = Dune::ALUGrid<dim,dim,Dune::simplex,Dune::conforming>;
        std::string  dgfFileName = ConfigParser::get<std::string>("Grid.filename");
        // // // create grid pointer
        Dune::GridPtr<HOSTGRID> gridPtr( dgfFileName);
        
        if(helper.rank()==0)
            Dune::gridinfo((*gridPtr));
#endif

        // ========= Grid Transformation ======== //
        typedef GridTransformation<dim> GRID_TRAFO;
        GRID_TRAFO gTrafo(helper.rank());

        typedef typename Dune::GeometryGrid<HOSTGRID,GRID_TRAFO> GRID;
#if USE_YASPGRID        
        GRID grid(hgrid,gTrafo);
#else
        GRID grid(*gridPtr,gTrafo);
#endif 
        if(helper.rank() == 0)
            std::cout << "Grid transformation complete" << std::endl;

        //Define Grid view: geometry grid wrapper Grid View
        typedef typename GRID::LeafGridView GV;
        const GV gv = grid.leafGridView();
        if(helper.rank() == 0)
            std::cout << "Grid view set up" << std::endl;

        typedef typename GRID::ctype CT;
        typedef Dune::PDELab::P0LocalFiniteElementMap<CT,double,dim> P0FEM;
        typedef typename Dune::PDELab::EntityBlockedOrderingTag OrderingTag;
        typedef Dune::PDELab::VectorGridFunctionSpace<GV,P0FEM,dim,Dune::PDELab::ISTL::VectorBackend<>,Dune::PDELab::ISTL::VectorBackend<>,Dune::PDELab::NoConstraints,OrderingTag> P0GFS;
        //typedef Dune::PDELab::GridFunctionSpace<HGV,P0FEM,Dune::PDELab::NoConstraints,VBE1> P0GFS;
        //typedef Dune::PDELab::PowerGridFunctionSpace<P0GFS,dim,VBE1,OrderTag> TGFS;
        using U0 = Dune::PDELab::Backend::Vector<P0GFS,double>;
        // fine grid part
        P0FEM p0fem(gv.template begin<0>()->type()); //Dune::GeometryTypes::simplex(dim));
        P0GFS p0gfs(gv,p0fem);
        p0gfs.name("displacement");

        U0 sigma(p0gfs,0.0);
        // now create the Tensor data 
#if USE_ALUGRID 
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
        

        
        //read the conductivity tensor data before release the grid pointer 
        readConductivityTensorsFromGrid((*gridPtr),gridPtr,Dune::PDELab::Backend::native(sigma));
        std::cout<<" sigma sizes "<<Dune::PDELab::Backend::native(sigma).size()<<" "<<Dune::PDELab::Backend::native(sigma)[0].size()<<std::endl;
        //vtkout_tensors((*gridPtr).leafGridView(),p0gfs,sigma,"tensors-1");exit(1);
        //std::shared_ptr<HOSTGRID> hgrid( gridPtr.release() ); // finally the grid pointer and store in shared pointer   
#endif     
        // instantiate the simulation object
        ElectroMechanicalSimulation<GV, P0GFS, U0> emsimulation(gv);	

        // run preprocess; load the parameters and some other stuff
        emsimulation.preprocess();
        
        // now run the main optimization
        emsimulation.run(p0gfs,sigma, elemtype);

        // run the post process
        emsimulation.postprocess();
        
        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
