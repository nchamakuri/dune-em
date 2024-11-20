#ifndef DUNE_ELASTICITY_MODEL_HH
#define DUNE_ELASTICITY_MODEL_HH

template<class GV, typename IonicModel = void>
class FiniteElasticity
{
public:
    enum{dim = GV::dimension, isParallel=!Dune::MPIHelper::isFake};
#ifdef CUBE
    enum {degreeu=2}; // displacement polynomial degree
    enum {degreep=1}; // pressure polynomial degree

    enum { mvol=16}; // number of quadrature points on element 
    enum { mbnd=3};  // number of quadrature points on face
    enum { faces=2*dim};  // number of faces

    enum { nu=Dune::QkStuff::QkSize<degreeu,dim>::value};  // # displacement degrees of freedom
    enum { np=Dune::QkStuff::QkSize<degreep,dim>::value};  // # pressure degrees of freedom
#else
    enum {degreeu=2}; // velocity polynomial degree
    enum {degreep=1}; // pressure polynomial degree

    enum { mvol=12}; // number of quadrature points on element 
    enum { mbnd=3};  // number of quadrature points on face
    enum { faces=dim+1};  // number of faces

    enum { nu=Dune::PB::PkSize<degreeu,dim>::value};  // # velocity degrees of freedom
    enum { np=Dune::PB::PkSize<degreep,dim>::value};  // # pressure degrees of freedom
#endif

    typedef typename GV::Grid::ctype CT; // type for coordinates
    typedef double Real;          // representation type for solution values
    typedef Dune::FieldVector<double,dim> FieldVec;

    const GV& gv;
    const IonicModel& ionicmodel;
    // coefficient of elasticity problem
    Real v_min, v_max;         // Ionic Model minimum and maximum potential
    Real mu, E, rho;           // Material Parameters
    Real beta;                 // Parameter used in finding/scaling the active stress 
    const bool active_strain;  // 
    Real gil, git;             // conductivity factors
    
    // # al and at are depend of the position so need to moved into a function
    // # Make this part general for all dimension
    const Dune::FieldMatrix<Real,1,dim> al = {{1.0,0.0}}; // Fiber direction longitudinal
    const Dune::FieldMatrix<Real,1,dim> at = {{0.0,1.0}}; // Fiber direction transverse
    
       
   // finite element map types
#ifdef CUBE
   typedef Dune::PDELab::QkLocalFiniteElementMap<GV,CT,Real,degreeu> FEMU; // finite element map for displacement
   typedef Dune::PDELab::QkLocalFiniteElementMap<GV,CT,Real,degreep> FEMP; // finite element map for pressure
#else
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,CT,Real,degreeu> FEMU; // finite element map for displacement
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,CT,Real,degreep> FEMP; // finite element map for pressure
#endif
   // finite element map objects
   FEMU femu;
   FEMP femp;

  typedef Dune::PDELab::ISTL::VectorBackend<> VectorBackend;
 // typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MatrixBackend;
  // grid function space for displacement
  typedef Dune::PDELab::VectorGridFunctionSpace
    <GV,FEMU,dim,VectorBackend,VectorBackend,CON> GFSU;
    GFSU gfsu;

  // grid function space for pressure
  typedef Dune::PDELab::GridFunctionSpace
    <GV,FEMP,CON,VectorBackend> GFSP;
    GFSP gfsp;

  // combined displacement-pressure grid function space
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VectorBackend,
    Dune::PDELab::LexicographicOrderingTag,
    GFSU,GFSP> GFSUP;
    GFSUP gfsup;
    // set up the coefficient vector
    using UP = Dune::PDELab::Backend::Vector<GFSUP,Real>;
    UP up;
    
    // Constractor for electric-mechanical intraction problem (passive and active stress)
    FiniteElasticity(const GV& gv_, const IonicModel& ionicmodel_): gv(gv_), ionicmodel(ionicmodel_), active_strain(true),
		femu(gv_), femp(gv_), gfsu(gv_,femu), gfsp(gv_,femp), gfsup(gfsu,gfsp), up(gfsup,0.0)
    {
        preprocess();
    }

    // Constructor for only mechanical problem (with out any active stress)
    FiniteElasticity(const GV& gv_): gv(gv_),  active_strain(true),
                femu(gv_), femp(gv_), gfsu(gv_,femu), gfsp(gv_,femp), gfsup(gfsu,gfsp), up(gfsup,0.0)
    {
        preprocess();
    }

    // Function to set the parameters and initilize the variables
    void preprocess()
    {
	// Getting the parameter values from the ini file
	mu =  ConfigParser::get("ElasticityModel.mu",4.0);
	beta = ConfigParser::get("ElasticityModel.beta",0.3);
        rho = ConfigParser::get("ElasticityModel.density", 0.0);
        //read_initial_from_file = ConfigParser::get("FiniteElasticity.read_initial_from_file",false);
        
        // Set the initial/ boundary conditions for d and p displacement and pressure
        auto ulambda = [&](const auto& x){
          Dune::FieldVector<Real,dim> rv(0.0);
          return rv;
        };
        auto uf = Dune::PDELab::makeGridFunctionFromCallable(gv,ulambda);
        auto plambda = [&](const auto& x){
          Real s = 0.0; return s;
        };
        auto pf = Dune::PDELab::makeGridFunctionFromCallable(gv,plambda);
        auto upf = Dune::PDELab::CompositeGridFunction<decltype(uf),decltype(pf)>(uf,pf);
        // initialize coefficient vector
        Dune::PDELab::interpolate(upf, gfsup,up);
      /*  using Dune::PDELab::Backend::native;
        std::string sfile = ConfigParser::get<std::string>("Data.finsoultion");
        std::stringstream strfilebuff;
        strfilebuff<<sfile<<"_"<<gv.comm().rank()<<".txt";
        std::string fname = strfilebuff.str();
        readVPDatafromFile(gv,native(up),fname);
       */ 
        // For electromechanical coupled problem get the ionic model parameters
        if(active_strain)
	{
	    v_min = ionicmodel.ionicpar.v_rest;
	    v_max = ionicmodel.ionicpar.v_peak;
	    gil   = ionicmodel.par.gil;
            git   = ionicmodel.par.git;
	} 
    }

   // Function to get the longitudinal fiber direction at each point
   //
   //
   //
   //
   //
   
    //! A function which returns whether a point x lies on a Dirichlet boundary for the degree of freedom dof (0,1,2) and for pressure if dof is -1.
    //! If the function returns false it assumes this node is on the Neumann boundary or an interior node.
    bool inline isDirichlet(FieldVec& x, int dof)
    {
        if (dof == -1) return false;
	if (x[1]<1e-5) return true;
        return false;
    }
    
    //! returns the Dirichlet boundary value at a point for the pressure or  degree of freedom dof.
    double inline evaluateDirichlet(const FieldVec& x, int dof) const
    {
        return 0.0;
    }

    //! function to evaluate the source function at each point
    inline void evaluateWeight(FieldVec& f, int id) const
    {
        f = 0.0;
        // f[2] = - gravity * rho; 
    }

    //! function which evaluates the Neumann boundary conditions at the point x
    inline void evaluateNeumann(const FieldVec &x, FieldVec& h,const FieldVec& normal) const
    {
        h = 0.0;
    }


    void postprocess(){}
	
    ~FiniteElasticity(){}

};
#endif
