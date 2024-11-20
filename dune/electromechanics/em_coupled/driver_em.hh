//===============================================================
// a driver for solving the Elastomechanical problem.
//===============================================================

template<typename GV, typename Scheme>
void driver_em (const GV& gv, const Scheme& scheme, Dune::MPIHelper& helper)
{
  // constants and types
  const int dim = GV::dimension;
  typedef double RF;                   // type for computations
  const int degreeu = Scheme::degreeu;
  const int degreep = Scheme::degreep;
  const int degreev = Scheme::degreev;
  const int degreew = Scheme::degreew;
  const int nu = Scheme::nu;
  const int np = Scheme::np;
  const int nv = Scheme::nv;
  const int nw = Scheme::nw;
  const int mvol = Scheme::mvol;
  const int mbnd = Scheme::mbnd;
  const int faces = Scheme::faces;

  // Printing info
  if(helper.rank() == 0){
      std::cout <<"Mechanical model\n";
      std::cout << "dim=" << dim << " degreeu=" << degreeu
	             << " degreep=" << degreep << " nu=" << nu
	             << " np=" << np << std::endl;  
      std::cout <<"\nElectrical model\n";
      std::cout << "dim=" << dim << " degreev=" << degreev
	             << " degreew=" << degreew << " nv=" << nv
	             << " nw=" << nw << std::endl;  
  }
  using Dune::PDELab::Backend::native;

  typedef MitchellShaefferRegularizedModel::IonicModel<GV,RF> IonicModel;
  IonicModel ionicmodel; // The code is writtem for ncomp =1 only
  enum {ncomp = IonicModel::ncomp};
  typedef Dune::FieldMatrix<double,6,6> MatrixBlock;
  //using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<MatrixBlock>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  

  // Define the initial/boundary condition for w (w0) gating variable
  auto wlambda = [](const auto& e, const auto& x ){
    auto global = e.geometry().global(x);
    RF s = 1.0;
    return s;
  };
  auto wf = Dune::PDELab::makeGridFunctionFromCallable(gv,wlambda);
  // Define the initial/boundary condition for v (v0) transmembrane potential
  auto vlambda = [](const auto& e, const auto& x ){
    auto global = e.geometry().global(x);
    RF s = -80;//200.0/(1+exp(20.0*sqrt(global[0]*global[0]+(global[1]-0.5)*(global[1]-0.5))))-84.0;
    //RF s = 1.0-(1.0+exp(-50.0*sqrt(global[0]*global[0]+(global[1]-0.5)*(global[1]-0.5))));
    //RF s = 1.0/(1+0.0001*exp(sqrt(0.5)*(global[0])));
    return s;
  };
  auto vf = Dune::PDELab::makeGridFunctionFromCallable(gv,vlambda);

  // Define boundary conditions for d and p displacement and pressure
  auto ulambda = [&](const auto& x){
      Dune::FieldVector<RF,dim> rv(0.0);
      return rv;
  };
  auto uf = Dune::PDELab::makeGridFunctionFromCallable(gv,ulambda);
  auto plambda = [&](const auto& x){
      RF s = 0.0; return s;
  };
  auto pf = Dune::PDELab::makeGridFunctionFromCallable(gv,plambda);
  auto upvwf = Dune::PDELab::CompositeGridFunction<
  decltype(uf),decltype(pf),decltype(vf),decltype(wf)>(uf,pf,vf,wf);

  // Define boundary constraints for w (Homogeneous Neumann boundary conditions)
  auto bwlambda = [](const auto& x){return false;};
  auto bw = Dune::PDELab::makeBoundaryConditionFromCallable(gv,bwlambda);
  auto bvlambda = [](const auto& x){return false;};
  auto bv = Dune::PDELab::makeBoundaryConditionFromCallable(gv,bvlambda);
  
  // Define boundary constraints for u and p
    auto buxlambda = [&](const auto& x){ // Dirichlet for x-component of velocity at y = 0
      if (x[1]<1e-5) return true;
      return false;
    };
    auto bux = Dune::PDELab::makeBoundaryConditionFromCallable(gv,buxlambda);
    auto buylambda = [&](const auto& x){ // Dirichlet for y-component of velocity at y = 0
      if (x[1]<1e-5) return true;
      return false;
    };
    auto buy = Dune::PDELab::makeBoundaryConditionFromCallable(gv,buylambda);
    auto bu = Dune::PDELab::CompositeConstraintsParameters<decltype(bux),decltype(buy)>(bux,buy);
    auto bplambda = [](const auto& x){return false;};
    auto bp = Dune::PDELab::makeBoundaryConditionFromCallable(gv,bplambda);
    auto bupvw = Dune::PDELab::CompositeConstraintsParameters<
          decltype(bu),decltype(bp),decltype(bv),decltype(bw)>(bu,bp,bv,bw);

  // grid function space for gating variable
  typedef typename Scheme::FEMW FEMW;
  typedef Dune::PDELab::ISTL::VectorBackend<> VectorBackend;
//  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMW,CON,VectorBackend> GFSW;
  GFSW gfsw(gv,scheme.femw);
  gfsw.name("gating variable");

  // grid function space for transmembrane potential
  typedef typename Scheme::FEMV FEMV;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMV,CON,VectorBackend> GFSV;
  GFSV gfsv(gv,scheme.femv);
  gfsv.name("transmembrane potential");

  // grid function space for displacement
  typedef typename Scheme::FEMU FEMU;
  typedef Dune::PDELab::VectorGridFunctionSpace
    <GV,FEMU,dim,VectorBackend,VectorBackend,CON> GFSU;
  GFSU gfsu(gv,scheme.femu);
  gfsu.name("displacement");

  // grid function space for pressure
  typedef typename Scheme::FEMP FEMP;
  typedef Dune::PDELab::GridFunctionSpace
    <GV,FEMP,CON,VectorBackend> GFSP;
  GFSP gfsp(gv,scheme.femp);
  gfsp.name("pressure");

  // combined grid function space
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VectorBackend,
    Dune::PDELab::LexicographicOrderingTag,
    GFSU,GFSP,GFSV,GFSW> GFSUPVW;
  GFSUPVW gfsupvw(gfsu,gfsp,gfsv,gfsw);
  gfsupvw.name("coupled-upvw");

  // assemble the constraints using the given constraints function
  typedef typename GFSUPVW::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bupvw,gfsupvw,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " from " << gfsupvw.globalSize() << std::endl;

  // set up the coefficient vector
  using UPVW = Dune::PDELab::Backend::Vector<GFSUPVW,RF>;
  UPVW upvw(gfsupvw);

  using Dune::PDELab::Backend::native;
  // // initialize coefficient vector
  // bool readFromFile = ConfigParser::get("Data.readfromfile", (bool) false);
  // if ( readFromFile)   
  // {
  //    // INPUT data file
  //    std::string sfname = ConfigParser::get<std::string>("Data.readsolution");
  //    std::stringstream strfilebuff;
  //    strfilebuff<<sfname<<"/txt/"<<sfname<<helper.rank()<<".txt";
  //    std::string fname = strfilebuff.str();
  //    readUPVWDataFromFile(gv,native(upvw),fname);
  // }
  // else
     Dune::PDELab::interpolate(upvwf, gfsupvw,upvw);


  using LOP = ElectromechanicalLOP<Scheme,IonicModel,ncomp>;
   LOP lop(scheme,ionicmodel,gv);
  using TLOP = ElectromechanicsMassLOP<FEMV>;
  TLOP tlop;
  MBE mbe(66);
  using GO0 = Dune::PDELab::GridOperator<GFSUPVW,GFSUPVW,LOP,MBE,RF,RF,RF,CC,CC>;
  GO0 go0(gfsupvw,cc,gfsupvw,cc,lop,mbe);
  using GO1 = Dune::PDELab::GridOperator<GFSUPVW,GFSUPVW,TLOP,MBE,RF,RF,RF,CC,CC>;
  GO1 go1(gfsupvw,cc,gfsupvw,cc,tlop,mbe);
  using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1>;
  IGO igo(go0,go1);

  using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU< GFSUPVW, CC >;
  LS ls(gfsupvw,cc, 5000,1);

  // Solver setup for Finite Elasticity
  typedef Dune::PDELab::NewtonMethod<IGO,LS> NewtonEM;
  NewtonEM newtonem(igo,ls);
  newtonem.setMinLinearReduction(5e-2);
  newtonem.setFixedLinearReduction(true);
  newtonem.setReduction(1e-5);
  newtonem.setAbsoluteLimit(1e-6);

  // select and prepare time-stepping scheme
  Dune::PDELab::OneStepThetaParameter<RF> method1(1.0);
  Dune::PDELab::Alexander2Parameter<RF> method2;
  Dune::PDELab::Alexander3Parameter<RF> method3;
  Dune::PDELab::TimeSteppingParameterInterface<RF>* pmethod=&method1;
  // if (torder==1) pmethod = &method1;
  // if (torder==2) pmethod = &method2;
  // if (torder==3) pmethod = &method3;
  Dune::PDELab::OneStepMethod<RF,IGO,NewtonEM,UPVW,UPVW> osmG(*pmethod,igo,newtonem);
  osmG.setVerbosityLevel(1);
  
  // prepare VTK writer
       int subsampling=1;//ptree.get("output.subsampling",(int)1);
       using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
       VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
       std::string filename=ConfigParser::get<std::string>("Data.fnsolution");
       std::stringstream vtksubfolder;
       vtksubfolder<<filename<<"/vtk";
       std::string subfilename = vtksubfolder.str();
       char cwd[PATH_MAX];
       getcwd(cwd, sizeof(cwd));
       struct stat st;
       if( stat( filename.c_str(), &st ) != 0 )
       {
          int stat = 0;
          stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
          stat = mkdir(subfilename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
          vtksubfolder.str("");
	  vtksubfolder<<filename<<"/txt";
          subfilename = vtksubfolder.str();
          stat = mkdir(subfilename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
          if( stat != 0 && stat != -1)
             std::cout << "Error: Cannot create directory "
                  << filename << std::endl;
       }
       using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
       VTKSEQUENCEWRITER vtkSequenceWriter(gv,filename,cwd + ("/" +  filename + "/vtk"),"",Dune::VTK::conforming);
       Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfsupvw,upvw);

   // prepare the OUTPUT data files
       std::string sfname = ConfigParser::get<std::string>("Data.fnsolution");
       std::stringstream strfilebuff;
       strfilebuff<<sfname<<"/txt/"<<sfname<<helper.rank()<<".txt";
       std::string fname = strfilebuff.str();

  // initialize simulation time and write the initial solution to vtk
  RF t0 = ConfigParser::get("Data.tini", (RF) 0.0);
  vtkSequenceWriter.write(t0,Dune::VTK::appendedraw);
  double stimulus1_time = ionicmodel.S1_time; 
  double stimulus2_time = ionicmodel.S2_time; 
  double stimulus_delta = ionicmodel.S1S2_delta;

  // time loop
  RF T = ConfigParser::get("global.tenddns", (RF) 20.0);
  RF dt = ionicmodel.maxtau;
  RF writetxtT = ConfigParser::get("Data.writetxtT", (RF) 20.0);
  int modulo  = ConfigParser::get("vtkout.modulo",3);
  unsigned counter = 0;
  bool writetofile = ConfigParser::get("Data.writetofile", (bool) false);

  Dune::Timer watch;
  watch.reset();
  
  while (t0<T-1e-8)
    {
      t0 += dt;
          
      ionicmodel.setTime(t0);
      ionicmodel.setTimeStep(dt);
      UPVW znew(upvw);
      osmG.apply(t0,dt,upvw,znew);
      upvw = znew;
      counter++;
      
      if(helper.rank() == 0)
        std::cout <<"Time count := "<<counter<<"  solve time: =" << watch.elapsed()<<std::endl;

      // output to VTK file
      if(counter%modulo == 0)
        vtkSequenceWriter.write(t0,Dune::VTK::appendedraw);

      // if(t0 >= writetxtT && writetofile )
      // {
      //   writeInitSol2File(gv,native(upvw),fname);
      //  if(helper.rank() == 0) std::cout<<" Compled and results written to the file\n";
      // }
      
    }
}
