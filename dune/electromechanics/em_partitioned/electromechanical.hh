// Primal solve of the bidomain equations
template<class GV, class P0GFS, class TD, typename VCC, typename BidHelper, typename CCUP, 
     typename FiEl, typename GetV>
class ElectroMechanicalSolver
{

    enum{dim=GV::dimension, ncomp = BidHelper::ncomp,  isParallel=BidHelper::isParallel};
    enum{degreeu=FiEl::degreeu, degreep= FiEl::degreep}; 
    enum{mvol=FiEl::mvol,  mbnd=FiEl::mbnd,faces=FiEl::faces, nu=FiEl::nu, np=FiEl::np};

    typedef typename BidHelper::Real Real;
    typedef typename BidHelper::VGFS VGFS;
    typedef typename BidHelper::UGFS UGFS;
    typedef typename BidHelper::MBE MBE; 
    typedef typename BidHelper::V V; 
    typedef typename BidHelper::U U; 
    typedef typename BidHelper::IonicModel IonicModel;
    typedef typename BidHelper::ODESolver ODESolver;
    typedef typename BidHelper::LinearSolver LinearSolver;

    typedef typename FiEl::FEMU FEMU;
    typedef typename FiEl::FEMP FEMP;
    typedef typename FiEl::GFSUP GFSUP;
    typedef typename FiEl::MatrixBackend MatrixBackend;
    typedef typename FiEl::UP UP;

    //    typedef typename BidHelper::LinearSolver_e LinearSolver_e;
    //  using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
    // typedef Dune::PDELab::GridFunctionSubSpace<VGFS,Path0> GFSV;
    //  typedef Dune::PDELab::DiscreteGridFunction<GFSV,V> DGFSV;

    // Rosenbrock Jacobian matrix type
    typedef Dune::FieldMatrix<bool,ncomp,ncomp> FieldBoolMatrix;
    typedef Dune::PDELab::BiDomain::RBBidJacobianOperator<IonicModel,ODESolver,FieldBoolMatrix,ncomp> LOP; 
    typedef Dune::PDELab::GridOperator<VGFS,VGFS,LOP,MBE,Real,Real,Real,VCC,VCC> GOS;
    // Finite Elasticity operator
    // typedef Dune::PDELab::Elasticity::FiniteElasticityLOP<GetV,FEMU,FEMP,degreeu,
    //                            degreep,nu,np,mvol,mbnd,faces> FELOP;
    typedef Dune::PDELab::Elasticity::FiniteElasticityLOP<FiEl,GetV, P0GFS, TD> FELOP;
    typedef Dune::PDELab::GridOperator<GFSUP,GFSUP,FELOP,MatrixBackend,Real,Real,Real,CCUP,CCUP> FEGO;

public:

    ElectroMechanicalSolver(const GV& gv_, const P0GFS& p0gfs_, const TD& sigma_, VCC& vcc_, BidHelper& bidhelper_, CCUP& ccup_, FiEl& fiel_, GetV& getv_, int qorder_=2, int verbose_=0):
        gv(gv_), p0gfs(p0gfs_), sigma(sigma_), vcc(vcc_), bidhelper(bidhelper_), ccup(ccup_), fiel(fiel_), getv(getv_), qorder(qorder_), verbose(verbose_)
    {
      if(gv.comm().rank()==0)
	std::cout << "BidSolve: Grid size: elements =" << gv.size(0)<<" on rank " << gv.comm().rank()<<" size "<<gv.comm().size()<<std::endl; 
    }

    // prepare the solver part
    void preprocess()
    {
        parabolic_solveriterations = 0.0; 
        parabolic_cputime = 0.0; matass_cputime = 0.0;
        primalsolve_cputime = 0.0;
    }

    // run the main promal solver routine
  template<class MatrixType, typename AGO,typename MGOS>
    bool run(V& v, UP& up, MatrixType& m, MatrixType& ai, AGO& ago, MGOS& mgos, bool debug = false, bool flag_vtkwriter=false)
    {
        using Dune::PDELab::Backend::native;

       // FELOP felop(gv,getv,fiel.femu,fiel.femp);
        FELOP felop(fiel,getv, p0gfs, sigma);
        using FEGO = Dune::PDELab::GridOperator<GFSUP,GFSUP,FELOP,MatrixBackend,Real,Real,Real,CCUP,CCUP>;
        FEGO fego(fiel.gfsup,ccup,fiel.gfsup,ccup,felop,MatrixBackend(66));

        // typedef typename FEGO::Traits::Jacobian K;
        // K k(fego,0.0);
        // fego.jacobian(up,k);
        // writeMatrixToMatlabHelper(native(k), 0, 0, std::cout);
        // exit(0);

        //  using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0;
        //  LS ls(5000,1);

        //  using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR;
        //  LS ls(5000,1);

        // using LS = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<FEGO>;
        // LS ls(fiel.gfsup,5000,2);

        // using LS = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
        // LS ls(5000,1);

        // using LS = Dune::PDELab::ISTLBackend_OVLP_CG_UMFPack<GFSUP, CCUP >;
        // LS ls(fiel.gfsup,ccup, 5000,1);

        // using LS = Dune::PDELab::ISTLBackend_SEQ_SuperLU; LS ls;

        using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU< GFSUP, CCUP >;
        LS ls(fiel.gfsup,ccup, 5000,1);

        //using LS = Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFSUP, CCUP >;
        //LS ls(fiel.gfsup,ccup, 5000,2);

        // using LS = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR< FEGO >;
        // LS ls(5000,2);

        // using LS = Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<FEGO>;
        // LS ls(fiel.gfsup,5000,2);

        // using LS = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<FEGO>;
        // LS ls(5000,1);
        // using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0< GFSUP, CCUP >;
        // LS ls(fiel.gfsup,ccup, 5000,0);

        //nonlinear Solver setup for Finite Elasticity
        typedef Dune::PDELab::NewtonMethod<FEGO,LS> NewtonFE;
        NewtonFE newtonFE(fego,ls);
        newtonFE.setMinLinearReduction(5e-2);
        newtonFE.setFixedLinearReduction(true);
        newtonFE.setReduction(1e-5);
        newtonFE.setAbsoluteLimit(1e-6);

        //newtonFE.setParameters(ptree.sub("newtonFE"));

        Dune::Timer watch, watch2;
        watch2.reset();
        LinearSolver linearsolver(bidhelper.vgfs,vcc,bidhelper.maxlineariters,bidhelper.verbose_primal);
        //LinearSolver linearsolver(bidhelper.maxlineariters,bidhelper.verbose_primal);

	
        int gfssize = bidhelper.gfssize;
        Real maxtau = bidhelper.ionicmodel.maxtau;
        ODESolver odesolver(gv,maxtau,gfssize);
        odesolver.reduction = bidhelper.reduction;

	
        // Local operator
        FieldBoolMatrix bool_matrix(1);
        //CreateBooleanMatrixFromSpeciesList(bidhelper.ionicmodel,bool_matrix);
        LOP lop(bidhelper.ionicmodel,odesolver,bool_matrix,qorder);
        MBE mbe(5);
        GOS gos(bidhelper.vgfs,vcc,bidhelper.vgfs,vcc,lop,mbe);
        typename BidHelper::V f(bidhelper.vinit);
        typename BidHelper::V Itr(bidhelper.vinit),  Istim_supp(bidhelper.vinit);
        typename BidHelper::U u(bidhelper.uinit), g(bidhelper.uinit);
        u = 0.0;
	//Now assemble rhs as const value
        // using stimulus-1 (S1)
        Istim_supp = 0.0;
        mgos.residual(v,Istim_supp);
	
        double t0 = bidhelper.t0;
        double tend = bidhelper.tn_dns;
        double tau = bidhelper.tau;

	    bidhelper.ionicmodel.setTime(0);
        bidhelper.ionicmodel.setTimeStep(tau);

       	double stimulus1_time = bidhelper.ionicmodel.S1_time; 
        double stimulus2_time = bidhelper.ionicmodel.S2_time; 
        double stimulus_delta = bidhelper.ionicmodel.S1S2_delta;
	std::string sfname = ConfigParser::get<std::string>("Data.fnsolution");
        unsigned counter = 0;
	   
        if(gv.comm().rank()==0)
            std::cout<<" Before VTKwrite in bidomain"<<std::endl;
        // vTK output
        int modulo = bidhelper.modulo;
       VTKOutput_mechano<GV,BidHelper,FiEl>* vtkwriter = (flag_vtkwriter ? new VTKOutput_mechano<GV,BidHelper,FiEl>(gv,bidhelper,fiel):NULL);
        if(vtkwriter!=NULL) 
            {
	      (*vtkwriter).write_all(up,v,t0);
            }
	if(gv.comm().rank()==0)
            std::cout<<" After VTKwrite in bidomain"<<std::endl;
        bool stimulus_act = false;
	bool assemble_s2 = true;
        watch2.reset();
        while(t0 <= tend)
            {
                // now solve the parabolic and ODE equations
                watch.reset();
                bidhelper.ionicmodel.setTime(t0);
                bidhelper.ionicmodel.setTimeStep(tau);

                newtonFE.apply(up); // Solving the Mechanical Model
                MatrixType j(gos);
                j = 0.0;
                gos.jacobian(v,j);
                //std::cout<<"frob norm "<< j.frobenius_norm2()<<std::endl;
                //j *= -1.0;
                native(j) += native(ai);
                // complete the unassemble part for LHS matrix for Rosenbrock method
                double dtg = 1.0/(tau*odesolver.gamma[0]);
                // now compute the (1/dtg*M-J) 
                native(j).axpy(dtg,native(m));
                //Dune::printSparseMatrix(std::cout,j.base(),"global stiffness matrix","row",9,1); 
                //if(debug && gv.comm().rank() == 0)
                //std::cout<<"frob norm "<< native(j).frobenius_norm2()<<"  "<<native(m).frobenius_norm2()<<"  "<<native(ai).frobenius_norm2()<<" "<<v.two_norm()<<std::endl;exit(1);
                if(debug)
                    assembleintra_cputime += watch.elapsed();

                // S1-S2 protocal stimulus
                bool stim_activate = false;
                Itr = 0.0;
		
		if(t0 >= (stimulus1_time) && t0 < (stimulus1_time+stimulus_delta) )
		  {
		    stim_activate = true;
		    // now add the excitation stimulus
		    //native(m).mv(native(Itr1),native(Itr));
		    Itr = Istim_supp;

		    // set the Stimulus current(required in Tentusscher model)
		    bidhelper.ionicmodel.setIstim(bidhelper.ionicmodel.Istim);
		  }
                else if( (t0+tau) > (stimulus2_time) && t0 < (stimulus2_time + stimulus_delta))
		  {
		    stim_activate = true;
		    if(assemble_s2)
		      {     
			Istim_supp = 0.0;
			mgos.residual(v,Istim_supp);
                               
			// deactivate the assemble of the second stimulus part
			assemble_s2 = false;
		      }
                       
		    //now copy to old Itr vector
		    Itr = Istim_supp;
		    // set the Stimulus current(required in Tentusscher model)
		    bidhelper.ionicmodel.setIstim(bidhelper.ionicmodel.Istim);
		  }
		else
		  {
		    bidhelper.ionicmodel.setIstim(0);
		  }
                // now finish the rhs by adding extra cellular part
                //f = Itr;
                // for(int i=0;i<u.N();i++)
                //     f.base()[i][0] = u.base()[i];
                // ai.base().usmv(-1.0,f.base(),Itr.base());  
		    
                watch.reset();
                Real err;
                // if post shock simulation is true then activate adaptivity
                //linearsolver.prestage();
                bool adaptive = false; 
                Real oldtau = tau;
                // solve the parabolic and ODEs using RB method
                bool applydualpart = false;
			
                bool accept = odesolver.solve(gos,j,ai,m,v,Itr,linearsolver,odesolver,t0,tau,applydualpart,adaptive);
                if(debug)
                    parabolic_cputime += watch.elapsed();
		
                // now update the time step (by default accept all time steps!!!)
                if(accept)
                    {
                        t0 += oldtau;
                        counter++;
                        if( (t0+tau) > stimulus2_time && !stimulus_act)
                            { 
                                stimulus_act = true;
                                tau = stimulus2_time-t0+1E-2;
                            }	 
                      
                        
                        if(debug && gv.comm().rank() == 0)
                            std::cout<<counter<<"  time : "<<t0<<"  oldstep:"<<oldtau<<"   newstep:"<<tau<<"  norms: "<<v.two_norm()<<" "<<Itr.two_norm()<<std::endl<<std::endl;
                        // increment the counter and dump the output
                        if(counter%modulo == 0)
                            if(vtkwriter!=NULL) 
                                {
                                    (*vtkwriter).write_all(up,v,t0);
                                }
		    //   // Write the solution to file for sequential run
            //           if(counter%100 == 0)
            //           {
            //                std::stringstream strfilebuff;
            //                strfilebuff<<sfname<<"_"<<gv.comm().rank()<<".txt";
            //                std::string fname = strfilebuff.str();
            //                writeData2File(gv,native(up),native(v),fname);
            //             if(debug && gv.comm().rank() == 0)  std::cout<<"At t0 = "<<t0<< " the data is written to file "<<fname<<std::endl;
            //             }

                    }
            }
	    
        if(debug && gv.comm().rank() == 0)
            {
                primalsolve_cputime = watch2.elapsed();
                std::cout
                    <<"      CPU times: assembling intra "<<assembleintra_cputime 
                    <<"  PDE/ODE solver: "<<parabolic_cputime
                    <<"  Primal solve:       "<<primalsolve_cputime<<std::endl;
            }
	   
        //if every thing goes well, return true
        return true;
    }
    
    // now run the DNS for post shock simulation
    // run the main promal solver routine
    bool run_postshock()
    {
	    
        return true;
    }
    
    void postprocess()
    {

    }

public:
    const GV& gv;
    const P0GFS& p0gfs;
    const TD& sigma; // tensor data
    VCC& vcc;
    CCUP& ccup;
    const GetV& getv;
    const BidHelper& bidhelper;
    const FiEl& fiel;
    int verbose, qorder;
    double parabolic_solveriterations, elliptic_solveriterations, assembleintra_cputime;
    double parabolic_cputime, elliptic_cputime, matass_cputime, primalsolve_cputime;
};
