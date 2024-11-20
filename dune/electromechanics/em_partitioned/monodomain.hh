// Primal solve of the bidomain equations
template<class GV, typename VCC, typename BidHelper>
class MonodomainSolver
{

    enum{dim=GV::dimension, ncomp = BidHelper::ncomp,  isParallel=BidHelper::isParallel};
    typedef typename BidHelper::Real Real;
    typedef typename BidHelper::VGFS VGFS;
    typedef typename BidHelper::UGFS UGFS;
    typedef typename BidHelper::MBE MBE; 
    typedef typename BidHelper::V V; 
    typedef typename BidHelper::U U; 
    typedef typename BidHelper::IonicModel IonicModel;
    typedef typename BidHelper::ODESolver ODESolver;
    typedef typename BidHelper::LinearSolver LinearSolver;
    //    typedef typename BidHelper::LinearSolver_e LinearSolver_e;
   
    // Rosenbrock Jacobian matrix type
  typedef Dune::FieldMatrix<bool,ncomp,ncomp> FieldBoolMatrix;
  typedef Dune::PDELab::BiDomain::RBBidJacobianOperator<IonicModel,ODESolver,FieldBoolMatrix,ncomp> LOP; 
    typedef Dune::PDELab::GridOperator<VGFS,VGFS,LOP,MBE,Real,Real,Real,VCC,VCC> GOS;
   
public:
    MonodomainSolver(const GV& gv_, VCC& vcc_, BidHelper& bidhelper_, int qorder_=2, int verbose_=0):gv(gv_), vcc(vcc_),bidhelper(bidhelper_), qorder(qorder_), verbose(verbose_)
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
    bool run(MatrixType& m, MatrixType& ai, AGO& ago, MGOS& mgos, bool debug = false, bool flag_vtkwriter=false)
    {
        Dune::Timer watch, watch2;
        watch2.reset();
	

        LinearSolver linearsolver(bidhelper.vgfs,bidhelper.maxlineariters,bidhelper.verbose_primal);
        //LinearSolver linearsolver(bidhelper.maxlineariters,bidhelper.verbose_primal);

	using Dune::PDELab::Backend::native;
	
        int gfssize = bidhelper.gfssize;
        Real maxtau = bidhelper.ionicmodel.maxtau;
        ODESolver odesolver(gv,maxtau,gfssize);
        odesolver.reduction = bidhelper.reduction;

	
        // Local operator
        FieldBoolMatrix bool_matrix(0.0);
        //CreateBooleanMatrixFromSpeciesList(bidhelper.ionicmodel,bool_matrix);
        LOP lop(bidhelper.ionicmodel,odesolver,bool_matrix,qorder);
        MBE mbe(5);
        GOS gos(bidhelper.vgfs,vcc,bidhelper.vgfs,vcc,lop,mbe);
        typename BidHelper::V v(bidhelper.vinit), f(bidhelper.vinit);
        typename BidHelper::V Itr(bidhelper.vinit),  Istim_supp(bidhelper.vinit);
        typename BidHelper::U u(bidhelper.uinit), g(bidhelper.uinit);
        u = 0.0;
        //bidhelper.ionicmodel.InitialStimulus(gv,bidhelper.vgfs,Itr1,Itr2);
	//Now assemble rhs as const value
        // using stimulus-1 (S1)

        double t0 = bidhelper.t0;
        double tend = bidhelper.tn_dns;
        double tau = bidhelper.tau;

         bidhelper.ionicmodel.setTime(0);
         bidhelper.ionicmodel.setTimeStep(tau);
        Istim_supp = 0.0;
        mgos.residual(v,Istim_supp);

	  
       	double stimulus1_time = bidhelper.ionicmodel.S1_time; 
        double stimulus2_time = bidhelper.ionicmodel.S2_time; 
        double stimulus_delta = bidhelper.ionicmodel.S1S2_delta;
        unsigned counter = 0;
	   
        if(gv.comm().rank()==0)
            std::cout<<" Before VTKwrite in bidomain"<<std::endl;
        // vTK output
        int modulo = bidhelper.modulo;
        VTKOutput<GV,BidHelper>* vtkwriter = (flag_vtkwriter ? new VTKOutput<GV,BidHelper>(gv,bidhelper):NULL);
        if(vtkwriter!=NULL) 
            {
	      (*vtkwriter).write(v,t0);
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
		    //bidhelper.ionicmodel.setIstim(0);
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
                bool adaptive = true; 
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
                            std::cout<<counter<<"  time : "<<t0<<"  oldstep:"<<oldtau<<"   newstep:"<<tau<<"  norms: "<<v.two_norm()<<" "<<Itr.two_norm()<<std::endl;
                        // increment the counter and dump the output
                        if(counter%modulo == 0)
                            if(vtkwriter!=NULL) 
                                {
                                    (*vtkwriter).write(v,t0);
                                }
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
    VCC& vcc;
    const BidHelper& bidhelper;
    int verbose, qorder;
    double parabolic_solveriterations, elliptic_solveriterations, assembleintra_cputime;
    double parabolic_cputime, elliptic_cputime, matass_cputime, primalsolve_cputime;
};
