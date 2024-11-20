#ifndef IONICMODEL_MSR_HH
#define IONICMODEL_MSR_HH

#include "ms_regularized_base.hh"

namespace MitchellShaefferRegularizedModel
{

 // for species list implementation
  struct SVm {
    constexpr static const char name[] = "Vm";
    constexpr static const double DefaultDiffusionConstant = 1.0;
  };
    
  struct SGw {
    constexpr static const char name[] = "Gw";
    constexpr static const double DefaultDiffusionConstant = 0;
  };
    
// set the stimulus domain
template<typename GV, typename RF,int m=1>
class  Stimulus1_MSR : public Dune::PDELab::AnalyticGridFunctionBase<
   Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,m>,
   Stimulus1_MSR<GV,RF,m> > {
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,m> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Stimulus1_MSR<GV,RF,m> > B;

   Stimulus1_MSR (const GV& gv, int dimw) : B(gv) 
   {
	  int dim = dimw;
	  if(dim == 2)
		 {zr[0] = 0.5; zr[1] = 0.5; }
	  else
          //{zr[0] = 0.9644; zr[1] = 1.13592; zr[2] = -0.49846;} // TBunnyC2
		 //{zr[0] = 1.93; zr[1] = 0.164; zr[2] = -0.05;} 
          {zr[0] = 0.5; zr[1] = 0.5; zr[2] = 0.5;} //cube geo
   }
   inline void evaluateGlobal (const typename Traits::DomainType& X, 
							   typename Traits::RangeType& Itr) const
   {
	  typename Traits::DomainType z; 
	  z  = zr;
	  z -= X;
	  RF dist = z.two_norm();
	  // if(X[1]<1e-8)
	  //if( dist < 0.25) //for cube 0.08
	  Itr = 0.0;
	  if(GV::dimension==2)
		 {
			if(X[1]<1e-8)  
			   Itr[0] = 30.0;	 
		 }
	  else
		 {
			// if(X[1]>.53 && X[1] <0.80) // for TBunnyC2
			//    if(X[2]<-1.0)
			// 	  Itr[0] = 30.0; //200.0;
			// if(X[1]>-1.02 && X[1] <-0.85) // for TBunnyC2
			//    if(X[2]>-0.1 && X[2]<0.1)
			// 	  Itr[0] = 30.0; //200.0;
			// if( dist < 0.4)
			//    Itr[0] = 20.0; 
         	// if(X[0]>-0.0106 && X[0] <0.065) // for HeartVolume by Alesandro
			//    if(X[2]<-0.02)
			// 	  Itr[0] = 30.0; //200.0;
             // if(X[2]<1e-8)   // for cubic geom
             //     Itr[0] = 30.0;
             if(X[1]<-0.35)   // for Heart geom
                 Itr[0] = 30.0;
		 }
   }

public:  
   typename Traits::DomainType zr; 
};


// set the stimulus domain
template<typename GV, typename RF,int m=1>
class Stimulus2_MSR : public Dune::PDELab::AnalyticGridFunctionBase<
   Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,m>,
   Stimulus2_MSR<GV,RF,m> > {
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,m> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Stimulus2_MSR<GV,RF,m> > B;

   Stimulus2_MSR (const GV& gv, int dimw) : B(gv) 
   {
	  int dim = dimw;
	  if(dim == 2)
		 {zr[0] = 1.55; zr[1] = 1.0; }
	  else
		 //	 {zr[0] = 0.32; zr[1] = 0.21; zr[2] = 1.0;} 
		 //{zr[0] = -0.185; zr[1] = 0.61; zr[2] = -0.55;} // TBUnnyC2 
          {zr[0] = 1.0; zr[1] = 1.0; zr[2] = 1.0;} //cubic geom  
   }
   inline void evaluateGlobal (const typename Traits::DomainType& X, 
							   typename Traits::RangeType& Itr) const
   {
	  typename Traits::DomainType z; 
	  z  = zr;
	  z -= X;
	  RF dist = z.two_norm();
	  //if(X[1]<1e-8)
	  if( dist < 0.3) //for cube 0.08
		 Itr[0] = 20.0;
	  else
		 Itr = 0.0; 
   }

public:  
   typename Traits::DomainType zr; 
};

   template<typename DomainType>
    double S1Support(const DomainType x, const int dim)
    {
        DomainType zr;
        if(dim == 2)
		 {zr[0] = 0.5; zr[1] = 0.5; }
	  else
          //{zr[0] = 0.9644; zr[1] = 1.13592; zr[2] = -0.49846;} // TBunnyC2
		 //{zr[0] = 1.93; zr[1] = 0.164; zr[2] = -0.05;} 
          {zr[0] = 0.5; zr[1] = 0.5; zr[2] = 0.5;} //cube geo
      
        zr -= x;
        double dist = zr.two_norm();
        double Itr = 0.0;
        // if(x[1]<0.021) // for cubic geometris
        //     Itr = 1.0;   
	// if(X[1]>.53 && X[1] <0.80) // for TBunnyC2
			//    if(X[2]<-1.0)
			// 	  Itr[0] = 30.0; //200.0;
			// if(X[1]>-1.02 && X[1] <-0.85) // for TBunnyC2
			//    if(X[2]>-0.1 && X[2]<0.1)
			// 	  Itr[0] = 30.0; //200.0;
			// if( dist < 0.4)
			//    Itr[0] = 20.0; 
         	// if(X[0]>-0.0106 && X[0] <0.065) // for HeartVolume by Alesandro
			//    if(X[2]<-0.02)
			// 	  Itr[0] = 30.0; //200.0;
             // if(X[2]<1e-8)   // for cubic geom
             //     Itr[0] = 30.0;

	if(dim == 2)
	  {
	  if(x[1]<0.15) // for cubic geometris
	    Itr = 1.0;
	  }
	else
	  {
	    ///if(x[1]<-0.35)   // for Heart geom
	    if(x[0] > 1.80)   // for Heart geom: TBunnyC  
	      Itr = 1.0;
	  }
        return Itr;
    }
    
    template<typename DomainType>
    double S2Support(const DomainType x, const int dim)
    {
        DomainType zr;    
        if(dim == 2)
		 {zr[0] = 1.55; zr[1] = 1.0; }
	  else
		 //	 {zr[0] = 0.32; zr[1] = 0.21; zr[2] = 1.0;} 
		 //{zr[0] = -0.185; zr[1] = 0.61; zr[2] = -0.55;} // TBUnnyC2 
	    //{zr[0] = 1.0; zr[1] = 1.0; zr[2] = 1.0;} //cubic geom
	  {zr[0] = 0.0281; zr[1] = 0.5744; zr[2] = 1.1987;} //cubic geom  
        zr -= x;
        double dist = zr.two_norm();
        double Itr = 0.0;

        // if( x[0] >0.1 && x[0] < 0.7) //for cube 0.08
        //     if( x[1] < 1.0) //for cube 0.08   
        //         Itr = 1.0;
        if( dist < 0.5) //for cube 0.08
            Itr = 1.0;
        else
            Itr = 0.0; 
        return Itr;
    }




    
template<typename V>
void initialsolutionMSRmodel(V& y)
{
   y[0] = -80.0;
   y[1] = 1.0;
}

// grid function for initial solution for gating variables
template<typename GV, typename RF,int ncomp=2>
class InitialSolution_MSR : public Dune::PDELab::AnalyticGridFunctionBase<
   Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp>,
   InitialSolution_MSR<GV,RF,ncomp> > {
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,InitialSolution_MSR<GV,RF,ncomp> > B;
   InitialSolution_MSR (const GV& gv) : B(gv) {}
   inline void evaluateGlobal (const typename Traits::DomainType& x, 
							   typename Traits::RangeType& y) const
   {
      // if(x[1]<1E-8)
	  // y[0] = 101.0;
      // else
	  // y[0] = 0.0;
	  // y[1] = 0.0;
	  y= 0.0;
	  initialsolutionMSRmodel(y);
   }
};

#if EXPLICIT_SCHEMES
// grid function for initial solution for Vm
template<typename GV, typename RF,int ncomp=2>
class InitialSolution_MSR_V : public Dune::PDELab::AnalyticGridFunctionBase<
   Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp>,
   InitialSolution_MSR_V<GV,RF,ncomp> > {
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,InitialSolution_MSR_V<GV,RF,ncomp> > B;
   InitialSolution_MSR_V (const GV& gv) : B(gv) {}
   inline void evaluateGlobal (const typename Traits::DomainType& x, 
							   typename Traits::RangeType& y) const
   {
       y = -80.0;
   }
};
// grid function for initial solution for gating variables
template<typename GV, typename RF,int ncomp=2>
class InitialSolution_MSR_W : public Dune::PDELab::AnalyticGridFunctionBase<
   Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp>,
   InitialSolution_MSR_W<GV,RF,ncomp> > {
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,InitialSolution_MSR_W<GV,RF,ncomp> > B;
   InitialSolution_MSR_W (const GV& gv) : B(gv) {}
   inline void evaluateGlobal (const typename Traits::DomainType& x, 
							   typename Traits::RangeType& y) const
   {
       y = 1.0;
   }
};
#endif
// grid function for initial solution for gating variables
template<typename GV, typename RF,int ncomp=2>
class InitialSolution_MSROptim : public Dune::PDELab::AnalyticGridFunctionBase<
   Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp>,
   InitialSolution_MSROptim<GV,RF,ncomp> > {
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,ncomp> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,InitialSolution_MSROptim<GV,RF,ncomp> > B;
   InitialSolution_MSROptim (const GV& gv) : B(gv) 
   {
	  const int dim = GV::dimension;
	  if(dim == 2)
		 {
			zr[0] = 0.5; zr[1] = 0.5;
			//zr[0] = 0.0; zr[1] = 0.0;	
		 }
	  else
		 {
			zr[0] = 0.9644; zr[1] = 1.13592; zr[2] = -0.49846;
		 }
   }
   inline void evaluateGlobal (const typename Traits::DomainType& x, 
							   typename Traits::RangeType& y) const
   { 
      typename Traits::DomainType z; 
      z  = zr;
      z -= x;
      double dist = z.two_norm();
      if( dist < 0.04) //for cube 0.08
		 y[0] = 101.0;
      else
		 y[0] = 0.0;
      y[1] = 0.0;
   }
public:  
   typename Traits::DomainType zr; 
};



// test model problem for 1 component
template<class GV, class RF>
class IonicModel
{
public:
   // for species list implementation
   using SpeciesList = Dune::LSM::Typelist<SVm, SGw>;

   using VmIL = Dune::LSM::Typelist<SVm, SGw>; //1,2
   using GwIL = Dune::LSM::Typelist<SVm, SGw>; //1,2
   // coupled system
  using SpeciesInteractionList = Dune::LSM::Typelist<VmIL, GwIL>;
                                                      
   enum {dim=GV::dimension, dimworld = GV::dimensionworld};
#if EXPLICIT_SCHEMES
  enum {ncomp=1,wncomp=1, cncomp=0};
  typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,RF,ncomp+wncomp> Traits;
  MSRModelParameters<GV,RF,ncomp+wncomp> ionicpar;
#else
    enum {ncomp=2};
  typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,RF,ncomp> Traits;
  MSRModelParameters<GV,RF,ncomp> ionicpar;
#endif
   
   typedef Dune::PDELab::BiDomain::ElectrodeBoundary::Type ElectrodeType;
   BidParameters<GV, RF,ncomp> par;
   
   enum{NODEs = 1};
  
   IonicModel (): par(), ionicpar()
   {
	  for (std::size_t i=0; i<Traits::dimDomain; i++)
		 for (std::size_t j=0; j<Traits::dimDomain; j++)
			sigma_i[i][j] = (1.0/ionicpar.Am)*par.sigma_i[i][j];
	  
	  // time step restriction for the time adaptive computations
      Istim =  ConfigParser::get("MSRModel.Istim",40.0);
      S1_time = ConfigParser::get("MSRModel.S1time",0.0);
      S2_time = ConfigParser::get("MSRModel.S2time",220.0);
      S1S2_delta = ConfigParser::get("MSRModel.S1S2delta",2.0);

      maxtau = ConfigParser::get("MSRModel.maxtau",1.0);
      time_dep = ConfigParser::get("MSRModel.dep_time",20.0);
      maxtau_dep = ConfigParser::get("MSRModel.dep_maxtau",0.25);
   }
  
   //! source/reaction term
   void
   rhs (const typename Traits::DomainType& x, 
		typename Traits::RangeType& v,typename Traits::RangeType& res) const
   {
	  res[0] = -ionicpar.Iion(v)/ionicpar.Cm;
	  res[1] =  ionicpar.Gating(v);
   }

#if EXPLICIT_SCHEMES
    //! source/reaction term
   RF rhs_Iion (const typename Traits::DomainType& x, 
		Dune::FieldVector<double,2> v) const
		//typename Traits::RangeType& v) const
   {
       return -ionicpar.Iion(v)/ionicpar.Cm;
   }
    //! source/reaction term
  std::array<std::array<RF,2>,wncomp> Gating_lhs_rhs (Dune::FieldVector<double,ncomp+wncomp>& v) const //typename Traits::RangeType& v) const
    {
       std::array<std::array<RF,2>,wncomp> res;
       res[0][0] = ionicpar.Gating_lhs(v);
       res[0][1] = ionicpar.Gating_rhs(v);
        return res;
    }

    //! source/reaction term
    RF jacobian_Iion (const typename Traits::DomainType& local, 
		      Dune::FieldVector<double,2>& v) const
		      //typename Traits::RangeType& v) const
   {
	  // derivative of I_ion w.r.t v,w
	  return -ionicpar.Iion_v(v)/ionicpar.Cm;
   }
  //! source/reaction term
    std::array<RF,wncomp> Concentrations_rhs (const Dune::FieldVector<double,ncomp+wncomp+cncomp>& sol) const
   {
       std::array<RF,wncomp> res;
       return res;
   }
  
#endif    
   //! source/reaction term
   void 
   jacobian (const typename Traits::DomainType& local, 
			 typename Traits::RangeType& v, typename Traits::JacobianType& res) const
   {
	  // derivative of I_ion w.r.t v,w
	  res[0][0] = -ionicpar.Iion_v(v)/ionicpar.Cm;
	  res[0][1] = -ionicpar.Iion_w(v)/ionicpar.Cm;

	  // derivative of G(v,w) w.r.t. v,w
	  res[1][0] = ionicpar.Gating_v(v);
	  res[1][1] = ionicpar.Gating_w(v);
   }

   //! diffusion tensor
   void
   Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
	   typename Traits::TensorType& sigma) const
   {
	  for (std::size_t i=0; i<Traits::dimDomain; i++)
		 for (std::size_t j=0; j<Traits::dimDomain; j++)
			sigma[i][j] = sigma_i[i][j];
   }

  //! diffusion tensor
   void
   Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x, const typename Traits::DomainType& al,
	   typename Traits::TensorType& sigma) const
   {
     sigma = 0.0;
     for (int i=0; i<dim; i++)
       for (int j=0; j<dim; j++)
	 if(i==j)
	   sigma[i][j] = par.gil*al[i]*al[j] + (1.0-al[i]*al[j])*par.git;
	 else
	   sigma[i][j] = par.gil*al[i]*al[j] + (0.0-al[i]*al[j])*par.git;
     sigma *= (1.0/ionicpar.Am);
   }
   //! diffusion tensor; extracellular for compatibility
   void De (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
			typename Traits::TensorType& sigma) const
   {
	  sigma = 0.0;
   } 
 
   
   //initialize the solution for parabolic and ode solution
   template<typename GFS, typename V>
   void InitialSolution(const GV& gv, GFS& vgfs, V& vinit)
   {
	  InitialSolution_MSR<GV,RF,ncomp> init(gv);  
	  Dune::PDELab::interpolate(init,vgfs,vinit);
	  return;
   }

    //initialize the solution for elliptic, parabolic and ode solution
    template<typename GFS, typename V>
    void InitialSolution(const GV& gv, GFS& vgfs, V& vinit, bool is_bid)
   {
          InitialSolution_MSR<GV,RF,ncomp+1> init(gv);
          Dune::PDELab::interpolate(init,vgfs,vinit);
          return;
   }

#if EXPLICIT_SCHEMES
 //initialize the solution for parabolic and ode solution
    template<typename VGFS, typename WGFS, typename V, typename W>
    void InitialSolution(const GV& gv, VGFS& vgfs, WGFS& wgfs, V& vinit, W& winit)
   {
	  InitialSolution_MSR_V<GV,RF,1> initv(gv);  
	  Dune::PDELab::interpolate(initv,vgfs,vinit);

          InitialSolution_MSR_W<GV,RF,1> initw(gv);  
	  Dune::PDELab::interpolate(initw,wgfs,winit);
      
	  return;
   }
#endif
   // this is onle for ODE system solver
   template<typename V>
   void initialsolution_model(V& y)
   {
	  initialsolutionMSRmodel(y);
   }
   
   //initialize the solution for parabolic and ode solution
   template<typename GFS, typename V>
   void InitialSolutionOptim(const GV& gv, GFS& vgfs, V& vinit)
   {
	  InitialSolution_MSROptim<GV,RF,ncomp> init(gv);  
	  Dune::PDELab::interpolate(init,vgfs,vinit);
	  return;
   }

#if EXPLICIT_SCHEMES 
   // initializing the stimulus for S1-S2 protocal
   template<typename GFS, typename V>
   void InitialStimulus(const GV& gv, GFS& vgfs, V& Itr1, V& Itr2) const
   {
	  Stimulus1_MSR<GV,RF,1> init1(gv,dim);  
	  Dune::PDELab::interpolate(init1,vgfs,Itr1);
	    
	  Stimulus2_MSR<GV,RF,1> init2(gv,dim);  
	  Dune::PDELab::interpolate(init2,vgfs,Itr2);
	  return;
   }
#else
 // initializing the stimulus for S1-S2 protocal
   template<typename GFS, typename V>
   void InitialStimulus(const GV& gv, GFS& vgfs, V& Itr1, V& Itr2) const
   {
	  Stimulus1_MSR<GV,RF,ncomp> init1(gv,dim);  
	  Dune::PDELab::interpolate(init1,vgfs,Itr1);
	    
	  Stimulus2_MSR<GV,RF,ncomp> init2(gv,dim);  
	  Dune::PDELab::interpolate(init2,vgfs,Itr2);
	  return;
   }
#endif

    RF getStimulusSupport(const typename Traits::DomainType& x) const
    {
        RF ret_val = 0.0;
        if(time>=S1_time && time<=S1_time + S1S2_delta)
            ret_val = Istim*S1Support(x,dim);
        else if(time+dt>=S2_time-1E-2 && time<=S2_time+2.0)
            ret_val = Istim*S2Support(x,dim);
        else
            ret_val = 0.0;
	 // if(ret_val > 0.0)
	 //   std::cout<<" ret_val "<<ret_val<<" "<<time<<" "<<Istim<<" "<<x[0]<<" "<<x[1]<<std::endl;  
        return ret_val;
    }

  //! set time for subsequent evaluation
  void setIstim (RF value) const
  {
  }  
   //! set time for subsequent evaluation
   void setTime (RF t) const
   {
	  time = t;
   }
 
   //! set time step for subsequent evaluation
   void setTimeStep (RF tau) const
   {
	  dt = tau;
   }
 
   std::string Name()
   {
	  std::string name = "==== Revised Mitchell-Shaeffer IonicModel ";
	  return name;
   }

   std::string component_name(int numb) const
   {
	  if(numb==0)
		 return "v";
	  else if(numb==1)
		 return "w";
	  else
	    return "u";
   }
private:
   typename Traits::TensorType sigma_i;

public:
    mutable RF time;
    mutable RF dt;
    RF maxtau, Istim;
    RF S1_time, S2_time, S1S2_delta;
    RF time_dep, maxtau_dep;
    double gil,git,gel,get;
   double lambda;
   RF G, vth, vp, n1, n2, n3, beta, xi;
};

}
#endif
