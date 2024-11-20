#ifndef IONICMODEL_MSR_BASE_HH
#define IONICMODEL_MSR_BASE_HH


namespace MitchellShaefferRegularizedModel
{

    // test model problem for 1 component
template<class GV, class Real, int ncomp=1>
class MSRModelParameters
{
   enum{Vm=0,Gw=1};
    
   typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,Real,ncomp> Traits;
public:
    
   MSRModelParameters()
   {
	  // tau_in   = 16.0;        tau_out    = 360.0;
	  // tau_open = 300.0;       tau_close  = 100.0;
	  tau_in = 0.3; // proportiional to APD increase
	  tau_out = 6.0;
	  tau_open = 120; 
	  tau_close = 150.0;
	  tau_close_RV = 120.0;   tau_close_m= 140.0;
	  tau_close_endo = 130.0; tau_close_epi = 90.0;
	  v_gate = -67.0;         v_rest     = -80.0;
	  v_peak = 20.0;

	  tau_in = ConfigParser::get("MSRModel.tau_in",0.3);
	  tau_out = ConfigParser::get("MSRModel.tau_out",6.0);
	  tau_open = ConfigParser::get("MSRModel.tau_open",120.0);
	  tau_close = ConfigParser::get("MSRModel.tau_close",150.0);
	  v_gate = ConfigParser::get("MSRModel.v_gate",-67.0);
	  v_rest = ConfigParser::get("MSRModel.v_rest",-80.0);
	  v_peak = ConfigParser::get("MSRModel.v_peak",20.0);
	  kappa = ConfigParser::get("MSRModel.kappa",1000.0);
	  Am = ConfigParser::get("MSRModel.Am",200.0);
	  Cm = ConfigParser::get("MSRModel.Cm",0.001);

	  g_Na = 1.0/tau_in; // [ms^-1]
	  g_K  = 1.0/tau_out; // [ms^-1]
	  //kappa = 100.0; // if kappa goes to infty then this is equivalent to original MS model
   }

   Real mv(typename Traits::RangeType& v) const
   {
	  if(v[Vm]<v_rest)
		 return 0;
	  else if(v[Vm]<v_peak)
		 return (v[Vm]-v_rest)/(v_peak-v_rest);
	  else
		 return 1.0;//v_peak;
	    
   }
 
   Real mv_v(typename Traits::RangeType& v) const
   {
	  if(v[Vm]<v_rest)
		 return 0;
	  else if(v[Vm]<v_peak)
		 return 1.0/(v_peak-v_rest);
	  else
		 return 0.0;
	    
   }
   Real mv_vv(typename Traits::RangeType& v) const
   {
	  return 0.0;
   }
   
   // function f
   Real fv(typename Traits::RangeType& v) const
   {
       Real v_phy = (v[Vm]-v_rest)/(v_peak-v_rest)-(v_gate-v_rest)/(v_peak-v_rest);
       return 1.0/2.0*(1.0+std::tanh(kappa*v_phy));
   }
    
    Real fv_v(typename Traits::RangeType& v) const
    {
	Real v_phy = (v[Vm]-v_rest)/(v_peak-v_rest)-(v_gate-v_rest)/(v_peak-v_rest);
	Real v_v = 1.0/(v_peak-v_rest);
	return 1.0/2.0*kappa*v_v*(1.0-std::tanh(kappa*v_phy)*std::tanh(kappa*v_phy));
    }
    
    Real fv_vv(typename Traits::RangeType& v) const
    {
	Real v_phy = (v[Vm]-v_rest)/(v_peak-v_rest)-(v_gate-v_rest)/(v_peak-v_rest);
	Real v_v = 1.0/(v_peak-v_rest);
	Real tanh_der = kappa*v_v*(1.0-std::tanh(kappa*v_phy)*std::tanh(kappa*v_phy));
	return -kappa*v_v*std::tanh(kappa*v_phy) * tanh_der;
    } 

   //function a
   inline Real av(typename Traits::RangeType& v) const
   {
	  Real f_con = fv(v);
	  return (1.0-f_con)/(tau_open+(tau_close-tau_open)*f_con);
   } 
   inline Real av_v(typename Traits::RangeType& v) const
   {
	  Real f_con = fv(v);
	  Real f_v = fv_v(v);
	  Real a = 1.0-f_con;
	  Real a_v = -f_v;
	  Real b = (tau_open+(tau_close-tau_open)*f_con);
	  Real b_v = (tau_close-tau_open)*f_v;
	  return (b*a_v-a*b_v)/(b*b);
   }  
   inline Real av_vv(typename Traits::RangeType& v) const
   {
	  Real f_con = fv(v);
	  Real f_v = fv_v(v);
	  Real f_vv = fv_vv(v);
	  Real a = 1.0-f_con;
	  Real a_v = -f_v;
	  Real a_vv = -f_vv;
	  Real b = (tau_open+(tau_close-tau_open)*f_con);
	  Real b_v = (tau_close-tau_open)*f_v;
	  Real b_vv = (tau_close-tau_open)*f_vv;
	  
	  Real part_a = b*b*(b_v*a_v + b*a_vv - a_v*b_v - a*b_vv);
	  Real part_b = 2.0*b*b_v*(b*a_v-a*b_v);
	  return (part_a-part_b)/std::pow(b,4.0);
   } 
  
   // function b
   inline Real bv(typename Traits::RangeType& v) const
   {
	  Real f_con = fv(v);
	  return f_con/(tau_open+(tau_close-tau_open)*f_con);
   } 
   inline Real bv_v(typename Traits::RangeType& v) const
   {
	  Real a = fv(v);
	  Real b = (tau_open+(tau_close-tau_open)*a);
	  Real a_v = fv_v(v);
	  Real b_v = (tau_close-tau_open)*a_v;
	  return (b*a_v-a*b_v)/(b*b);	
   }
   inline Real bv_vv(typename Traits::RangeType& v) const
   {
	  Real a    = fv(v);
	  Real a_v  = fv_v(v);
	  Real a_vv = fv_vv(v);
	  Real b    = (tau_open+(tau_close-tau_open)*a);
	  Real b_v  = (tau_close-tau_open)*a_v;
	  Real b_vv = (tau_close-tau_open)*a_vv;
	  Real part_a = b*b*(b_v*a_v + b*a_vv - a_v*b_v - a*b_vv);
	  Real part_b = 2.0*b*b_v*(b*a_v-a*b_v);
	  return ((part_a-part_b)/std::pow(b,4.0));	
   }
   
   // ionic part
   Real Iion(typename Traits::RangeType& v) const
   {
	  Real m_v = mv(v);
	  Real a = -g_Na*v[Gw]*m_v*m_v*(v_peak-v[Vm]);
	  Real b = g_K*(v[Vm]-v_rest);///(v_peak-v_rest);
	  return a+b;
   }
   
   // derivative of Iion
   Real Iion_v(typename Traits::RangeType& v) const
   {
	  Real m_v = mv(v);
	  Real mv_der = mv_v(v);
	  Real a = -g_Na*v[Gw]*(2.0*m_v*mv_der*(v_peak-v[Vm])-m_v*m_v);
	  Real b = g_K;
	  return a+b;
   }

   Real Iion_vv(typename Traits::RangeType& v) const
   {
	  Real m_v = mv(v);
	  Real mv_der = mv_v(v);
	  Real mv_der_der = mv_vv(v); // anyway it is zero
	  Real a = -g_Na*v[Gw]*(2.0*mv_der*(mv_der*(v_peak-v[Vm])-m_v)-2.0*m_v*mv_der);
	  Real b = 0.0;
	  return a+b;
   }  
   
   Real Iion_vw(typename Traits::RangeType& v) const
   {
	  Real m_v = mv(v);
	  Real mv_der = mv_v(v);
	  Real a = -g_Na*(2.0*m_v*mv_der*(v_peak-v[Vm])-m_v*m_v);
	  Real b = 0.0;
	  return a+b;
   }  
   
   
   // derivative of Iion
   Real Iion_w(typename Traits::RangeType& v) const 
   {
	  Real m_v = mv(v);
	  return (-g_Na*m_v*m_v*(v_peak-v[Vm] ));
   }
   
   Real Iion_wv(typename Traits::RangeType& v) const 
   {
	  Real m_v = mv(v);
	  Real mv_der = mv_v(v);
	  return (-g_Na*(2.0*m_v*mv_der*(v_peak-v[Vm])-m_v*m_v));
   }
   // derivative of Iion
   Real Iion_ww(typename Traits::RangeType& v) const 
   {
	  return 0.0;
   }
   // Gating variable RHS
   Real Gating(typename Traits::RangeType& v) const
   {
	  return av(v)*(1.0-v[Gw])-bv(v)*v[Gw];    
   }
#if EXPLICIT_SCHEMES
    // Gating variable LHS
    Real Gating_lhs(typename Traits::RangeType& v) const
   {
       return (av(v)+bv(v));    
   }
    // Gating variable RHS
   Real Gating_rhs(typename Traits::RangeType& v) const
   {
	  return av(v);    
   }  
#endif
    
   // derivative w.r.t v
   Real Gating_v(typename Traits::RangeType& v) const
   {
	  return av_v(v)*(1.0-v[Gw])-bv_v(v)*v[Gw];    
   }
 
   // derivative w.r.t v and again w.r.t. v
   Real Gating_vv(typename Traits::RangeType& v) const
   {
	  return av_vv(v)*(1.0-v[Gw])-bv_vv(v)*v[Gw];    
   }

   // derivative w.r.t v and w
   Real Gating_vw(typename Traits::RangeType& v) const
   {
	  return -av_v(v)-bv_v(v);    
   }
   
   // derivative w.r.t w
   Real Gating_w(typename Traits::RangeType& v) const
   {
	  return -av(v)-bv(v);
   }

   // derivative w.r.t w
   Real Gating_wv(typename Traits::RangeType& v) const
   {
	  return -av_v(v)-bv_v(v);
   }
   // derivative w.r.t w
   Real Gating_ww(typename Traits::RangeType& v) const
   {
	  return 0.0;
   }  
   Real tau_in, tau_out, tau_open, tau_close;
   Real tau_close_RV, tau_close_m, tau_close_endo, tau_close_epi;
   Real v_gate, v_rest, v_peak;
   Real g_Na, g_K, kappa;
   Real Am, Cm;
};


}
#endif
