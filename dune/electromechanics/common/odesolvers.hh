/* ***************************************** */
/* File    : main function                   */
/*                                           */  
/* Date    : 19-05-2007                      */
/*                                           */
/* History : created by Chamakuri Nagaiah    */
/*                                           */
/* Remarks :                                 */
/*********************************************/


//ODEsolver
class ODESolver
{
public:
   typedef double Real;
   typedef double RT;
  
 //Constructor
    ODESolver(Real ATol_ = 1e-3) :  ATol(ATol_), RTol(1.0)
    {
        maxIterations=1000; 
        reduction=1E-6;
        error_old = 1e-5;
        tau_prev = 1e-5;
        facmax = 2.0;
        facmin = 0.1;
        fac = 0.95; 
        TTol = ConfigParser::get("RBSolver.tautol",1e-3);//1e-3;//5E-6;//5e-3;//1.0e-6;
        ATol = ConfigParser::get("RBSolver.atol",1e-3);
        RTol = ConfigParser::get("RBSolver.rtol",0.1);
        first = true;
    }


    // implementation from Lang-Teleaga ROS3PL paper
    bool PIController(int emb_order, Real maxtau, Real error_new, Real sol_norm, Real& tau_old, Real& tau_new) 
    {
        // is it first time step
        if(first)
            {
                tau_new = fac * tau_old * std::sqrt( TTol/error_new );
                tau_prev = tau_old; 
                error_old = error_new;
                first = false;
                tau_new = std::min(maxtau, std::min(facmax * tau_old,tau_new)); //restrict the time step
                return true;
            }
        
        // computing the error 
        Real error = std::sqrt(error_new*error_new/(ATol + RTol*sol_norm));

        // Now decide time step accept/reject
        bool accept = (error < TTol);
        // std::cout<<"                norms "<<error_new<<" "<<std::abs(std::isnan(error_new))<<" "<< (error_new!=error) <<" "<<isfinite(error_new)<<std::endl;


        if(!accept || !(error_new!=error))
            {
		if(Dune::MPIHelper::getCommunication().rank()==0)
                    std::cout<<" Rejected step : Error: "<< error_new <<" newError: " <<error <<" TTol: "<<TTol <<std::endl;
                                                                                                           
                tau_new = tau_old/2.0; // maximum can be reduced by factor of 2.0

		if(std::abs(tau_new) < 1e-12 || !isfinite(tau_new))
		  {
		    std::cout<<" Old time step is very small: exiting ..."<<tau_old<<std::endl;
		    std::cout<<" Error: "<< error <<" TTol: "<<TTol <<std::endl;
		    exit(1);
		  }
                return false;

            }
        
        // if time step accepts
        tau_new = tau_old/tau_prev * std::pow((TTol * error_old/(error * error)),1.0/(emb_order+1)) * tau_old;

        tau_new = std::min(maxtau, std::min(facmax * tau_old,tau_new));
        //std::cout<<"pi : "<<error_new<<" "<<error<<" "<<sol_norm<<" "<<tau_old<<" "<<tau_new<<std::endl;
	
        tau_prev = tau_old;
        error_old = error;

        if(std::abs(tau_new) < 1e-10)
            {
                std::cout<<" Old time step is very small: exiting ..."<<tau_old<<std::endl;
                std::cout<<" Error: "<< error <<" TTol: "<<TTol <<std::endl;
                exit(1);
            }
        
        return true;
    }  



    
public:

    //constants that are needed for the numerical calculation
    int maxIterations;
    double reduction;
    Real TTol, ATol, RTol;
    Real error_old, tau_prev ;
    Real facmax , facmin, fac, maxtau;
    bool first;
};

//#include "operatorsplittingsolvers.hh"
#include "rosenbrocksolver.hh"
