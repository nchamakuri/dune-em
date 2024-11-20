/* ***************************************** */
/* File    : main function                   */
/*                                           */  
/* Date    : 19-05-2007                      */
/*                                           */
/* History : created by Chamakuri Nagaiah    */
/*                                           */
/* Remarks :                                 */
/*********************************************/


//Linearly implicit schemes with 2 stages
class LIRK2stages : public ODESolver
{

public:
   //constructor
   LIRK2stages ()
   {
	  gamma[0] = 0.0; 
	  gamma[1] = 0.0; 
	  gamma[2] = 0.0; 
         
	  alpha[0] = 0.0; 
	  alpha[1] = 0.0; 
     	
	  m1[0] = 0.0;   
	  m1[1] = 0.0;
         
	  m1b[0] = 0.0;  
	  m1b[1] = 0.0;  
    
	  a[0][0] = 0.0; a[0][1] = 0.0; 
	  a[1][0] = 0.0; a[1][1] = 0.0; 
     
	  c[0][0] = 0.0; c[0][1] = 0.0; 
	  c[1][0] = 0.0; c[1][1] = 0.0; 
     
	  // One needs to tweak these parameters in order to speed up the calculations
	  // WARNING: Becareful with tolerance for time step, if chosen wrong tolerance may lead
	  // to spurious results in calculations
	  tautol = 2.0e-3; 
	  facmax = 2.0;
	  facmin = 0.1;
	  fac = 0.95; 
   }

   //solver interface 
   template<class GV,class GFS, class GOS, class J, class M, class A,class ODE, class S, class PSP, class V, class PARAM>
   bool odesolve_primal(const GV& gv, const GFS& gfs,const GOS& gos, const J& j, const M& m, const A& ai, const ODE& coeff, S& solver, PSP& psp, V& usol, V& sol, V& Ie, V& f, RT t, RT *tau, PARAM& params, RT *err, int gfssize = 1, bool adaptive=false) 
   {
	  RT t_temp;
	  V rhs(sol);
	  V k1(sol);
	  V k2(sol);
	  V temp1(sol);
	  V temp2(Ie);
	  RT dt = (*tau);
	  int ncomp=1; 
   
	  //////////////////////////////////////////////// 
	  //////////First step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
 
	  // here alpha_1 = 0  
	  t_temp = t + coeff.alpha[0]*dt;

	  rhs = 0.0;
    
	  // Evaluation of f(.,.) 
	  gos.residual(sol,rhs);
 
	  //temp1= 0.0; m.umv(rhs,temp1);
	  //    Dune::PDELab::FilenameHelper fn("monoODE_Q1");
	  //     vtkout(gv,gfs,rhs,fn.getName());exit(1); 
	  // apply control
	  temp2 = 0.0;
	  m.mv(Ie,temp2);
	  //ai.usmv(-1.0,usol,temp2); 
	  temp2 -= usol;
	  ai.usmv(-1.0,sol,rhs);
    
	  rhs += temp2;
 
	  // rhs *= (dt*coeff.gamma[0]);
	  //std::cout<<" norms "<< rhs.two_norm()<<std::endl;
	  // initialization for the first vector
	  k1 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k1,rhs,t);
   
	  //solving J*k1 = rhs with Solver and storing the solution in k1
	  //for all stages the left hand side matrix is same
	  //solver.apply(k1,rhs,res1);
	  solver.apply(j,k1,rhs,reduction);
   
	  //temporary vector which is useful in next stage computation
	  temp1 = k1;

	  //////////////////////////////////////////////// 
	  //////////second step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[1] * dt);
  
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 *=  coeff.a[1][0];
	  temp1 += sol;
	  //std::cout<<" k1 norms "<< k1.two_norm()<<"  "<< sol.two_norm()<<"  "<< rhs.two_norm()<<std::endl;        
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0;
	  gos.residual(temp1,rhs);

	  ai.usmv(-1.0,temp1,rhs); 

	  // aplly control
	  rhs += temp2;  
	  //rhs *= (dt*coeff.gamma[0]);
   
   
	  //evaluation of the rest of the right hand side:
	  //calculate M*tau*Gammaa*c21/tau*k1
	  temp1 = k1;
	  temp1 *= (coeff.c[1][0]/dt);          //Gammaa*c21*k1
	  //temp1 *= (coeff.gamma[0]*coeff.c[1][0]);    
	  //needs mass matrix-vector-multiplication
	  //M.umv(temp1,rhs) means rhs = rhs + M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the second stage vector
	  k2 = 0.0; 
   
	  //std::cout<<" norms "<< rhs.two_norm()<<std::endl;exit(1);
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k2,rhs,t);
   
	  //solving J*k2 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //solver.apply(k2,rhs,res2);
   
	  solver.apply(j,k2,rhs,reduction);
  
	  ////////////////////////////////////////////
	  ////////////// New solution ////////////
	  ////////////////////////////////////////////   
       
	  //sol = (oldsol) + (m1*k1) + (m2*k2)
    
	  if(adaptive)
		 {   

			temp1 = sol;
			temp1.axpy(coeff.m1[0],k1);
			temp1.axpy(coeff.m1[1],k2); 
	
			//////the second new solution constructed for adaptive time step
			//rhs = (oldsol) + (m1b*k1) + (m2b*k2)
			rhs = sol;
			rhs.axpy(coeff.m1b[0],k1);
			rhs.axpy(coeff.m1b[1],k2);
	
			//find the error between two solutions
			rhs -= temp1;
	
			// using L2-norm
			double error = solver.norm(rhs)/gfssize;
	
	
			bool accept;
			if(error < tautol)
			   accept = true;
			else
			   accept = false;
	
			double taunew;
			// calculate the new time step
			taunew = fac * dt * std::sqrt( tautol/error );
	
			if (taunew > facmax * dt)
			   taunew = facmax * dt;
			if (taunew < facmin * dt)
			   taunew = facmin * dt;
			if(taunew > maxtau)
			   taunew = maxtau;
	
			if(accept)
			   {
				  sol = temp1;
				  *tau = taunew;
				  return accept;
			   }
			else 
			   {
				  *tau = taunew; *err = error;
				  if(gv.comm().rank()==0)
					 std::cout<<" Step is rejcted: time := "<<t<<" tauold = "<<dt <<" taunew := "<<taunew <<" error = "<<error<<std::endl;
				  return accept;
			   }
	
		 }
	  else
		 {
			sol.axpy(coeff.m1[0],k1);
			sol.axpy(coeff.m1[1],k2);
		 }
   
	  return true;
      
   }

   //  //solver interface 
   template<class GV,class GOS, class J, class M, class A, class ODE, class S, class PSP, class V, class PARAM>
   bool odesolve_dual(const GV& gv, const GOS& gos, const J& j, const M& m, const A& ai, const ODE& coeff, S& solver, PSP& psp, V& p, V& sol, V& b, RT t, RT *tau, PARAM& params, RT *err)   
   {
	  RT t_temp;
	  V rhs(sol);
	  V k1(sol);
	  V k2(sol);
	  V temp1(sol);
	  V temp2(sol);
	  RT dt = (*tau);
	  int ncomp=1; 
  
	  //////////////////////////////////////////////// 
	  //////////First step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////

	  // here alpha_1 = 0  
	  t_temp = t + coeff.alpha[0]*dt;

	  rhs = 0.0;
    
	  gos.residual(sol,rhs);
    
	  ai.umv(sol,rhs);
    
	  // apply rest of the vectors
	  //   if(dual_act == 0)
	  //       {
	  temp2 = 0.0;
	  m.usmv(-1.0,b,temp2);//m.usmv(-1.0,b,temp2);
	  //	ai.umv(p,temp2); // from diffusion of \phi_e 
	  temp2 += p;
	  //temp2 -= der; // add derivative contribution 
	  rhs += temp2;
	  //   }

	
	  //   Dune::PDELab::FilenameHelper fn("monodual_Q1");
	  //     vtkout_2comps(gv,gfs,temp1,rhs,fn.getName());exit(1);
    
	  // initialization for the first vector
	  k1 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k1,rhs,t);
   
	  //solving J*k1 = rhs with Solver and storing the solution in k1
	  //for all stages the left hand side matrix is same
	  //solver.apply(k1,rhs,res1);
	  solver.apply(j,k1,rhs,reduction);
   
	  //temporary vector which is useful in next stage computation
	  temp1 = k1;
    
	  //////////////////////////////////////////////// 
	  //////////second step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[1] * dt);
  
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 *=  coeff.a[1][0];
	  temp1 += sol;
            
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0;
	  gos.residual(temp1,rhs);
  
	  ai.umv(temp1,rhs); 
	  // aplly control
	  rhs += temp2;
  
	  //evaluation of the rest of the right hand side:
	  //calculate M*tau*Gammaa*c21/tau*k1
	  temp1 = k1;
	  temp1 *= (coeff.c[1][0]/dt);          //Gammaa*c21*k1
	  //temp1 *= (coeff.gamma[0]*coeff.c[1][0]);  
	  //needs mass matrix-vector-multiplication
	  //M.umv(temp1,rhs) means rhs = rhs + M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
   
    
	  //initialization for the second stage vector
	  k2 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k2,rhs,t);
   
	  //solving J*k2 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //solver.apply(k2,rhs,res2);  
    
	  solver.apply(j,k2,rhs,reduction);
  
	  ////////////////////////////////////////////
	  ////////////// New solution ////////////
	  ////////////////////////////////////////////   
       
	  //sol = (oldsol) + (m1*k1) + (m2*k2)
	  sol.axpy(coeff.m1[0],k1);
	  sol.axpy(coeff.m1[1],k2);
  
	  return true;
   }

  
public:
   //coefficients
   typedef Dune::FieldVector<double,2> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,2,2> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,3> LocalVector3;
   LocalVector3 gamma;  
   // tolerance for the new time step
   double tautol; 
   // factors the pick up the new time step
   double facmax , facmin, fac, maxtau;   
};



//RoS2 method 
class ROS2 : public LIRK2stages
{
 
public:
   //constructor
   ROS2() : LIRK2stages()
   { 
        
	  gamma[0]=  1.707106781186547;
	  gamma[1]=  1.707106781186547; 
	  gamma[2]= -1.707106781186547;
	
	  alpha[0]= 0.0;  
	  alpha[1]= 1.0; 
	
	  m1[0]= 0.8786796564403575;   
	  m1[1]= 0.2928932188134525;    
	
	  m1b[0]= 0.5857864376269050;
	  m1b[1]= 0.0;
       
	  a[0][0] = 0.0;                 a[0][1] = 0.0;              
	  a[1][0] = 0.5857864376269050;  a[1][1] = 0.0;               
	
	  c[0][0] = 0.5857864376269050;  c[0][1] = 0.0;                
	  c[1][0] = 1.171572875253810;   c[1][1] = 0.5857864376269050; 
   }

public:
   //coefficients
   typedef Dune::FieldVector<double,2> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,2,2> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,3> LocalVector3;
   LocalVector3 gamma;
};





//Linearly implicit schemes with 3 internal stages
// at every stage one has to solve linear system
template<class GV>
class LIRK3stages : public ODESolver
{
   enum{stages=3};
public:
   //constructor   
   LIRK3stages(const GV& gv_, double maxtau_, unsigned gfssize_):
	  gv(gv_),maxtau(maxtau_),gfssize(gfssize_)
   {
	  gamma[0] = 0.0; 
	  gamma[1] = 0.0; 
	  gamma[2] = 0.0; 
	  gamma[3] = 0.0;
    
	  alpha[0] = 0.0; 
	  alpha[1] = 0.0; 
	  alpha[2] = 0.0;
	
	  m1[0] = 0.0;   
	  m1[1] = 0.0;
	  m1[2] = 0.0;
    
	  m1b[0] = 0.0; 
	  m1b[1] = 0.0;  
	  m1b[2] = 0.0; 
    
	  a[0][0] = 0.0; a[0][1] = 0.0; a[0][2] = 0.0; 
	  a[1][0] = 0.0; a[1][1] = 0.0; a[1][2] = 0.0;
	  a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 0.0;
    
	  c[0][0] = 0.0; c[0][1] = 0.0; c[0][2] = 0.0; 
	  c[1][0] = 0.0; c[1][1] = 0.0; c[1][2] = 0.0;
	  c[2][0] = 0.0; c[2][1] = 0.0; c[2][2] = 0.0; 

	  // I hope this should work
      first_iteration = true;
   }
   
   // run the ode solver; coefficients for ode solver is available via constructor
   template<class GOS, class MatrixType, class V, class Solver, class COEFF>
   bool solve(GOS& gos, MatrixType& j, const MatrixType& ai, const MatrixType& m, V& sol, V& rhs_const, Solver& solver, const COEFF& coeff, Real t, Real& tau, bool applydualpart=false, bool adaptive=false)
   {
	  Real t_temp;
	  V rhs(sol);
	  V k1(sol);
	  V k2(sol);
	  V k3(sol);
	  V temp1(sol);
	  Real dt = tau;
	  
	  using Dune::PDELab::Backend::native;
	  
	  static V last_k1(sol);
	  static V last_k2(sol); 
	  static V last_k3(sol);

	  if(first_iteration)
		 {
			last_k1 = 0.0; last_k2 = 0.0; last_k3 = 0.0;
			first_iteration = false;
		 }
	  
	  //////////////////////////////////////////////// 
	  //////////First step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////

	  // here alpha_1 = 0  
	  t_temp = t + coeff.alpha[0]*dt;

	  //f=calculate_f(tn+alpha1*tau,c)
	  rhs = 0.0;
    
	  // Evaluation of f(.,.) 
	  gos.residual(sol,rhs);
	  
	  // double p_norm = rhs.two_norm(); p_norm = gv.comm().sum(p_norm);
	  // double q_norm = sol.two_norm(); q_norm = gv.comm().sum(q_norm);
	  // //if(applydualpart) 
	  // 	 if(gv.comm().rank()==0) std::cout<<"   rhs norm " <<p_norm<<" "<<q_norm<<std::endl; 
	  //temp1= 0.0; m.umv(rhs,temp1);
	  //    Dune::PDELab::FilenameHelper fn("monoODE_Q1");
	  //     vtkout(gv,gfs,rhs,fn.getName());exit(1); 
	  // apply control
	  if(applydualpart) // dual solve
		 {
		   native(ai).umv(native(sol),native(rhs));
		 }
	  else // Primal solve
		 {
			//ai.usmv(-1.0,usol,temp2);
		   native(ai).usmv(-1.0,native(sol),native(rhs)); 	   
		 }
	  
	  // p_norm = rhs.two_norm(); p_norm = gv.comm().sum(p_norm);
	  // q_norm = j.frobenius_norm(); q_norm = gv.comm().sum(q_norm);
	  // //if(applydualpart) 
	  // 	 if(gv.comm().rank()==0) std::cout<<"   rhs norm1 " <<p_norm<<" "<<q_norm<<std::endl; 
		 
	  // finally add the const contribution of RHS 
	  rhs += rhs_const;
	
	  // initialization for the first vector
	  //if(adaptive)
		 k1 = last_k1;
	  //else
	  //	 k1 = 0.0;
// double	   p_norm = rhs.two_norm(); p_norm = gv.comm().sum(p_norm);
// double	   q_norm = k1.two_norm(); q_norm = gv.comm().sum(q_norm);
// 	  // //if(applydualpart) 
// 	   if(gv.comm().rank()==0) std::cout<<"   rhs norm1 " <<p_norm<<" "<<q_norm<<std::endl; 
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k1,rhs,t);
   
	  //solving J*k1 = rhs with Solver and storing the solution in k1
	  //for all stages the left hand side matrix is same
	  //solver.apply(k1,rhs,res1);
	  //Dune::PDELab::set_constrained_dofs(cc,0.0,rhs);
	  solver.apply(j,k1,rhs,reduction);
	  
	  //temporary vector which is useful in next stage computation
	  temp1 = k1;
	  //if(adaptive)
		 last_k1 = k1;
	  
	  //////////////////////////////////////////////// 
	  //////////second step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[1] * dt);
  
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 *=  coeff.a[1][0];
	  temp1 += sol;
            
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);
	  // if(applydualpart) 
	  // if(gv.comm().rank()==0) std::cout<<" rhs norm2 " <<rhs.two_norm()<<" "<<k1.two_norm()<<" "<<temp1.two_norm()<<std::endl;
	  if(applydualpart)
		 {
		   native(ai).umv(native(temp1),native(rhs)); 
		 }
	  else
		 {
		   native(ai).usmv(-1.0,native(temp1),native(rhs)); 
		 }
	  //ai.usmv(-1.0,temp1,rhs); 
	  // aplly control
	  rhs += rhs_const;
      //if(applydualpart) 
	  //	  if(gv.comm().rank()==0) std::cout<<" rhs norm2 " <<rhs.two_norm()<<" "<<rhs_const.two_norm()<<" "<<temp1.two_norm()<<" "<<ai.base().frobenius_norm()<<std::endl;
	  // temporary vector to store the rhs for the third stage
	  // at present k3 is temp vector
	  k3 = rhs;
    
	  //evaluation of the rest of the right hand side:
	  //calculate M*tau*Gammaa*c21/tau*k1
	  temp1 = k1;
	  temp1 *= (coeff.c[1][0]/dt);          //Gammaa*c21*k1
      
	  //needs mass matrix-vector-multiplication
	  //M.umv(temp1,rhs) means rhs = rhs + M*temp1
	  native(m).usmv(-1.0,native(temp1),native(rhs));
	  // if(applydualpart) if(gv.comm().rank()==0) std::cout<<" rhs norm2:final " <<rhs.two_norm()<<" "<<dt<<" "<<coeff.c[1][0]<<" "<<m.frobenius_norm()<<" "<<temp1.two_norm()<<" "<<k1.two_norm()<<std::endl;
	  //initialization for the second stage vector
	  //if(adaptive)
		 k2 = last_k2; 
	  // else
	  // 	 k2 = 0.0;
	  //if(applydualpart) 
		 //if(gv.comm().rank()==0) std::cout<<"   k2: rhs norm " <<rhs.two_norm()<<" "<<k2.two_norm()<<std::endl;
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k2,rhs,t);
   
	  //solving J*k2 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //solver.apply(k2,rhs,res2);
	  //Dune::PDELab::set_constrained_dofs(cc,0.0,rhs);
	  solver.apply(j,k2,rhs,reduction);
	  //if(adaptive) 
		 last_k2 = k2;
	    
	  // 	 p_norm = k3.two_norm(); p_norm = gv.comm().sum(p_norm);
	  // 	 q_norm = k2.two_norm(); q_norm = gv.comm().sum(q_norm);
	  // // //if(applydualpart) 
	  // 	 if(gv.comm().rank()==0) std::cout<<"   2rhs norm " <<p_norm<<" "<<q_norm<<std::endl; 
	  
	  //////////////////////////////////////////////// 
	  //////////third step of LIRK 3-stages///////////
	  ////////////////////////////////////////////////
    
	  // in third stage we don't need to calculate rhs b'cause a[3][2] = 0 
	  // and alpha[2] = alpha[3], k3 is released
	  rhs = k3;
    
	  //evaluation of the rest of the right hand side
	  // M * tau*Gamma *(c31*k1+c32*k2)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[2][0]/dt);               //Gammaa*c31*k1
	  temp1.axpy((coeff.c[2][1]/dt),k2);
    
	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
	  native(m).usmv(-1.0,native(temp1),native(rhs));
    
	  //initialization for the third stage vector
	  //if(adaptive) 
		 k3 = last_k3;
	  // else
	  // 	 k3 = 0.0;
	 
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k3,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //Dune::PDELab::set_constrained_dofs(cc,0.0,rhs);
	  solver.apply(j,k3,rhs,reduction);
	  //if(adaptive) 
		 last_k3 = k3;

	  // Real norm1, norm2, norm3, norm4;
	  // norm1 = solver.norm(k1);
	  // norm2 = solver.norm(k2);
	  // norm3 = solver.norm(k3);
	  // norm4 = solver.norm(sol); 
	  // if(gv.comm().rank()==0) 
	  // 	 std::cout<<"   on: norm " <<norm1<<" "<<norm2<<" "<<norm3<<" "<<norm4<<std::endl;
	  

	  ////////////////////////////////////////////
	  ////////////// New solution ////////////
	  ////////////////////////////////////////////   
       
	  //sol = (oldsol) + (m1*k1) + (m2*k2)

	  if(adaptive)
		 {   

			temp1 = sol;
			temp1.axpy(coeff.m1[0],k1);
			temp1.axpy(coeff.m1[1],k2); 
			temp1.axpy(coeff.m1[2],k3); 
	
			//////the second new solution constructed for adaptive time step
			//rhs = (oldsol) + (m1b*k1) + (m2b*k2)
			rhs = sol;
			rhs.axpy(coeff.m1b[0],k1);
			rhs.axpy(coeff.m1b[1],k2);
			rhs.axpy(coeff.m1b[2],k3);

			//find the error between two solutions
			rhs -= temp1;
	
			// using L2-norm; to avoid by zero add small number
			Real error = solver.norm(rhs)/gfssize + 1E-15;
            Real sol_norm = solver.norm (sol);

            // calculate the new timestep according to controller
            Real taunew;
            bool accept = PIController(1,maxtau,error,sol_norm,dt,taunew);

            
			if(accept)
			   {
				  sol = temp1;
				  tau = taunew;
				  return accept;
			   }
			else 
			   {
				  tau = taunew;
				  if(gv.comm().rank()==0)
                      std::cout<<" Step is rejcted: time := "<<t<<" tauold = "<<dt <<" taunew := "<<taunew <<" error = "<<error<<std::endl;
				  return accept;
			   }
	
		 }
	  else
		 {
			sol.axpy(coeff.m1[0],k1);
			sol.axpy(coeff.m1[1],k2);
			sol.axpy(coeff.m1[2],k3);
			// double p_norm = k1.two_norm(); p_norm = gv.comm().sum(p_norm);
			// double q_norm = k2.two_norm(); q_norm = gv.comm().sum(q_norm);
			// double r_norm = k3.two_norm(); r_norm = gv.comm().sum(r_norm);
			// double s_norm = sol.two_norm(); s_norm = gv.comm().sum(s_norm);
			// if(applydualpart) 
			//  if(gv.comm().rank()==0) 
			//     std::cout<<"   sol norm " <<p_norm<<" "<<q_norm<<" "<<r_norm<<" "<<s_norm<<std::endl; 
		 }
   
	  return true;
      
   }
  
    
private:
   // import stuff from PDE solver
   const GV& gv;
   // Solver& solver;
   unsigned gfssize;
public:  
   //coefficients
   typedef Dune::FieldVector<Real,stages> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<Real,stages,stages> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<Real,stages+1> LocalVectorGamma;
   LocalVectorGamma gamma;  

    // factors the pick up the new time step
   Real maxtau;
   bool first_iteration;
};




//RoS3P
template<class GV>
class ROS3p : public LIRK3stages<GV>
{
 
public:
   //constructor
   ROS3p(const GV& gv_,double maxtau_,unsigned gfssize_):LIRK3stages<GV>(gv_,maxtau_,gfssize_)
   { 
        
	  gamma[0]=  7.886751345948129e-1;
	  gamma[1]=  7.886751345948129e-1; 
	  gamma[2]= -2.11324865405187e-1;
	  gamma[3]= -1.077350269189626;
	
	  alpha[0]= 0.0;  
	  alpha[1]= 1.0; 
	  alpha[2]= 1.0;
	
	  m1[0]= 2.0;   
	  m1[1]= 5.773502691896258e-1;    
	  m1[2]= 4.226497308103742e-1;
	
	  m1b[0]= 2.113248654051871; /*2.22649730810374243;*/
	  m1b[1]= 1.000000000000000;   /*1.42264973081037426;*/
	  m1b[2]= 4.226497308103742e-1; /*0.42264973081037424;*/
	
	  a[0][0]=0.0;                a[0][1]=0.0;               a[0][2]=0.0; 
	  a[1][0]=1.267949192431123;  a[1][1]=0.0;               a[1][2]=0.0;
	  a[2][0]=1.267949192431123;  a[2][1]=0.0;               a[2][2]=0.0;
	
	  c[0][0]=1.267949192431123;  c[0][1]=0.0;               c[0][2]=0.0; 
	  c[1][0]=1.607695154586736;  c[1][1]=1.267949192431123; c[1][2]=0.0;
	  c[2][0]=3.464101615137755;  c[2][1]=1.732050807568877; c[2][2]=1.267949192431123;
   }

public:
   //coefficients
   typedef Dune::FieldVector<double,3> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,3,3> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,4> LocalVector4;
   LocalVector4 gamma;
};

//RoS3PW
template<class GV>
class ROS3PW: public LIRK3stages<GV>
{
   
public:
   //constructor
   ROS3PW(const GV& gv_, double maxtau_,unsigned gfssize_):LIRK3stages<GV>(gv_,maxtau_,gfssize_)
   { 
        
	  gamma[0] = 7.886751345948129e-1;
	  gamma[1]=  7.886751345948129e-1;
	  gamma[2]= -2.11324865405187e-1;    
	  gamma[3]= -1.077350269189626;
	
	  alpha[0]= 0.0;  
	  alpha[1]= 1.0; 
	  alpha[2]= 1.0;
	
	  m1[0]= 2.0;   
	  m1[1]= -0.16987298107781; 
	  m1[2]= 1.07179676972449;
	
	  m1b[0]= 1.06624327025936;
	  m1b[1]= 0.19059892324150;
	  m1b[2]= 1.07179676972449;
	
	  a[0][0]=0.0;               a[0][1]=0.0;               a[0][2]=0.0; 
	  a[1][0]=2.0;               a[1][1]=0.0;               a[1][2]=0.0;
	  a[2][0]=0.63397459621556;  a[2][1]=0.0;               a[2][2]=0.0;
	
	  c[0][0]=1.267949192431123; c[0][1]=0.0;               c[0][2]=0.0; 
	  c[1][0]=2.535898384862250; c[1][1]=1.267949192431123; c[1][2]=0.0;
	  c[2][0]=0.52932852445504;  c[2][1]=-0.27451905283833; c[2][2]=1.267949192431123;
   }

public:
   //coefficients
   typedef Dune::FieldVector<double,3> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,3,3> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,4> LocalVector4;
   LocalVector4 gamma;
 
};

//WMETHOD
template<class GV>
class WMETHOD : public LIRK3stages<GV>
{
public:
   //constructor
   WMETHOD(const GV& gv_, double maxtau_,unsigned gfssize_):LIRK3stages<GV>(gv_,maxtau_,gfssize_)
   { 
	  gamma[0] = 0.2928932188134524; 
	  gamma[1] = 0.0; 
	  gamma[2] = 0.0; 
	  gamma[3] = 0.0;
	
	  alpha[0] = 0.0; 
	  alpha[1] = 1.0; 
	  alpha[2] = 1.0;
	
	  m1[0] = 1.0;   
	  m1[1] = -0.20710678118654757;
	  m1[2] = 0.5;
	
	  m1b[0] =  0.8292893218813453; 
	  m1b[1] = -0.3278174593052024;  
	  m1b[2] =  0.6207106781186548;
	
	  a[0][0] = 0.0; a[0][1] = 0.0; a[0][2] = 0.0; 
	  a[1][0] = 1.0; a[1][1] = 0.0; a[1][2] = 0.0;
	  a[2][0] = 1.0; a[2][1] = 0.0; a[2][2] = 0.0;
	
	  c[0][0]=0.0;               c[0][1]=0.0;                  c[0][2]=0.0; 
	  c[1][0]=3.414213562373095; c[1][1]=0.0;                  c[1][2]=0.0;
	  c[2][0]=1.0;               c[2][1]=-0.41421356237309515; c[2][2]=0.0;
   }

public:
   //coefficients
   typedef Dune::FieldVector<double,3> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,3,3> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,4> LocalVector4;
   LocalVector4 gamma;
};

//ROWDA3
template<class GV>
class ROWDA3 : public LIRK3stages<GV>
{
public:
   //constructor
   ROWDA3(const GV& gv_, double maxtau_, unsigned gfssize_):LIRK3stages<GV>(gv_,maxtau_,gfssize_)
   { 
	  gamma[0] = 0.4358665215084590; 
	  gamma[1] = 0.4358665215084590; 
	  gamma[2] = 0.6044552840655588; 
	  gamma[3] = 6.3797887993448800;
	
	  alpha[0] = 0.0; 
	  alpha[1] = 0.7; 
	  alpha[2] = 0.7;
	
	  m1[0] =  2.236727045296589;   
	  m1[1] =  2.250067730969645;
	  m1[2] = -0.2092514044390320;
	
	  m1b[0] =  2.059356167645941; 
	  m1b[1] =  0.1694014319346527;  
	  m1b[2] =  0.0;
	
	  a[0][0] = 0.0;               a[0][1] = 0.0; a[0][2] = 0.0; 
	  a[1][0] = 1.605996252195329; a[1][1] = 0.0; a[1][2] = 0.0;
	  a[2][0] = 1.605996252195329; a[2][1] = 0.0; a[2][2] = 0.0;
	
	  c[0][0] =   2.294280360279042;  c[0][1] =  0.0;               c[0][2]=0.0; 
	  c[1][0] = - 0.8874044410657823; c[1][1] =  2.294280360279042; c[1][2]=0.0;
	  c[2][0] = -23.98747971635035;   c[2][1] = -5.263722371562130; c[2][2]=2.294280360279042;
   }

public:
   //coefficients
   typedef Dune::FieldVector<double,3> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,3,3> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,4> LocalVector4;
   LocalVector4 gamma;
};




//Linearly implicit schemes with 4 stages
template<class GV>
class LIRK4stages : public ODESolver
{
    enum{stages=3};
public:
   //constructor   
   LIRK4stages(const GV& gv_, double maxtau_, unsigned gfssize_): gv(gv_),maxtau(maxtau_),gfssize(gfssize_)
   {
	  gamma[0] = 0.0; 
	  gamma[1] = 0.0; 
	  gamma[2] = 0.0; 
	  gamma[3] = 0.0;
	  gamma[4] = 0.0;
	
	  alpha[0] = 0.0; 
	  alpha[1] = 0.0; 
	  alpha[2] = 0.0;
	  alpha[3] = 0.0;
	
	  m1[0] = 0.0;   
	  m1[1] = 0.0;
	  m1[2] = 0.0;
	  m1[3] = 0.0;
	
	  m1b[0] =  0.0; 
	  m1b[1] =  0.0;  
	  m1b[2] =  0.0;
	  m1b[3] =  0.0;
	
	  a[0][0] = 0.0; a[0][1] = 0.0; a[0][2] = 0.0; a[0][3] = 0.0;  
	  a[1][0] = 0.0; a[1][1] = 0.0; a[1][2] = 0.0; a[1][3] = 0.0;
	  a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 0.0; a[2][3] = 0.0;
	  a[3][0] = 0.0; a[3][1] = 0.0; a[3][2] = 0.0; a[3][3] = 0.0;
	  
	  c[0][0] = 0.0; c[0][1] = 0.0; c[0][2] = 0.0; c[0][3] = 0.0;    
	  c[1][0] = 0.0; c[1][1] = 0.0; c[1][2] = 0.0; c[1][3] = 0.0;
	  c[2][0] = 0.0; c[2][1] = 0.0; c[2][2] = 0.0; c[2][3] = 0.0;
	  c[3][0] = 0.0; c[3][1] = 0.0; c[3][2] = 0.0; c[3][3] = 0.0;

      first_iteration = true; 
   }

    //solver interface 
   // run the ode solver; coefficients for ode solver is available via constructor
    template<class GOS, class MatrixType, class V, class Solver, class COEFF>
   bool solve(GOS& gos, MatrixType& j, const MatrixType& ai, const MatrixType& m, V& sol, V& rhs_const, Solver& solver, const COEFF& coeff, Real t, Real& tau, bool applydualpart=false, bool adaptive=false) 
   {
	  RT t_temp;
	  V rhs(sol);
	  V k1(sol);
	  V k2(sol);
	  V k3(sol);
	  V k4(sol);
	  V temp1(sol);
	  V temp2(sol); 
	  RT dt = tau;

      using Dune::PDELab::Backend::native;
	  
	  static V last_k1(sol);
	  static V last_k2(sol); 
	  static V last_k3(sol);
      static V last_k4(sol);

      if(first_iteration)
          {
              last_k1 = 0.0; last_k2 = 0.0; last_k3 = 0.0;
              first_iteration = false;
          }

      
      //////////////////////////////////////////////// 
	  //////////First step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////

	  // here alpha_1 = 0  
	  t_temp = t + coeff.alpha[0]*dt;

	  //f=calculate_f(tn+alpha1*tau,c)
	  rhs = 0.0;
	  // Evaluation of f(.,.) 
	  gos.residual(sol,rhs); 
    
	  //temp1= 0.0; m.umv(rhs,temp1);
	  //    Dune::PDELab::FilenameHelper fn("monoODE_Q1");
	  //     vtkout(gv,gfs,rhs,fn.getName());exit(1);

      if(applydualpart) // dual solve
		 {
		   native(ai).umv(native(sol),native(rhs));
		 }
	  else // Primal solve
		 {
			//ai.usmv(-1.0,usol,temp2);
		   native(ai).usmv(-1.0,native(sol),native(rhs)); 	   
		 }
	 
    
      // finally add the const contribution of RHS 
	  rhs += rhs_const;
	
	  // initialization for the first vector
      k1 = last_k1;
      
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k1,rhs,t);
   
	  //solving J*k1 = rhs with Solver and storing the solution in k1
	  //for all stages the left hand side matrix is same
	  //solver.apply(k1,rhs,res1);
      
	  solver.apply(j,k1,rhs,reduction);
   
	  //temporary vector which is useful in next stage computation
	  temp1 = k1;
      last_k1 = k1;
      
	  //////////////////////////////////////////////// 
	  //////////second step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[1] * dt);
  
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 *=  coeff.a[1][0];
	  temp1 += sol;
            
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);

      if(applydualpart)
		 {
		   native(ai).umv(native(temp1),native(rhs)); 
		 }
	  else
		 {
		   native(ai).usmv(-1.0,native(temp1),native(rhs)); 
		 }

      // aplly control
	  rhs += rhs_const;
      
	  //evaluation of the rest of the right hand side:
	  //calculate M*tau*Gammaa*c21/tau*k1
	  temp1 = k1;
	  temp1 *= (coeff.c[1][0]/dt);          //Gammaa*c21*k1
      
	  //needs mass matrix-vector-multiplication
	  //M.umv(temp1,rhs) means rhs = rhs + M*temp1
      native(m).usmv(-1.0,native(temp1),native(rhs));
	  //initialization for the second stage vector
	  //k2 = 0.0;
      k2 = last_k2;
      
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k2,rhs,t);
   
	  //solving J*k2 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //solver.apply(k2,rhs,res2);
 
	  solver.apply(j,k2,rhs,reduction);
      last_k2 = k2;
      
	  //////////////////////////////////////////////// 
	  //////////third step of LIRK 3-stages///////////
	  ////////////////////////////////////////////////
    
	  //first argument of the function f
	  t_temp = t + (coeff.alpha[2] * dt);
    
	  //second argument of the function f
	  //temp1 = sol + a31 * k1 + a32 * k2 
	  temp1 = k1;
	  temp1 *= (coeff.a[2][0]);               
	  temp1.axpy(coeff.a[2][1],k2);
	  temp1 += sol;
    

	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);

      if(applydualpart)
		 {
		   native(ai).umv(native(temp1),native(rhs)); 
		 }
	  else
		 {
		   native(ai).usmv(-1.0,native(temp1),native(rhs)); 
		 }

      // aplly control
	  rhs += rhs_const;
     
	  // temporary vector to store the rhs for the third stage
	  // at present k3 is temp vector
	  k4 = rhs;
    
	  //evaluation of the rest of the right hand side
	  // M * tau*Gamma *(c31*k1+c32*k2)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[2][0]/dt);               //Gammaa*c31*k1
	  temp1.axpy((coeff.c[2][1]/dt),k2);
    
	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
      native(m).usmv(-1.0,native(temp1),native(rhs));
	  //initialization for the third stage vector
	  k3 = last_k2;
    
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k3,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  solver.apply(j,k3,rhs,reduction);
      last_k3 = k3;
      
	  //////////////////////////////////////////////// 
	  //////////fourth step of LIRK 4-stages///////////
	  //////////////////////////////////////////////// 
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[3] * dt); 

     
	  // in fourth stage we don't need to calculate rhs b'cause a43 = 0 
	  // and alpha[2] = alpha[3], k3 is released
	  rhs = k4;
    
	  //evaluation of the rest of the right hand side
	  // m * (c31*k1+c32*k2+c33*k3)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[3][0]/dt);               //Gammaa*c41*k1
	  temp1.axpy(coeff.c[3][1]/dt,k2);
	  temp1.axpy(coeff.c[3][2]/dt,k3);

	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
      native(m).usmv(-1.0,native(temp1),native(rhs));
	  //initialization for the third stage vector
	  k4 = last_k4;
    
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k4,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  solver.apply(j,k4,rhs,reduction);
      last_k4 = k4;
	  ////////////////////////////////////////////
	  ////////////// New solution ////////////
	  ////////////////////////////////////////////   
       
	  //sol = (oldsol) + (m1*k1) + (m2*k2)

	  if(adaptive)
		 {   

			temp1 = sol;
			temp1.axpy(coeff.m1[0],k1);
			temp1.axpy(coeff.m1[1],k2); 
			temp1.axpy(coeff.m1[2],k3); 
			temp1.axpy(coeff.m1[3],k4); 
	
			//////the second new solution constructed for adaptive time step
			//rhs = (oldsol) + (m1b*k1) + (m2b*k2)
			rhs = sol;
			rhs.axpy(coeff.m1b[0],k1);
			rhs.axpy(coeff.m1b[1],k2);
			rhs.axpy(coeff.m1b[2],k3);
			rhs.axpy(coeff.m1b[3],k4);
	
			//find the error between two solutions
			rhs -= temp1;
			
			// using L2-norm
			double error = solver.norm(rhs)/gfssize + 1E-15;
            Real sol_norm = solver.norm(sol);

            // calculate tau_new accoding controller
            Real taunew;
            bool accept = PIController(2,maxtau,error,sol_norm,dt,taunew);

            
			if(accept)
			   {
				  sol = temp1;
				  tau = taunew;
				  return accept;
			   }
			else 
			   {
				  tau = taunew; 
				  if(gv.comm().rank()==0)
                      std::cout<<" Step is rejcted: time := "<<t<<" tauold = "<<dt <<" taunew := "<<taunew <<" error = "<<error<<std::endl;
				  // for(int i=0;i<18;i++)
				  //   {
				  //     using Dune::PDELab::Backend::native; 
				  //     temp1 = 0.0;
				  //     for(int j=0;j<sol.N();j++)
				  // 	native(temp1)[i][0] = native(rhs)[j][i];
				  //     double norm_m = solver.norm(temp1);
				  //     if(norm_m>1e-1)
				  // 	if(gv.comm().rank()==0)
				  // 	  std::cout<<"  "<<i<<" "<<norm_m<<std::endl;
				  //   }
				  return accept;
			   }
	
	
		 }
	  else
		 {
			sol.axpy(coeff.m1[0],k1);
			sol.axpy(coeff.m1[1],k2);
			sol.axpy(coeff.m1[2],k3);
			sol.axpy(coeff.m1[3],k4);
		 }
   
	  return true;
      
   }
  


    
   //solver interface 
   template<class GFS, class GOS, class J, class M, class A, class ODE, class S, class PSP, class V, class PARAM>
   bool odesolve_primal(const GFS& gfs,const GOS& gos, const J& j, const M& m, const A& ai, const ODE& coeff, S& solver, PSP& psp, V& usol, V& sol, V& Ie, V& f, RT t, RT *tau, PARAM& params, RT *err, int gfssize = 1, bool adaptive=false) 
   {
	  RT t_temp;
	  V rhs(sol);
	  V k1(sol);
	  V k2(sol);
	  V k3(sol);
	  V k4(sol);
	  V temp1(sol);
	  V temp2(Ie); 
	  RT dt = (*tau);
	  int ncomp=1; 

	  //////////////////////////////////////////////// 
	  //////////First step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////

	  // here alpha_1 = 0  
	  t_temp = t + coeff.alpha[0]*dt;

	  //f=calculate_f(tn+alpha1*tau,c)
	  rhs = 0.0;
	  // Evaluation of f(.,.) 
	  gos.residual(sol,rhs); 
    
	  //temp1= 0.0; m.umv(rhs,temp1);
	  //    Dune::PDELab::FilenameHelper fn("monoODE_Q1");
	  //     vtkout(gv,gfs,rhs,fn.getName());exit(1); 
	  // apply control
	  temp2 = 0.0;
	  m.mv(Ie,temp2); 
	  //   ai.usmv(-1.0,usol,temp2); 
	  temp2 -= usol;
	  ai.usmv(-1.0,sol,rhs);
	  rhs += temp2;
    
	  // initialization for the first vector
	  k1 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k1,rhs,t);
   
	  //solving J*k1 = rhs with Solver and storing the solution in k1
	  //for all stages the left hand side matrix is same
	  //solver.apply(k1,rhs,res1);
  
	  solver.apply(j,k1,rhs,reduction);
   
	  //temporary vector which is useful in next stage computation
	  temp1 = k1;
    
	  //////////////////////////////////////////////// 
	  //////////second step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[1] * dt);
  
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 *=  coeff.a[1][0];
	  temp1 += sol;
            
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);
	  ai.usmv(-1.0,temp1,rhs); 
	  // aplly control
	  rhs += temp2;
     
	  //evaluation of the rest of the right hand side:
	  //calculate M*tau*Gammaa*c21/tau*k1
	  temp1 = k1;
	  temp1 *= (coeff.c[1][0]/dt);          //Gammaa*c21*k1
      
	  //needs mass matrix-vector-multiplication
	  //M.umv(temp1,rhs) means rhs = rhs + M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the second stage vector
	  k2 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k2,rhs,t);
   
	  //solving J*k2 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //solver.apply(k2,rhs,res2);
 
	  solver.apply(j,k2,rhs,reduction);

	  //////////////////////////////////////////////// 
	  //////////third step of LIRK 3-stages///////////
	  ////////////////////////////////////////////////
    
	  //first argument of the function f
	  t_temp = t + (coeff.alpha[2] * dt);
    
	  //second argument of the function f
	  //temp1 = sol + a31 * k1 + a32 * k2 
	  temp1 = k1;
	  temp1 *= (coeff.a[2][0]);               
	  temp1.axpy(coeff.a[2][1],k2);
	  temp1 += sol;
    

	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);
	  ai.usmv(-1.0,temp1,rhs); 
	  // aplly control
	  rhs += temp2;
     
	  // temporary vector to store the rhs for the third stage
	  // at present k3 is temp vector
	  k4 = rhs;
    
	  //evaluation of the rest of the right hand side
	  // M * tau*Gamma *(c31*k1+c32*k2)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[2][0]/dt);               //Gammaa*c31*k1
	  temp1.axpy((coeff.c[2][1]/dt),k2);
    
	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the third stage vector
	  k3 = 0.0;
    
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k3,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  solver.apply(j,k3,rhs,reduction);
 
	  //////////////////////////////////////////////// 
	  //////////fourth step of LIRK 4-stages///////////
	  //////////////////////////////////////////////// 
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[3] * dt); 

     
	  // in fourth stage we don't need to calculate rhs b'cause a43 = 0 
	  // and alpha[2] = alpha[3], k3 is released
	  rhs = k4;
    
	  //evaluation of the rest of the right hand side
	  // m * (c31*k1+c32*k2+c33*k3)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[3][0]/dt);               //Gammaa*c41*k1
	  temp1.axpy(coeff.c[3][1]/dt,k2);
	  temp1.axpy(coeff.c[3][2]/dt,k3);

	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the third stage vector
	  k4 = 0.0;
    
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k4,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  solver.apply(j,k4,rhs,reduction);
    
	  ////////////////////////////////////////////
	  ////////////// New solution ////////////
	  ////////////////////////////////////////////   
       
	  //sol = (oldsol) + (m1*k1) + (m2*k2)

	  if(adaptive)
		 {   

			temp1 = sol;
			temp1.axpy(coeff.m1[0],k1);
			temp1.axpy(coeff.m1[1],k2); 
			temp1.axpy(coeff.m1[2],k3); 
			temp1.axpy(coeff.m1[3],k4); 
	
			//////the second new solution constructed for adaptive time step
			//rhs = (oldsol) + (m1b*k1) + (m2b*k2)
			rhs = sol;
			rhs.axpy(coeff.m1b[0],k1);
			rhs.axpy(coeff.m1b[1],k2);
			rhs.axpy(coeff.m1b[2],k3);
			rhs.axpy(coeff.m1b[3],k4);
	
			//find the error between two solutions
			rhs -= temp1;
	
			// using L2-norm
			double error = psp.norm(rhs)/gfssize;
	
	
			// bool accept;
			// if(error < tautol)
			//    accept = true;
			// else
			//    accept = false;
	
			// double taunew;
			// // calculate the new time step
			// taunew = fac * dt * std::sqrt( tautol/error );
	
			// if (taunew > facmax * dt)
			//    taunew = facmax * dt;
			// if (taunew < facmin * dt)
			//    taunew = facmin * dt;

			// if(taunew > maxtau)
			//    taunew = maxtau;

			// if(accept)
			//    {
			// 	  sol = temp1;
			// 	  *tau = taunew;
			// 	  return true;
			//    }
			// else 
			//    {
			// 	  *tau = taunew; *err = error;
			// 	  if(gv.comm().rank()==0)
			// 		 std::cout<<" Step is rejcted: time := "<<t<<" tauold = "<<dt <<" taunew := "<<taunew <<" error = "<<error<<std::endl;
			// 	  return false;
			//    }
	
		 }
	  else
		 {
			sol.axpy(coeff.m1[0],k1);
			sol.axpy(coeff.m1[1],k2);
			sol.axpy(coeff.m1[2],k3);
			sol.axpy(coeff.m1[3],k4);
		 }
   
	  return true;
      
   }
   //solver interface 
   template<class GOS, class J, class M, class A, class ODE, class S, class PSP, class V, class PARAM>
   bool odesolve_dual(const GOS& gos, const J& j, const M& m, const A& ai, const ODE& coeff, S& solver, PSP& psp, V& p, V& sol, V& b, RT t, RT *tau, PARAM& params, RT *err) 
   {
	  RT t_temp;
	  V rhs(sol);
	  V k1(sol);
	  V k2(sol);
	  V k3(sol); 
	  V k4(sol);
	  V temp1(sol);
	  V temp2(b);
	  RT dt = (*tau);
	  int ncomp=1; 
  
	  //////////////////////////////////////////////// 
	  //////////First step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////

	  // here alpha_1 = 0  
	  t_temp = t + coeff.alpha[0]*dt;

	  rhs = 0.0;
    
	  // Evaluation of f(.,.) 
	  gos.residual(sol,rhs);       
    
	  ai.umv(sol,rhs);
    
	  // apply rest of the vectors for dual solve
	  //if(lin_dual == false)
	  // {
	  temp2 = 0.0;
	  m.usmv(-1.0,b,temp2);
	  ////////////////////////////m.usmv(params.n2/params.vp,r,temp2);	
	  //ai.umv(p,temp2); 
	  temp2 += p;
	  rhs += temp2;
	  //}

	  //   Dune::PDELab::FilenameHelper fn("monodual_Q1");
	  //     vtkout_2comps(gv,gfs,temp1,rhs,fn.getName());exit(1);
    
	  // initialization for the first vector
	  k1 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k1,rhs,t);
   
	  //solving J*k1 = rhs with Solver and storing the solution in k1
	  //for all stages the left hand side matrix is same
	  //solver.apply(k1,rhs,res1);
	  solver.apply(j,k1,rhs,reduction);
    
	  //temporary vector which is useful in next stage computation
	  temp1 = k1;
    
	  //////////////////////////////////////////////// 
	  //////////second step of LIRK 2-stages///////////
	  ////////////////////////////////////////////////
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[1] * dt);
  
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 *=  coeff.a[1][0];
	  temp1 += sol;
            
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);
	  //params.calculateRHS(gv,mapper,a,m,v,w,temp1,rhs);
	  // // aplly control
	  //if(lin_dual == false)
	  rhs += temp2;
    
	  ai.umv(temp1,rhs);
    
	  //evaluation of the rest of the right hand side:
	  //calculate M*tau*Gammaa*c21/tau*k1
	  temp1 = k1;
	  temp1 *= (coeff.c[1][0]/dt);          //Gammaa*c21*k1
      
	  //needs mass matrix-vector-multiplication
	  //M.umv(temp1,rhs) means rhs = rhs + M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the second stage vector
	  k2 = 0.0;
   
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k2,rhs,t);
   
	  //solving J*k2 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  //solver.apply(k2,rhs,res2);
	  solver.apply(j,k2,rhs,reduction);
   
	  //////////////////////////////////////////////// 
	  //////////third step of LIRK 3-stages///////////
	  ////////////////////////////////////////////////
	  //calculate second argument of f(.,.)
	  // temp1 = oldsol + k1 * a21 
	  temp1 = k1;
	  temp1 *=  coeff.a[2][0]; 
	  temp1.axpy(coeff.a[2][1],k2);
	  temp1 += sol;
    
	  //evaluation of the function f
	  //f = calculate_f (t_temp, temp1);   //f=f(t + alpha2*dt, c + dt*a21*k1)
	  //params.calculateRHS(grid,mapper,grid.leafView(),A,M,temp1,rhs,t_temp);
	  rhs = 0.0; 
	  gos.residual(temp1,rhs);
	  //params.calculateRHS(gv,mapper,a,m,v,w,temp1,rhs);
	  // // aplly control
	  //if(lin_dual == false)
	  rhs += temp2;
    
	  ai.umv(temp1,rhs);
      
	  // temporary vector to store the rhs for the third stage
	  // at present k3 is temp vector
	  k4 = rhs;
      
	  //evaluation of the rest of the right hand side
	  // M * tau*Gamma *(c31*k1+c32*k2)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[2][0]/dt);               //Gammaa*c31*k1
	  temp1.axpy((coeff.c[2][1]/dt),k2);
    
	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the third stage vector
	  k3 = 0.0;
     
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k3,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  solver.apply(j,k3,rhs,reduction);
    

	  //////////////////////////////////////////////// 
	  //////////fourth step of LIRK 4-stages///////////
	  //////////////////////////////////////////////// 
    
	  //update the time for the 2-stage
	  t_temp = t + (coeff.alpha[3] * dt);  

	  // in fourth stage we don't need to calculate rhs b'cause a43 = 0 
	  // and alpha[2] = alpha[3], k3 is released
	  rhs = k4;
    
	  //evaluation of the rest of the right hand side
	  // m * (c31*k1+c32*k2+c33*k3)/tau
	  temp1 = k1;
	  temp1 *= (coeff.c[3][0]/dt);               //Gammaa*c41*k1
	  temp1.axpy(coeff.c[3][1]/dt,k2);
	  temp1.axpy(coeff.c[3][2]/dt,k3);
      
	  //needs matrix-vector-multiplication
	  //M.usmv(-1.0,temp1,rhs) means rhs=rhs - M*temp1
	  m.usmv(-1.0,temp1,rhs);
    
	  //initialization for the third stage vector
	  k4 = 0.0;
    
	  //assemble boundary conditions before calculate solution
	  //asd.AssembleDirichletVec(grid,mapper,k4,rhs,t);
    
	  //solving J*k3 = rhs with Solver and storing the solution in k2
	  //for all stages the left hand side matrix is same
	  solver.apply(j,k4,rhs,reduction);
    
	  ////////////////////////////////////////////
	  ////////////// New solution ////////////
	  ////////////////////////////////////////////   
    
	  //sol = (oldsol) + (m1*k1) + (m2*k2)
	  sol.axpy(coeff.m1[0],k1);
	  sol.axpy(coeff.m1[1],k2);
	  sol.axpy(coeff.m1[2],k3);
	  sol.axpy(coeff.m1[3],k4);
	  return true;
   }              


private:
   // import stuff from PDE solver
   const GV& gv;
   // Solver& solver;
    unsigned gfssize; 

public:
   //coefficients
   typedef Dune::FieldVector<double,4> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,4,4> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,5> LocalVector5;
   LocalVector5 gamma;  
    
    double maxtau;
    bool first_iteration;

};


//ROS3PL
template<class GV>
class ROS3PL : public LIRK4stages<GV>
{
public:
   //constructor
    ROS3PL (const GV& gv_,double maxtau_,unsigned gfssize_) : LIRK4stages<GV>(gv_,maxtau_,gfssize_)
   { 
	  gamma[0] = 0.435866521508459; 
	  gamma[1] = 0.435866521508459; 
	  gamma[2] =-0.064133478491541; 
	  gamma[3] = 0.111028172512505;
	  gamma[4] = 0.000000000000000;
	
	  alpha[0] = 0.000000000000000; 
	  alpha[1] = 0.500000000000000; 
	  alpha[2] = 1.000000000000000;
	  alpha[3] = 1.000000000000000;
	
	  m1[0] = 2.463070773030053;   
	  m1[1] = 1.147140180139521;
	  m1[2] = 0.000000000000000;
	  m1[3] = 1.000000000000000;
	
	  m1b[0] =  2.346947683513665; 
	  m1b[1] =  0.456530569451895;  
	  m1b[2] =  0.056949243945495;
	  m1b[3] =  0.738684936166224;
	
	  a[0][0] = 0.000000000000000;  a[0][1] = 0.000000000000000; 
	  a[0][2] = 0.000000000000000;  a[0][3] = 0.000000000000000;
	
	  a[1][0] = 1.147140180139521;  a[1][1] = 0.000000000000000; 
	  a[1][2] = 0.000000000000000;  a[1][3] = 0.000000000000000;

	  a[2][0] = 2.463070773030053;  a[2][1] = 1.147140180139521; 
	  a[2][2] = 0.000000000000000;  a[2][3] = 0.000000000000000;

	  a[3][0] = 2.463070773030053;  a[3][1] = 1.147140180139521; 
	  a[3][2] = 0.000000000000000;  a[3][3] = 0.000000000000000;
	  
	  c[0][0] = 2.294280360279042;  c[0][1] = 0.000000000000000; 
	  c[0][2] = 0.000000000000000;  c[0][3] = 0.000000000000000; 
   
	  c[1][0] = 2.631861185781065;  c[1][1] = 2.294280360279042; 
	  c[1][2] = 0.000000000000000;  c[1][3] = 0.000000000000000;

	  c[2][0] = 1.302364158113095;  c[2][1] =-2.769432022251304; 
	  c[2][2] = 2.294280360279042;  c[2][3] = 0.000000000000000;

	  c[3][0] = 1.552568958732400;  c[3][1] =-2.587743501215153; 
	  c[3][2] = 1.416993298352020;  c[3][3] = 2.294280360279042;
   }

public:
   //coefficients
   typedef Dune::FieldVector<double,4> LocalVector;
   LocalVector alpha,m1,m1b;	   
   typedef Dune::FieldMatrix<double,4,4> LocalMatrix;
   LocalMatrix a,c;
   typedef Dune::FieldVector<double,5> LocalVector5;
   LocalVector5 gamma;
};
