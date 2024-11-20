#ifndef DUNE_BIDOMAIN_FEM_STIFFNESSMATRIX_EXTRA_HH
#define DUNE_BIDOMAIN_FEM_STIFFNESSMATRIX_EXTRA_HH

namespace Dune {
   namespace PDELab {
	  namespace BiDomain {
	
		 template<class P0GFS, class U0,class P>
		 class StiffnessMatrixLocalOperator_extra : 
			public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };
			typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
             typedef typename P0LFSU::template Child<0>::Type TLFSU;
             typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
             //typedef typename P0LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
			// typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
			// typedef typename P0LFSU::template Child<0>::Type TLFSU;
			// typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;

             
			StiffnessMatrixLocalOperator_extra (const P0GFS& p0gfs_, const U0& sigma_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

		   //this inly for MEA setting: rhs vector
		   // volume integral depending on test and ansatz functions
		   template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
		   void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,  R& r) const
		   {
		     // domain and range field type
		     typedef typename LFSU::Traits::FiniteElementType::
		       Traits::LocalBasisType::Traits::DomainFieldType DF;
		     typedef typename LFSU::Traits::FiniteElementType::
		       Traits::LocalBasisType::Traits::RangeFieldType RF;
		     typedef typename LFSU::Traits::FiniteElementType::
		       Traits::LocalBasisType::Traits::JacobianType JacobianType;
		     typedef typename LFSU::Traits::FiniteElementType::
		       Traits::LocalBasisType::Traits::RangeType RangeType;
    
		     typedef typename LFSU::Traits::SizeType size_type;
        
		     // dimensions
		     const int dim = EG::Entity::dimension;
		     // const int dimw = EG::Geometry::dimensionworld;
     
		     typedef FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
    
		     // select quadrature rule
		     Dune::GeometryType gt = eg.geometry().type();
		     const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
		     
		     int counter =0;
		     for (size_type i=0; i<lfsu.size(); ++i)
		       if(x(lfsu,i)>0) // if electrode position
			 counter++;
		     
		     if(counter<3) //if not a MEA element
		       return;
		     
		     // loop over quadrature points
		     for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
		       {
			 // evaluate shape functions 
			 std::vector<RangeType> phi(lfsu.size());
			 FESwitch::basis(lfsu.finiteElement()).evaluateFunction(it->position(),phi);
	
			 // integrate (K grad u)*grad phi_i + a_0*u*phi_i
			 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
			 for (size_type i=0; i<lfsu.size(); i++)
			   r.accumulate(lfsu,i, phi[i]*factor); 
		       }
		   }
  
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // domain and range field type
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSU::Traits::SizeType size_type;
    
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
			   //typedef BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
    
			   // dimensions
			   const int dim = EG::Entity::dimension;
			   

			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
			   typename P::Traits::TensorType tensor_i, tensor_e;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			   params.Di(eg.entity(),localcenter,tensor_i);
			   params.De(eg.entity(),localcenter,tensor_e);
			   
			   // do it only when heterogenity is set 
			   if(hetero)
			     {
			       typedef Dune::PDELab::LFSIndexCache<P0LFSU> LFSU0Cache;
			       typedef typename U0::template ConstLocalView<LFSU0Cache> U0View;
			       P0LFSU p0lfsu(p0gfs);
			       LFSU0Cache lfsu0Cache(p0lfsu);
			       U0View xview(sigma);
			       // local coefficients
			       std::vector<typename U0::ElementType> y(p0lfsu.maxSize());
			       p0lfsu.bind(eg.entity());
			       lfsu0Cache.update();
			       xview.bind(lfsu0Cache);
			       xview.read(y);
			       xview.unbind();
			       
			       DF mask_factor = y[0];//y(p0lfsu,0);	 
			       if(mask_factor>0.5)
				 mask_factor = 3.0; // 0.8 old 0.6
			       else if(mask_factor>=0.0)
				 mask_factor = 1.0; // in this case no multiplication is required
			       else
				 mask_factor = 0.0; // for negative values, i.e in the bath domain

			       
			       //if(mask_factor>1E-5) // only white spots will be affected
				 for(int i=0;i<dim;i++)
				   for(int j=0;j<dim;j++)
				     tensor_i[i][j] *= mask_factor;
			     }
			   //tensor = tensor_i;
			   tensor_i += tensor_e; // i.e. tensor_ie = tensor_i + tensor_e
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(lfsu.size());
					 lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
					 typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
					 for (size_type i=0; i<lfsu.size(); i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
					 // // do it only when heterogenity is set 
					 // if(hetero)
					 // 	{
					 // 	   if(mask_factor>=0.0)
					 // 		  for(int i=0;i<dim;i++)
					 // 			 for(int j=0;j<dim;j++)
					 // 				tensor_i[i][j] *= mask_factor;
					 // 	}
					
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
					 for (size_type i=0; i<lfsu.size(); i++)
						tensor_i.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<lfsu.size(); j++)
						for (size_type i=0; i<lfsu.size(); i++) 
						   mat.accumulate(lfsu,i,lfsu,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const P& params;
			unsigned int intorder;
			bool hetero;
		 };



          template<class P0GFS, class U0,class P>
		 class StiffnessMatrixLocalOperator_extra_tissue : 
			public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };
			typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
             typedef typename P0LFSU::template Child<0>::Type TLFSU;
             typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
         
             
             StiffnessMatrixLocalOperator_extra_tissue (const P0GFS& p0gfs_, const U0& sigma_, const P& params_, std::array<double,4>& t_domain_, unsigned int intorder_=2, bool hetero_=false)
                 : p0gfs(p0gfs_),sigma(sigma_), params(params_), t_domain(t_domain_), intorder(intorder_),hetero(hetero_)
		    {}

		    
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // domain and range field type
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSU::Traits::SizeType size_type;
    
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
			   //typedef BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
    
			   // dimensions
               const int dim = EG::Entity::dimension;
               
               Dune::FieldVector<DF,dim> globalcenter = eg.geometry().center();

               if(globalcenter[0]< t_domain[0] || globalcenter[1]< t_domain[1] || globalcenter[0]>t_domain[2]  || globalcenter[1] >t_domain[3])
                   return;

			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
			   typename P::Traits::TensorType tensor_i, tensor_e, tensor;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			   params.Di(eg.entity(),localcenter,tensor_i);
			   params.De(eg.entity(),localcenter,tensor_e);
			   
			   typedef Dune::PDELab::LFSIndexCache<P0LFSU> LFSU0Cache;
			   typedef typename U0::template ConstLocalView<LFSU0Cache> U0View;
			   P0LFSU p0lfsu(p0gfs);
			   LFSU0Cache lfsu0Cache(p0lfsu);
			   U0View xview(sigma);
			   // local coefficients
			   std::vector<typename U0::ElementType> y(p0lfsu.maxSize());
			   p0lfsu.bind(eg.entity());
			   lfsu0Cache.update();
			   xview.bind(lfsu0Cache);
			   xview.read(y);
			   xview.unbind();
			 
			   DF mask_factor = y[0];//y(p0lfsu,0);	 
			   if(mask_factor>0.5)
			   	  mask_factor = 0.01;
			   else
			   	  mask_factor = 1.0;
			   // do it only when heterogenity is set 
			   if(hetero)
				  {
					 if(mask_factor>=0.0)
						for(int i=0;i<dim;i++)
						   for(int j=0;j<dim;j++)
							  tensor_i[i][j] *= mask_factor;
				  }
               
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(lfsu.size());
					 lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
                     typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
					 for (size_type i=0; i<lfsu.size(); i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
				
					 tensor = tensor_i;
					 tensor += tensor_e;
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
					 for (size_type i=0; i<lfsu.size(); i++)
						tensor.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<lfsu.size(); j++)
						for (size_type i=0; i<lfsu.size(); i++) 
						   mat.accumulate(lfsu,i,lfsu,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const P& params;
            std::array<double,4>& t_domain;
			unsigned int intorder;
			bool hetero;
		 };


          
		 template<typename FiniteElementMap, class P0GFS, class U0,class P>
		 class StiffnessMatrixLocalOperator_sigma_extra : 
			public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };
			//typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
			//typedef typename P0LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
			typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
			typedef typename P0LFSU::template Child<0>::Type TLFSU;
			typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
			
			StiffnessMatrixLocalOperator_sigma_extra (const P0GFS& p0gfs_, const U0& sigma_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

  
  
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // domain and range field type
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSU::Traits::SizeType size_type;
    
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
			   //typedef BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
    
			   // dimensions
               const int dim = EG::Entity::dimension;
			   
			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
			   typename P::Traits::TensorType tensor_i, tensor_e;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			
			   // P0LFSU p0lfsu(p0gfs);
			   // p0lfsu.bind(eg.entity());
			   // LocalVector<typename U0::ElementType, TrialSpaceTag> yy(p0lfsu.size());
			   // // read coefficents
			   // p0lfsu.vread(sigma,yy);
			   typedef Dune::PDELab::LFSIndexCache<P0LFSU> LFSU0Cache;
			   typedef typename U0::template ConstLocalView<LFSU0Cache> U0View;
			   P0LFSU p0lfsu(p0gfs);
			   LFSU0Cache lfsu0Cache(p0lfsu);
			   U0View xview(sigma);
			   // local coefficients
			   LocalVector<typename U0::ElementType, TrialSpaceTag> yy(p0lfsu.maxSize());
			   p0lfsu.bind(eg.entity());
			   lfsu0Cache.update();
			   xview.bind(lfsu0Cache);
			   xview.read(yy);
			   xview.unbind();
			   Dune::FieldVector<DF,dim> y(0.0); 
			   for (size_type i=0; i<dim; i++)
				  y[i] = yy(p0lfsu.child(i),0);
			   //params.Di(eg.entity(),localcenter,y,tensor_i); // for incorporating conductivity tensosrs
               params.Di(eg.entity(),localcenter,tensor_i);
			   params.De(eg.entity(),localcenter,tensor_e);
               Dune::FieldMatrix<double,2,2> tens(0.0);
			   tensor_i += tensor_e;
               //tens = 0.1;
               //std::cout<<" global "<<eg.geometry().global(localcenter)<<std::endl;

			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 // std::vector<JacobianType> js(lfsu.size());
					 // lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
                    
                    const std::vector<JacobianType>& js = cache.evaluateJacobian(it->position(),lfsu.finiteElement().localBasis());

					 // Dune::FieldMatrix<DF,dimw,dim> jac;
					 // // transform gradient to real element
					 // jac = eg.geometry().jacobianInverseTransposed(it->position());
                     const typename EG::Geometry::JacobianInverseTransposed jac
                       = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
                     for (size_type i=0; i<lfsu.size(); i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
                     for (size_type i=0; i<lfsu.size(); i++)
                       tensor_i.mv(gradphi[i],Kgradphi[i]);
                     
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<lfsu.size(); j++)
						for (size_type i=0; i<lfsu.size(); i++) 
                          {
                            double value = (Kgradphi[j]*gradphi[i]);
                            mat.accumulate(lfsv,i,lfsv,j, value*factor);
                            //std::cout<<" value "<<value <<" rank "<<Dune::MPIHelper::getCollectiveCommunication().rank()<<std::endl;
                          //mat.accumulate(lfsv,i,lfsv,j, (Kgradphi[j]*gradphi[i])*factor); 
                          }
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const P& params;
			unsigned int intorder;
			bool hetero;
           typedef typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;
           Dune::PDELab::LocalBasisCache<LocalBasisType> cache;

		 };

	 template<class P0GFS, class U0,class P>
		 class StiffnessMatrixLocalOperator_sigma_fibro_extra : 
			public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };
			//typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
			//typedef typename P0LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
			typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
			typedef typename P0LFSU::template Child<0>::Type TLFSU;
			typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
			
			StiffnessMatrixLocalOperator_sigma_fibro_extra (const P0GFS& p0gfs_, const U0& sigma_, const U0& fibro_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), fibro(fibro_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

  
  
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // domain and range field type
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSU::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSU::Traits::SizeType size_type;
    
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
			   //typedef BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
    
			   // dimensions
               const int dim = EG::Entity::dimension;
			   
			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
			   typename P::Traits::TensorType tensor_i, tensor_e;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			
			   // P0LFSU p0lfsu(p0gfs);
			   // p0lfsu.bind(eg.entity());
			   // LocalVector<typename U0::ElementType, TrialSpaceTag> yy(p0lfsu.size()), zz(p0lfsu.size());;
			   // // read coefficents
			   // p0lfsu.vread(sigma,yy);
			   // p0lfsu.vread(fibro,zz);
			   typedef Dune::PDELab::LFSIndexCache<P0LFSU> LFSU0Cache;
			   typedef typename U0::template ConstLocalView<LFSU0Cache> U0View;
			   P0LFSU p0lfsu(p0gfs);
			   LFSU0Cache lfsu0Cache(p0lfsu);
			   U0View yview(sigma);
			   U0View zview(fibro);
			   // local coefficients
			   LocalVector<typename U0::ElementType, TrialSpaceTag> yy(p0lfsu.maxSize()),zz(p0lfsu.maxSize());
			   p0lfsu.bind(eg.entity());
			   lfsu0Cache.update();
			   yview.bind(lfsu0Cache);zview.bind(lfsu0Cache);
			   yview.read(yy);zview.read(zz);
			   yview.unbind();zview.unbind();
			   double rand_value = zz(p0lfsu.child(0),0);
			   Dune::FieldVector<DF,dim> y(0.0); 
			   for (size_type i=0; i<dim; i++)
				  y[i] = yy(p0lfsu.child(i),0);
			   params.Di(eg.entity(),localcenter,y,tensor_i);
			   params.De(eg.entity(),localcenter,tensor_e);
			 
			   tensor_i *= rand_value;
			   tensor_i += tensor_e;
			   
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(lfsu.size());
					 lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
                     typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
					 for (size_type i=0; i<lfsu.size(); i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
					 for (size_type i=0; i<lfsu.size(); i++)
						tensor_i.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<lfsu.size(); j++)
						for (size_type i=0; i<lfsu.size(); i++) 
						   mat.accumulate(lfsu,i,lfsu,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;const U0& fibro;
			const P& params;
			unsigned int intorder;
			bool hetero;
		 };


	  }// namespace BIDOPTIM
   } // namespace PDELab
} // namespace Dune
#endif
