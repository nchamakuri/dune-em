#ifndef DUNE_BIDOMAIN_FEM_STIFFNESSMATRIX_INTRA_HH
#define DUNE_BIDOMAIN_FEM_STIFFNESSMATRIX_INTRA_HH

namespace Dune {
   namespace PDELab {
	  namespace BiDomain {


	    // this iss needed in extracellular potential calcullation for intra cellular stiff matrix
	    // due to scaling in ionic model
	 template<class P0GFS, class U0,class P,int ncomp>
		 class StiffnessMatrixLocalOperator_intra_unscaling : 
			public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };
		   typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS> P0LFSU;
		   //typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
		   //typedef typename P0LFSU::template Child<0>::Type TLFSU;
		   typedef typename P0LFSU::template Child<0>::Type TLFSU;
		   typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
		   
		    
  
			StiffnessMatrixLocalOperator_intra_unscaling (const P0GFS& p0gfs_, const U0& sigma_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

  
  
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // we assume all solution components are the same
			   typedef typename LFSU::template Child<0>::Type LFSU1;
			   const LFSU1& lfsu1 = lfsu.template child<0>();
			   typedef typename LFSV::template Child<0>::Type LFSV1;
			   const LFSV1& lfsv1 = lfsv.template child<0>();
		    
			   // domain and range field type
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSV1::Traits::SizeType size_type;
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV1::Traits::FiniteElementType> FESwitch;
			   // dimensions
			   const int dim = EG::Entity::dimension;
			   //const int dimw = EG::Geometry::dimensionworld;

			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
		    
			   const int usize = lfsu1.size();
			   const int vsize = lfsv1.size();
			   
			   typename P::Traits::TensorType tensor_i;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			   params.Di(eg.entity(),localcenter,tensor_i);
			    
			   
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
			  
			       DF mask_factor = y[0];	 
			       if(mask_factor>0.5)
				 mask_factor = 3.0; //0.8 //0.6
			       // else
			       // 	 mask_factor = 1.0; // in this case no multiplication is required
			       
			       if(mask_factor > 0.0)
				 for(int i=0;i<dim;i++)
				   for(int j=0;j<dim;j++)
				     tensor_i[i][j] *= mask_factor;
			     }                     
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(usize);
					 lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js);
					 //Dune::FieldMatrix<DF,dimw,dim> jac;
					 typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(usize);
					 for (size_type i=0; i<usize; i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
					 //tensor = tensor_i;
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(usize);
					 for (size_type i=0; i<usize; i++)
						tensor_i.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<usize; j++)
						for (size_type i=0; i<usize; i++) 
						   mat.accumulate(lfsv1,i,lfsv1,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const P& params;
			unsigned int intorder;
			bool hetero;
		 };



		 template<class P0GFS, class U0,class P,int ncomp>
		 class StiffnessMatrixLocalOperator_intra : 
              //public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };
             template<typename LFSU, typename LFSV>
             void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                                  LocalSparsityPattern& pattern) const
             {
                 int vsize = lfsv.child(0).size();
                 int usize = lfsu.child(0).size();
                 //std::cout<<" sizes "<<vsize<<" "<<usize<<" "<<lfsv.child(0).size()<<" "<<lfsu.child(0).size()<<std::endl;
                 for (size_t i=0; i<vsize; ++i)
                     for (size_t j=0; j<usize; ++j)
                         pattern.addLink(lfsv.child(0),i,lfsu.child(0),j);
             }
   
             
		     typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS> P0LFSU;
		     //typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
			//typedef typename P0LFSU::template Child<0>::Type TLFSU;
             typedef typename P0LFSU::template Child<0>::Type TLFSU;
             typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
             //	typedef typename P0LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;
  
			StiffnessMatrixLocalOperator_intra (const P0GFS& p0gfs_, const U0& sigma_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

  

             // jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // we assume all solution components are the same
			   typedef typename LFSU::template Child<0>::Type LFSU1;
			   const LFSU1& lfsu1 = lfsu.template child<0>();
			   typedef typename LFSV::template Child<0>::Type LFSV1;
			   const LFSV1& lfsv1 = lfsv.template child<0>();
		    
			   // domain and range field type
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSV1::Traits::SizeType size_type;
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV1::Traits::FiniteElementType> FESwitch;
			   // dimensions
               const int dim = EG::Entity::dimension;
			   //const int dimw = EG::Geometry::dimensionworld;

			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
		    
			   const int usize = lfsu1.size();
			   const int vsize = lfsv1.size();
			   
			   typename P::Traits::TensorType tensor_i, tensor;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			   params.Di(eg.entity(),localcenter,tensor_i);
			    
	       
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(usize);
					 lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js);
					 //Dune::FieldMatrix<DF,dimw,dim> jac;
                     typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(usize);
					 for (size_type i=0; i<usize; i++)
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
					 tensor = tensor_i;
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(usize);
					 for (size_type i=0; i<usize; i++)
						tensor.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<usize; j++)
						for (size_type i=0; i<usize; i++) 
						   mat.accumulate(lfsv1,i,lfsv1,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const P& params;
			unsigned int intorder;
			bool hetero;
		 };

		template<class P0GFS, class U0,class P, typename GetF, int ncomp>
                class StiffnessMatrixLocalOperator_mechano_intra :
              	//public Dune::PDELab::FullVolumePattern,
              	 public Dune::PDELab::LocalOperatorDefaultFlags
                 {
                 public:

                        enum { doPatternVolume = true }; // pattern assembly flags
                        enum { doAlphaVolume = true }; //residual assembly flags

                        template<typename LFSU, typename LFSV>
                        void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                                  LocalSparsityPattern& pattern) const
                        {
                                int vsize = lfsv.child(0).size();
                                int usize = lfsu.child(0).size();
                                for (size_t i=0; i<vsize; ++i)
                                        for (size_t j=0; j<usize; ++j)
                                                pattern.addLink(lfsv.child(0),i,lfsu.child(0),j);
                        }

                        typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS> P0LFSU;
                        typedef typename P0LFSU::template Child<0>::Type TLFSU;
                        typedef typename TLFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType P0RF;

                        StiffnessMatrixLocalOperator_mechano_intra (const P0GFS& p0gfs_, const U0& sigma_, const P& params_,GetF& getF_, unsigned int intorder_=2, bool hetero_=false)
                           : p0gfs(p0gfs_),sigma(sigma_), params(params_), getF(getF_), intorder(intorder_),hetero(hetero_)
                    {}

			// jacobian of volume term
                        template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
                        void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                                                                  LocalMatrix& mat) const
                    {
                           // we assume all solution components are the same
                           typedef typename LFSU::template Child<0>::Type LFSU1;
                           const LFSU1& lfsu1 = lfsu.template child<0>();
                           typedef typename LFSV::template Child<0>::Type LFSV1;
                           const LFSV1& lfsv1 = lfsv.template child<0>();

                            //domain and range field type
                           typedef typename LFSV1::Traits::FiniteElementType::
                                  Traits::LocalBasisType::Traits::DomainFieldType DF;
                           typedef typename LFSV1::Traits::FiniteElementType::
                                  Traits::LocalBasisType::Traits::RangeFieldType RF;
                           typedef typename LFSV1::Traits::FiniteElementType::
                                  Traits::LocalBasisType::Traits::JacobianType JacobianType;
                           typedef typename LFSV1::Traits::FiniteElementType::
                                  Traits::LocalBasisType::Traits::RangeType RangeType;
                           typedef typename LFSV1::Traits::SizeType size_type;
                            //domain and range field type
                           typedef FiniteElementInterfaceSwitch<typename LFSV1::Traits::FiniteElementType> FESwitch;
                           // dimensions
                           const int dim = EG::Entity::dimension;
                           //const int dimw = EG::Geometry::dimensionworld;

                           //select quadrature rule
                           Dune::GeometryType gt = eg.geometry().type();
                           const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

                           const int usize = lfsu1.size();
                           const int vsize = lfsv1.size();

                           typename P::Traits::TensorType tensor_i, tensor;
                           Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
                           params.Di(eg.entity(),localcenter,tensor_i);

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
                           DF mask_factor = y[0];
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
                                Dune::FieldMatrix<RF,dim,dim> Finv,FinvT;
                           // loop over quadrature points
                           for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                                  {
                                          // evaluate deformation gradient
                                          auto F = getF(eg.entity(),it->position());
                                          // Evaluate (I+GradV)^-1 and (I+GradV)^-T
                                        Dune::FMatrixHelp::invertMatrix(F,Finv);
                                        Dune::FMatrixHelp::invertMatrix_retTransposed(F,FinvT);
                                         // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                                         std::vector<JacobianType> js(usize);
                                         lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js);
                                         //Dune::FieldMatrix<DF,dimw,dim> jac;
                     typename EG::Geometry::JacobianInverseTransposed jac;
                                         // transform gradient to real element
                                         jac = eg.geometry().jacobianInverseTransposed(it->position());
                                         std::vector<Dune::FieldVector<RF,dim> > gradphi(usize);
                                        for (size_type i=0; i<usize; i++)
                                                {
                                                   gradphi[i] = 0.0;
                                                   jac.umv(js[i][0],gradphi[i]);
                                                }
                                         // do it only when heterogenity is set
                                          if(hetero)
                                                {
                                                   if(mask_factor>=0.0)
                                                          for(int i=0;i<dim;i++)
                                                                 for(int j=0;j<dim;j++)
                                                                        tensor_i[i][j] *= mask_factor;
                                                }
                                         tensor = tensor_i;
                                         // evaluate (F^-1)D(F^-T)
                                        Dune::FieldMatrix<RF,dim,dim> ret1, ret2;
                                        Dune::FMatrixHelp::multMatrix(Finv,tensor,ret1);
                                        Dune::FMatrixHelp::multMatrix(ret1,FinvT,ret2);
                                         // compute K * gradient of shape functions
                                         std::vector<Dune::FieldVector<RF,dim> > Kgradphi(usize);
                                         for (size_type i=0; i<usize; i++)
                                                ret2.mv(gradphi[i],Kgradphi[i]);
                                                // tensor.mv(gradphi[i],Kgradphi[i]);

                                         // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
                                         RF factor = it->weight() * eg.geometry().integrationElement(it->position());

                                         for (size_type j=0; j<usize; j++)
                                                for (size_type i=0; i<usize; i++)
                                                   mat.accumulate(lfsv1,i,lfsv1,j, (Kgradphi[j]*gradphi[i])*factor);
                                  }
                    }

                 private:
                        const P0GFS& p0gfs;
                        const U0& sigma;
                        const GetF& getF;
                        const P& params;
                        unsigned int intorder;
                        bool hetero;
                 };






	 template<class P0GFS, class U0,class P,int ncomp>
		 class StiffnessMatrixLocalOperator_sigma_intra : 
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
  
			StiffnessMatrixLocalOperator_sigma_intra (const P0GFS& p0gfs_, const U0& sigma_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

  
  
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // we assume all solution components are the same
			   typedef typename LFSU::template Child<0>::Type LFSU1;
			   const LFSU1& lfsu1 = lfsu.template child<0>();
			   typedef typename LFSV::template Child<0>::Type LFSV1;
			   const LFSV1& lfsv1 = lfsv.template child<0>();
		    
			   // domain and range field type
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSV1::Traits::SizeType size_type;
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV1::Traits::FiniteElementType> FESwitch;
			   // dimensions
               const int dim = EG::Entity::dimension;
			   			   
			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
		    
			   const int usize = lfsu1.size();
			   const int vsize = lfsv1.size();
		    
			   typename P::Traits::TensorType tensor_i;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			  
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
			   // P0LFSU p0lfsu(p0gfs);
			   // p0lfsu.bind(eg.entity());
			   // LocalVector<typename U0::ElementType, TrialSpaceTag> yy(p0lfsu.size());
			   // // read coefficients
			   // p0lfsu.vread(sigma,yy);
			   
			   Dune::FieldVector<DF,dim> y(0.0); 
			   for (size_type i=0; i<dim; i++)
				  y[i] = yy(p0lfsu.child(i),0);
			   
			   //params.Di(eg.entity(),localcenter,y,tensor_i); // for incorporating the sigma's
               params.Di(eg.entity(),localcenter,tensor_i);
			   
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(usize);
					 lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js);
                     typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(usize);
					 for (size_type i=0; i<usize; i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(usize);
					 for (size_type i=0; i<usize; i++)
						tensor_i.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<usize; j++)
						for (size_type i=0; i<usize; i++) 
						   mat.accumulate(lfsv1,i,lfsv1,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const P& params;
			unsigned int intorder;
			bool hetero;
		 };

		 
        template<class P0GFS, class U0,class P,int ncomp>
		 class StiffnessMatrixLocalOperator_sigma_fibro_intra : 
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
  
			StiffnessMatrixLocalOperator_sigma_fibro_intra (const P0GFS& p0gfs_, const U0& sigma_, const U0& fibro_, const P& params_,unsigned int intorder_=2, bool hetero_=false)
			   : p0gfs(p0gfs_),sigma(sigma_), fibro(fibro_), params(params_), intorder(intorder_),hetero(hetero_)
		    {}

  
  
			// jacobian of volume term
			template<typename EG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
			   // we assume all solution components are the same
			   typedef typename LFSU::template Child<0>::Type LFSU1;
			   const LFSU1& lfsu1 = lfsu.template child<0>();
			   typedef typename LFSV::template Child<0>::Type LFSV1;
			   const LFSV1& lfsv1 = lfsv.template child<0>();
		    
			   // domain and range field type
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::DomainFieldType DF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeFieldType RF;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::JacobianType JacobianType;
			   typedef typename LFSV1::Traits::FiniteElementType::
				  Traits::LocalBasisType::Traits::RangeType RangeType;
			   typedef typename LFSV1::Traits::SizeType size_type;
			   // domain and range field type
			   typedef FiniteElementInterfaceSwitch<typename LFSV1::Traits::FiniteElementType> FESwitch;
			   // dimensions
               const int dim = EG::Entity::dimension;
			   
			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
		    
			   const int usize = lfsu1.size();
			   const int vsize = lfsv1.size();
		    
			   typename P::Traits::TensorType tensor_i;
			   Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			  
			   // P0LFSU p0lfsu(p0gfs);
			   // p0lfsu.bind(eg.entity());
			   // LocalVector<typename U0::ElementType, TrialSpaceTag> yy(p0lfsu.size()), zz(p0lfsu.size());
			   // // read coefficients
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
			   
			   double rand_value = zz(p0lfsu.child(0),0);
			   if(rand_value<0.4)
				  rand_value = 1E-8;
			   Dune::FieldVector<DF,dim> y(0.0); 
			   for (size_type i=0; i<dim; i++)
				  y[i] = yy(p0lfsu.child(i),0);
			   
			   params.Di(eg.entity(),localcenter,y,tensor_i);
			   tensor_i *= rand_value;
               yview.unbind();
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
					 std::vector<JacobianType> js(usize);
					 lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js);
                     typename EG::Geometry::JacobianInverseTransposed jac;
					 // transform gradient to real element
					 jac = eg.geometry().jacobianInverseTransposed(it->position());
					 std::vector<Dune::FieldVector<RF,dim> > gradphi(usize);
					 for (size_type i=0; i<usize; i++)
						{
						   gradphi[i] = 0.0;
						   jac.umv(js[i][0],gradphi[i]);
						}
	 
					 // compute K * gradient of shape functions
					 std::vector<Dune::FieldVector<RF,dim> > Kgradphi(usize);
					 for (size_type i=0; i<usize; i++)
						tensor_i.mv(gradphi[i],Kgradphi[i]);
	
					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type j=0; j<usize; j++)
						for (size_type i=0; i<usize; i++) 
						   mat.accumulate(lfsv1,i,lfsv1,j, (Kgradphi[j]*gradphi[i])*factor); 
				  }
		    }

		 private:  
			const P0GFS& p0gfs;
			const U0& sigma;
			const U0& fibro;
			const P& params;
			unsigned int intorder;
			bool hetero;
		 };




	  }// namespace BIDOPTIM
   } // namespace PDELab
} // namespace Dune
#endif
