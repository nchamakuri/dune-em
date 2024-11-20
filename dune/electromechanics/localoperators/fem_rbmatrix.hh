#ifndef DUNE_BIDOMAIN_FEM_RBMATRIX_HH
#define DUNE_BIDOMAIN_FEM_RBMATRIX_HH

#include<dune/pdelab/finiteelement/localbasiscache.hh>

namespace Dune {
   namespace PDELab {
	  namespace BiDomain {
		 /** a local operator for solving the equation
		  *
		  *   - \Delta u + au   = f   in \Omega
		  *    \nabla u \cdot n = 0   on \partial\Omega
		  *
		  * with conforming finite elements on all types of grids in any dimension
		  *
		  */  



		 template<class P, class ODE, class BoolMat, int ncomp>
		 class RBBidJacobianOperator : 
			//	public NumericalJacobianVolume<RBBidJacobianOperator<P,ODE,ncomp> >,
			//public Dune::PDELab::FullVolumePattern,
			public Dune::PDELab::LocalOperatorDefaultFlags
		 {
		 public:
			// pattern assembly flags
			enum { doPatternVolume = true };

			// residual assembly flags
			enum { doAlphaVolume = true };

             // define sparsity pattern of operator representation
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

			RBBidJacobianOperator(const P& params_,const ODE& op_, const BoolMat& bm_, unsigned int intorder_=2)
			   : param(params_), op(op_), bm(bm_), intorder(intorder_)
		    {}

  
			// volume integral depending on test and ansatz functions
			template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
			void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
		    {
   
			   // we assume all solution components are the same
			   typedef typename LFSU::template Child<0>::Type LFSV1;
			   const LFSV1& lfsv1 = lfsu.template child<0>();
   
 
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
    
			   // dimensions
			   //const int dim = EG::Geometry::dimension;
               const int dim = EG::Entity::dimension;
			   //const int dimw = EG::Geometry::dimensionworld;

			   int vsize = lfsv1.size();
    
			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 const Dune::FieldVector<DF,dim> local = it->position();
	 
					 // evaluate basis functions
					 std::vector<RangeType> phi(vsize);
					 lfsv1.finiteElement().localBasis().evaluateFunction(local,phi);
	
					 // evaluate u
					 Dune::FieldVector<RF,ncomp> v(0.0);
					 for (size_type i=0; i<vsize; i++)
						for (size_type j=0; j<ncomp; j++)
						   v[j] += x(lfsu.child(j),i)*phi[i];
				
					 // evaluate ionic/reaction term
					 Dune::FieldVector<RF,ncomp> rhs(0.0);
					 param.rhs(eg.geometry().global(local),v,rhs);
	
					 // integrate (K grad u)*grad phi_i + a_0*u*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
					 for (size_type i=0; i<vsize; i++)
						{
						   RF value = phi[i]*factor;
						   for(int rcomp=0;rcomp<ncomp;rcomp++)
							  r.accumulate(lfsv.child(rcomp),i, rhs[rcomp] * value);
						}
				  }
			}
  
			// jacobian of volume term 
			template<typename EG, typename LFSU, typename X, typename LFSV, typename  LocalMatrix>
			void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
								  LocalMatrix& mat) const
		    {
	
			   // we assume all solution components are the same
			   typedef typename LFSU::template Child<0>::Type LFSV1;
			   const LFSV1& lfsv1 = lfsu.template child<0>();
    
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
    
			   // dimensions
               const int dim = EG::Entity::dimension;
			   //const int dimw = EG::Geometry::dimensionworld;

			   int vsize = lfsv1.size();
    
			   // select quadrature rule
			   Dune::GeometryType gt = eg.geometry().type();
			   const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
			   // loop over quadrature points
			   for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
				  {
					 const Dune::FieldVector<DF,dim> local = it->position();
	 
					 // evaluate basis functions
					 std::vector<RangeType> phi(vsize);
					 lfsv1.finiteElement().localBasis().evaluateFunction(local,phi);
	
					 // evaluate u
					 Dune::FieldVector<RF,ncomp> v(0.0);
					 for (size_type i=0; i<vsize; i++)
						for (size_type j=0; j<ncomp; j++)
							  v[j] += x(lfsu.child(j),i)*phi[i];
					 
					 Dune::FieldMatrix<RF,ncomp,ncomp> res(0.0);
					 param.jacobian(eg.geometry().global(local),v,res);

					 // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
					 RF factor = it->weight() * eg.geometry().integrationElement(it->position());
	
					 for (size_type i=0; i<vsize; i++)
						for (size_type j=0; j<vsize; j++)
						   {
							  RF mvalue = phi[j]*phi[i]*factor;
							  for(int rcomp=0;rcomp<ncomp;rcomp++)
								 for(int ccomp=0;ccomp<ncomp;ccomp++)
                                     if(bm[rcomp][ccomp])
                                         mat.accumulate(lfsv.child(rcomp),i,lfsv.child(ccomp),j,-res[rcomp][ccomp]*mvalue);
						   }
				  }
		    }

		 private:
			const P& param;
			const ODE& op;
             const BoolMat& bm;
			unsigned int intorder;
		 };



	  }// namespace BIDOPTIM
   } // namespace PDELab
} // namespace Dune
#endif
