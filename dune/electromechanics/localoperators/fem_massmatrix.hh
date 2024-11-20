#ifndef DUNE_BIDOMAIN_FEM_MASSMATRIX_HH
#define DUNE_BIDOMAIN_FEM_MASSMATRIX_HH

#include<dune/pdelab/finiteelement/localbasiscache.hh>

namespace Dune {
    namespace PDELab {
        namespace BiDomain {
	
            template<class P,int ncomp>
            class MassMatrixLocalOperator : 
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

                MassMatrixLocalOperator (const P& params_,unsigned int intorder_=2)
                    : params(params_),intorder(intorder_)
                {}

                // volume integral depending on test and ansatz functions
                template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
                void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,  R& r) const
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
    
                
                    typedef FiniteElementInterfaceSwitch<typename LFSV1::Traits::FiniteElementType> FESwitch;
    
                    // select quadrature rule
                    Dune::GeometryType gt = eg.geometry().type();
                    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
              
                    // loop over quadrature points
                    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                        {
                            auto xglobal = eg.geometry().global(it->position());
                            double Istim = params.getStimulusSupport(xglobal);
                            if(Istim < 1E-5)
                                continue;
                      
                            // //if (xglobal[1]<0.021) std::cout<<" value "<<xglobal[1]<<std::endl;
                            // if (xglobal[1]>0.021)
                            //     continue;
                            // evaluate shape functions 
                            std::vector<RangeType> phi(lfsv1.size());
                            FESwitch::basis(lfsv1.finiteElement()).evaluateFunction(it->position(),phi);
                     
                            // integrate (K grad u)*grad phi_i + a_0*u*phi_i
                            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
                            for (size_type i=0; i<vsize; i++)
                                r.accumulate(lfsv1,i, Istim*phi[i]*factor); //r[i] += phi[i]*factor;
                        }
                }
             
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
                    const int usize = lfsu1.size();
                    const int vsize = lfsv1.size();
                    // dimensions
                    const int dim = EG::Entity::dimension;
                    //const int dimw = EG::Geometry::dimensionworld;
	       
       
                    // select quadrature rule
                    Dune::GeometryType gt = eg.geometry().type();
                    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
			   
                    // loop over quadrature points
                    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                        {
                            // evaluate shape functions 
                            std::vector<RangeType> phi(usize);
                            FESwitch::basis(lfsu1.finiteElement()).evaluateFunction(it->position(),phi);
                            // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
                            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
                            for (size_type i=0; i<usize; i++)
                                for (size_type j=0; j<usize; j++)
                                    {
                                        RF value = ( phi[j]*phi[i] )*factor; 
                                        for(int comp=0;comp<ncomp;comp++)
                                            {
                                                const LFSV1& lfssol = lfsv.child(comp);
                                                mat.accumulate(lfssol,i,lfssol,j,value);   
                                            }
                                    }
                        }
			   
			   
                } 
  
            private:  
                const P& params;
                unsigned int intorder;
            };

            // template<class P>
            class MassMatrix_ForPenaltyMethod: 
                public Dune::PDELab::FullVolumePattern,
                public Dune::PDELab::LocalOperatorDefaultFlags
            {
            public:
                // pattern assembly flags
                enum { doPatternVolume = true };

                // residual assembly flags
                enum { doAlphaVolume = true };
	   
                // MassMatrix_ForPenaltyMethod (const P& params_,unsigned int intorder_=2)
                //    : Params(params_),intorder(intorder_)
                MassMatrix_ForPenaltyMethod (unsigned int intorder_=2)
                    : intorder(intorder_)
                {}

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

                    // loop over quadrature points
                    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                        {
                            // evaluate shape functions 
                            std::vector<RangeType> phi(lfsu.size());
                            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(it->position(),phi);
	
                            // integrate (K grad u)*grad phi_i + a_0*u*phi_i
                            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
                            for (size_type i=0; i<lfsu.size(); i++)
                                r.accumulate(lfsu,i, phi[i]*factor); //r[i] += phi[i]*factor;
                        }
                }
	   
                // jacobian of volume term :: do nothing here
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
                    //const int dimw = EG::Geometry::dimensionworld;
    
                    // select quadrature rule
                    Dune::GeometryType gt = eg.geometry().type();
                    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
    
                    // loop over quadrature points
                    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                        {
                            // evaluate shape functions 
                            std::vector<RangeType> phi(lfsu.size());
                            FESwitch::basis(lfsu.finiteElement()).evaluateFunction(it->position(),phi);
                            // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
                            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
                            for (size_type i=0; i<lfsu.size(); i++)
                                for (size_type j=0; j<lfsu.size(); j++)
                                    mat.accumulate(lfsu,i,lfsu,j,( phi[j]*phi[i] )*factor);
                        }
                } 
  
            private:  
                //	const P& Params;
                unsigned int intorder;
            };



        }// namespace BIDOPTIM
    } // namespace PDELab
} // namespace Dune
#endif
