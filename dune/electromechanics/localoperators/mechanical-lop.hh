#ifndef DUNE_ELASTICITY_MECHANICAL_LOP_HH
#define DUNE_ELASTICITY_MECHANICAL_LOP_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

namespace Dune {
  namespace PDELab {
    namespace Elasticity {

//////////// Finite elasticity LOP /////////////////////////////////////////
template<class GFSV, class V, typename FEMU, typename FEMP, int degreeu, int degreep, int nu, int np, int mvol, int mbnd, int faces>
    class FiniteElasticityLOP :
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::NumericalJacobianVolume<FiniteElasticityLOP<GFSV,V,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces> >,
      public Dune::PDELab::NumericalJacobianApplyVolume<FiniteElasticityLOP<GFSV,V,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces> >
    {
      // define useful types
      typedef typename FEMU::Traits::FiniteElementType FiniteElementTypeU;
      typedef typename FiniteElementTypeU::Traits::LocalBasisType LocalBasisTypeU;
      typedef typename FEMP::Traits::FiniteElementType FiniteElementTypeP;
      typedef typename FiniteElementTypeP::Traits::LocalBasisType LocalBasisTypeP;
      typedef typename LocalBasisTypeU::Traits::DomainType DomainType;
      typedef typename LocalBasisTypeU::Traits::DomainFieldType DF;
      typedef typename LocalBasisTypeU::Traits::RangeType RangeType;
      typedef typename LocalBasisTypeU::Traits::RangeFieldType RF;
      typedef typename LocalBasisTypeU::Traits::JacobianType JacobianType;

      enum {dim=LocalBasisTypeU::Traits::dimDomain};

      RF time = 0.0; // guess what
      const GFSV& gfsv;
      const V& vg;

      // quadrature rules
      Dune::FieldVector<DF,dim> qpvol[mvol];   // quadrature points on volume
      RF wvol[mvol];                           // quadrature weight on refelem
      // Dune::FieldVector<DF,dim-1> qpbnd[mbnd]; // quadrature points on boundary
      // RF wbnd[mbnd];                           // quadrature weight on refelem

      // evaluations of basis functions on the reference element at quadrature points
      RF phihatvol[nu][mvol];          // velocity
      RF psihatvol[np][mvol];          // pressure
      RF Ghat[dim][nu][mvol];          // gradients of velocity basis functions (volume only)
      // RF phihatbnd[faces][nu][mbnd];   // velocity

      // // coefficient functions
      //const GetV& getv;

      // physical parameters
      // *Check the parameter value in ionic model 
      RF v_min = -80.0, v_max = 40.0;
      // RF Sigma_l = 1e-3, Sigma_t = 1e-3;
      RF Sigma_l = 3.28e-1, Sigma_t = 6.99e-2; //gil and git
      RF mu1 = 4.0, beta = 0.3;
      const Dune::FieldMatrix<RF,1,dim> al = {{1.0,0.0}}; // Fiber direction longitudinal
      const Dune::FieldMatrix<RF,1,dim> at = {{0.0,1.0}}; // Fiber direction transverse

    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doLambdaVolume = true };
      enum { doAlphaVolume = true };

      typedef typename Dune::PDELab::LocalFunctionSpace<GFSV> LFSV;
      typedef typename LFSV::template Child<0>::Type LFSV0;
      typedef typename LFSV0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType V0RF;
      typedef typename LFSV0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeT;
      // constructor
      template<typename GV>
      FiniteElasticityLOP (const GV& gv, const GFSV& gfsv_,const V& vg_, const FEMU& femu, const FEMP& femp)
        : gfsv(gfsv_), vg(vg_)
      {
        // // The idea is to do all computations on the reference element only once here in the constructor.
        // // This implies we assume all elements are of the same type, e.g. simplices,
        // // and we use the first element of the mesh as a template.
        // // Moreover, we want to use compile-time known loop bounds. This is accomplished
        // // by checking that the run-time given basius coomplis with the compile-time
        // // given numbers.
        
        // // get finite element basis using given rpresentative element
        auto felu = femu.find(*gv.template begin<0>());
        auto felp = femp.find(*gv.template begin<0>());

        // // check size of the supplied basis
        if (felu.localBasis().size()!=nu) {
          std::cout << "Basis size mismatch for velocity!" << std::endl;
          exit(1);
        }
        if (felp.localBasis().size()!=np) {
          std::cout << "Basis size mismatch for pressure!" << std::endl;
          exit(1);
        }

        // find quadrature rules with the given number of points and
        // evaluate quadrature rules
        Dune::GeometryType gt = felu.type();
        //std::cout << "New Navier-Stokes LOP on " << gt << std::endl;
        int ordervol=-1; // volume rule
        // for (int order=1; order<=20; order++)
        //   {
        //     const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,order);
        //     std::cout << " order " << order << "qpoints " << rule.size() << std::endl;
        //   }
        // exit(1);
        for (int order=20; order>=1; order--)
          {
            const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,order);
            if (rule.size()==mvol)
              {
                ordervol = order;
                //std::cout << "order of volume quadrature with " << mvol << " points is " << order << std::endl;
                for (int i=0; i<mvol; i++) {
                  qpvol[i] = rule[i].position();
                  wvol[i] = rule[i].weight();
                }
                break;
              }
          }
        if (ordervol<0) {
          std::cout << "Could not find volume quadruture rule with that many points!" << std::endl;
          exit(1);
        }
        // int orderbnd=-1; // boundary rule
        // Dune::GeometryType gtbnd;
        // if (gt.isSimplex()) gtbnd = Dune::GeometryTypes::simplex(dim-1);
        // if (gt.isCube()) gtbnd = Dune::GeometryTypes::cube(dim-1);
        // for (int order=20; order>=1; order--)
        //   {
        //     const Dune::QuadratureRule<RF,dim-1>& rule = Dune::QuadratureRules<RF,dim-1>::rule(gtbnd,order);
        //     if (rule.size()==mbnd)
        //       {
        //         orderbnd = order;
        //         //std::cout << "order of boundary quadrature with " << mbnd << " points is " << order << std::endl;
        //         for (int i=0; i<mbnd; i++) {
        //           qpbnd[i] = rule[i].position();
        //           wbnd[i] = rule[i].weight();
        //         }
        //         break;
        //       }
        //   }
        // if (orderbnd<0) {
        //   std::cout << "Could not find boundary quadruture rule with that many points!" << std::endl;
        //   exit(1);
        // }

        // evaluate basis functions on reference element for velocity and pressure
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> phi(nu);
            felu.localBasis().evaluateFunction(qpvol[k],phi);
            for (int i=0; i<nu; i++) phihatvol[i][k] = phi[i];
          }
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> psi(np);
            felp.localBasis().evaluateFunction(qpvol[k],psi);
            for (int i=0; i<np; i++) psihatvol[i][k] = psi[i];
          }
        // for (const auto& is : intersections(gv,*gv.template begin<0>()))
        //   {
        //     auto fgeo_self = is.geometryInInside();
        //     int face = is.indexInInside();
        //     //std::cout << "face " << face << std::endl;
        //     for (int k=0; k<mbnd; k++) // loop over face quadrature points
        //       {
        //         auto qp = fgeo_self.global(qpbnd[k]);
        //         std::vector<RangeType> phi(nu);
        //         felu.localBasis().evaluateFunction(qp,phi);
        //         for (int i=0; i<nu; i++) phihatbnd[face][i][k] = phi[i];
        //       }
        //   }
        
        // evaluate gradients of basis functions on reference element
        std::vector<JacobianType> js(nu);
        for (int k=0; k<mvol; k++)
          {
            felu.localBasis().evaluateJacobian(qpvol[k],js);
            for (int i=0; i<nu; i++)
              for (int j=0; j<dim; j++)
                Ghat[j][i][k] = js[i][0][j];
          }
      }

      
      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSTrial, typename X, typename LFSTest, typename Residual>
      void alpha_volume (const EG& eg, const LFSTrial& lfstrial, const X& x, const LFSTest& lfstest, Residual& residual) const
      {
      //   // extract local function spaces
        const auto& lfsu = lfstrial.template child<0>();    // velocity node
        const auto& lfsp = lfstrial.template child<1>();  // pressure

	typedef Dune::PDELab::LFSIndexCache<LFSV> LFSVCache;
	typedef typename V::template ConstLocalView<LFSVCache> VView;
	LFSV lfsv(gfsv);
	const auto& lfsv0 = lfsv.template child<0>(); 
	int vsize = lfsv0.size();
	LFSVCache lfsvCache(lfsv);
	VView xview(vg);

	std::vector<typename V::ElementType> y(lfsv.maxSize());
	lfsv.bind(eg.entity());
	lfsvCache.update();
	xview.bind(lfsvCache);
	xview.read(y);
        xview.unbind();
std::cout<<"y0 = "<<y[0]<<std::endl;
	// extract degrees of freedom
        RF ZT[dim][nu]; // velocity degrees of freedom
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j) ZT[i][j] = x(lfsu_i,j);
          }
        RF YT[np]; // pressure degrees of freedom
        for (size_t j=0; j<np; ++j) YT[j] = x(lfsp,j);
       
        // evaluate geometry transformation at all quadrature points
        auto geo = eg.geometry();
        RF JT[dim][dim][mvol]; // Jacobian inverse transposed for gradient transformations
        RF factor[mvol];       // quadrature weight times determinant
        
        if (geo.type().isSimplex())
          {
            const auto S = geo.jacobianInverseTransposed(qpvol[0]);
            auto integrationelement = geo.integrationElement(qpvol[0]);
            for (size_t alpha=0; alpha<mvol; ++alpha)
                factor[alpha] = integrationelement*wvol[alpha];
            for (size_t i=0; i<dim; ++i)
              for (size_t j=0; j<dim; ++j)
                for (size_t alpha=0; alpha<mvol; ++alpha)
                  JT[i][j][alpha] = S[i][j];
          }
        if (geo.type().isCube())
          for (size_t alpha=0; alpha<mvol; ++alpha)
            {
              const auto S = geo.jacobianInverseTransposed(qpvol[alpha]);
              auto integrationelement = geo.integrationElement(qpvol[alpha]);
              factor[alpha] = integrationelement*wvol[alpha];
              for (size_t i=0; i<dim; ++i)
                for (size_t j=0; j<dim; ++j)
                  JT[i][j][alpha] = S[i][j];
            }

        // compute gradient of velocity basis functions at all quadrature points
        RF GT[dim][nu][mvol] = {{{0.0}}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<nu; ++j)
            for (size_t k=0; k<dim; ++k) // reduction over k
              for (size_t alpha=0; alpha<mvol; ++alpha)
                GT[i][j][alpha] += JT[i][k][alpha]*Ghat[k][j][alpha];

        // Now evaluate velocity, velocity gradient and pressure at quadrature points
        // Here it turns out to be good idea to have the quadrature loop inside!
        
        // compute velocity at all quadrature points
        RF UT[dim][mvol] = {{0.0}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<nu; ++j) // reduction over j
            for (size_t alpha=0; alpha<mvol; ++alpha)
              UT[i][alpha] += ZT[i][j]*phihatvol[j][alpha];

        // compute velocity gradient at all quadrature points
        RF GradUT[dim][dim][mvol] = {{{0.0}}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<dim; ++j)
            for (size_t k=0; k<nu; ++k) // reduction over k
              for (size_t alpha=0; alpha<mvol; ++alpha)
                GradUT[i][j][alpha] += ZT[i][k]*GT[j][k][alpha];

        // Compute deformation gradient F and det(F)
        RF J[mvol] = {0.0}; 
        Dune::FieldMatrix<RF,dim,dim> F[mvol], FinvT[mvol];
        for (size_t alpha=0; alpha<mvol; ++alpha)
        {
          for (size_t i=0; i<dim; ++i){
             for (size_t j=0; j<dim; ++j){
                F[alpha][i][j] = 1.0*(i==j)+GradUT[i][j][alpha];
             }
          }
          J[alpha] = Dune::FMatrixHelp::invertMatrix_retTransposed(F[alpha],FinvT[alpha]);
        }
        

      //   // compute mu1*J0^n(F^{n+1})(F0^n)^-1(F0^n)^-T
        Dune::FieldMatrix<RF,dim,dim> Inter[mvol];
        {
          RF J0 = 0.0;
          Dune::FieldMatrix<RF,dim,dim> F0, F0invT, KronAt, KronAl;
          for (size_t alpha=0; alpha<mvol; ++alpha)
          { 
             //std::cout<<" qpvol = "<<qpvol[0]<<std::endl;
            // auto v = getv(eg.entity(),qpvol[alpha]);//eg.entity().geometry().center());//qpvol[alpha]);
            std::vector<RangeT> phi(vsize);
            lfsv0.finiteElement().localBasis().evaluateFunction(qpvol[alpha],phi);
	    V0RF v = 0.0;
            for (size_t i=0; i<vsize; i++)
		v += y[i]*phi[i];
        
            Dune::FMatrixHelp::multTransposedMatrix(at,KronAt);
            Dune::FMatrixHelp::multTransposedMatrix(al,KronAl);
            RF Gl = -Sigma_l*beta*(v-v_min)/(v_max-v_min+v);
            RF Gt = -Sigma_t*beta*(v-v_min)/(v_max-v_min+v);
            for (size_t i=0; i<dim; ++i){
               for (size_t j=0; j<dim; ++j){
                  F0[i][j] = 1.0*(i==j)+Gl*KronAl[i][j]+Gt*KronAt[i][j];
               }
           }
            J0 = Dune::FMatrixHelp::invertMatrix_retTransposed(F0,F0invT);
            Dune::FieldMatrix<RF,dim,dim> TmpVar;
            Dune::FMatrixHelp::multTransposedMatrix(F0invT,TmpVar);
            Dune::FMatrixHelp::multMatrix(F[alpha],TmpVar,Inter[alpha]);
            Inter[alpha] *= (mu1*J0);
          }
        }
        // compute pressure at all quadratue points
        RF PT[mvol] = {0.0};
        for (size_t j=0; j<np; ++j) // reduction over j
          for (size_t alpha=0; alpha<mvol; ++alpha)
            PT[alpha] += YT[j]*psihatvol[j][alpha];

        // Now come the residual contributions term by term.
        // These are ultimately reductions over the quadrature points but
        // possibly involve some intermediate results beforehand
        RF RUT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]
        RF RPT[np]  = {0.0};       // residual contribution of this element to pressure test function [j]

        // integral mu1*J0^n(F^{n+1})(F0^n)^-1(F0^n)^-T:GradPhi
        {
          RF intermediate[dim][nu][mvol] = {{{0.0}}};
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              for (size_t k=0; k<dim; ++k) // reduction over k
                for (size_t alpha=0; alpha<mvol; ++alpha)
                  intermediate[i][j][alpha] += Inter[alpha][i][k]*GT[k][j][alpha];
          for (size_t i=0; i<dim; ++i) 
            for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
              for (size_t j=0; j<nu; ++j)
                RUT[i][j] += intermediate[i][j][alpha]*factor[alpha];
        }
        // -integral p^{n+1}(F^-T):GradPhi
        {
          RF intermediate[dim][nu][mvol] = {{{0.0}}};
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              for (size_t k=0; k<dim; ++k) // reduction over k
                for (size_t alpha=0; alpha<mvol; ++alpha)
                  intermediate[i][j][alpha] += (PT[alpha]*FinvT[alpha][i][k])*GT[k][j][alpha];
          for (size_t i=0; i<dim; ++i) 
            for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
              for (size_t j=0; j<nu; ++j)
                RUT[i][j] -= intermediate[i][j][alpha]*factor[alpha];
        }
  
        // integral (J^{n+1}-1)q
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t j=0; j<np; ++j) // pressure test functions here!
            RPT[j] -= (J[alpha])*psihatvol[j][alpha]*factor[alpha];

        // accumulate residual contributions for this element
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j)
               residual.accumulate(lfsu_i,j, RUT[i][j]);
          }
        for (size_t j=0; j<np; ++j)
          residual.accumulate(lfsp,j, RPT[j]);
      }
            //! residual contribution of volume source term
      template<typename EG, typename LFSTest, typename Residual>
      void lambda_volume (const EG& eg, const LFSTest& lfstest, Residual& residual) const
      {
        // define types
        using LFSV = typename LFSTest::template Child<0>::Type;   // the velocity node
        using LFSV0 = typename LFSV::template Child<0>::Type;     // first velocity component

        // extract local function spaces
        const auto& lfsv = lfstest.template child<0>();    // displacement node
        const auto& lfsp = lfstest.template child<1>();  // pressure
        
        // evaluate geometry transformation at all quadrature points
        auto geo = eg.geometry();
        RF factor[mvol];       // quadrature weight times determinant
        if (geo.type().isSimplex())
          {
            auto integrationelement = geo.integrationElement(qpvol[0]);
            for (size_t alpha=0; alpha<mvol; ++alpha)
                factor[alpha] = integrationelement*wvol[alpha];
          }
        if (geo.type().isCube())
          for (size_t alpha=0; alpha<mvol; ++alpha)
            {
              auto integrationelement = geo.integrationElement(qpvol[alpha]);
              factor[alpha] = integrationelement*wvol[alpha];
            }

        // Now evaluate velocity
        // Here it turns out to be good idea to have the quadrature loop inside!
        
        // evaluate force at all quadrature points
        RF FT[dim][mvol] = {{0.0}};
        for (size_t alpha=0; alpha<mvol; ++alpha)
          {
            auto x = eg.entity().geometry().global(qpvol[alpha]);
            //auto f = flocal(eg.entity(),qpvol[alpha]);
            //for (size_t i=0; i<dim; ++i)
             FT[0][alpha] = 0.0;
             FT[1][alpha] = 0.0;
          }

        // Now come the residual contributions term by term.
        RF RT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]
        RF RPT[np]  = {0.0};       // residual contribution of this element to pressure test function [j]

        // mass residual
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              RT[i][j] += -1.0*mu1*FT[i][alpha]*phihatvol[j][alpha]*factor[alpha];

        // integral (1)q
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t j=0; j<np; ++j) // pressure test functions here!
            RPT[j] += psihatvol[j][alpha]*factor[alpha];

        // accumulate residual contributions for this element
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsv_i = lfsv.child(i);
            for (size_t j=0; j<nu; ++j)
              residual.accumulate(lfsv_i,j,-RT[i][j]);
          }
        for (size_t j=0; j<np; ++j)
          residual.accumulate(lfsp,j, RPT[j]);  

      }

    };
    }
  }
}
#endif
