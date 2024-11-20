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


//////////// Finite elasticity LOP /////////////////////////////////////////
template<typename Scheme,class IonicModel,int ncomp> 
    class ElectromechanicalLOP :
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>,
      public Dune::PDELab::NumericalJacobianVolume<ElectromechanicalLOP<Scheme,IonicModel,ncomp> >,
      public Dune::PDELab::NumericalJacobianApplyVolume<ElectromechanicalLOP<Scheme,IonicModel,ncomp> >
    {
      // define useful types
      typedef typename Scheme::FEMU FEMU;
      typedef typename Scheme::FEMP FEMP;
      typedef typename Scheme::FEMV FEMV;
      typedef typename Scheme::FEMW FEMW;
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
      enum {degreeu = Scheme::degreeu};
      enum {degreep = Scheme::degreep};
      enum {degreev = Scheme::degreev};
      enum {degreew = Scheme::degreew};
      enum {nu = Scheme::nu};
      enum {np = Scheme::np};
      enum {nv = Scheme::nv};
      enum {nw = Scheme::nw};
      enum {mvol = Scheme::mvol};
      enum {mbnd = Scheme::mbnd};
      enum {faces = Scheme::faces};

      RF time = 0.0;
      using LocalBasisV = typename FEMV::Traits::FiniteElementType::Traits::LocalBasisType;
      Dune::PDELab::LocalBasisCache<LocalBasisV> cache;

      // quadrature rules
      Dune::FieldVector<DF,dim> qpvol[mvol];   // quadrature points on volume
      RF wvol[mvol];                           // quadrature weight on refelem
      Dune::FieldVector<DF,dim-1> qpbnd[mbnd]; // quadrature points on boundary
      RF wbnd[mbnd];                           // quadrature weight on refelem

      // evaluations of basis functions on the reference element at quadrature points
      RF phihatvol[nu][mvol];          // velocity
      RF psihatvol[np][mvol];          // pressure
      RF zetahatvol[nv][mvol];         // transmemberain potential
      RF etahatvol[nw][mvol];          // gating variable
      RF Ghat[dim][nu][mvol];          // gradients of velocity basis functions (volume only)
      RF phihatbnd[faces][nu][mbnd];   // velocity


      // physical parameters
      RF v_min, v_max;         // Ionic Model minimum and maximum potential
      RF mu1, E, rho;           // Material Parameters
      RF beta;                 // Parameter used in finding/scaling the active stress 
      RF gil, git;             // conductivity factors
      // # al and at are depend of the position so need to moved into a function
      // # Make this part general for all dimension
      const Dune::FieldMatrix<RF,1,dim> al = {{1.0,0.0}}; // Fiber direction longitudinal
      const Dune::FieldMatrix<RF,1,dim> at = {{0.0,1.0}}; // Fiber direction transverse
    
    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doLambdaVolume = true };
      enum { doAlphaVolume = true };

      // constructor
      template<typename GV>
      ElectromechanicalLOP (const Scheme& scheme_, const IonicModel& ionicmodel_, const GV& gv):
      scheme(scheme_), ionicmodel(ionicmodel_)
      {
        // The idea is to do all computations on the reference element only once here in the constructor.
        // This implies we assume all elements are of the same type, e.g. simplices,
        // and we use the first element of the mesh as a template.
        // Moreover, we want to use compile-time known loop bounds. This is accomplished
        // by checking that the run-time given basius coomplis with the compile-time
        // given numbers.
        
	//Set parameters
	v_min = ionicmodel.ionicpar.v_rest;
        v_max = ionicmodel.ionicpar.v_peak;
        gil   = ionicmodel.par.gil;
        git   = ionicmodel.par.git;
	// Getting the parameter values from the ini file
        mu1 =  ConfigParser::get("ElasticityModel.mu",4.0);
        beta = ConfigParser::get("ElasticityModel.beta",0.3);
        
        // get finite element basis using given rpresentative element
        auto felu = scheme.femu.find(*gv.template begin<0>());
        auto felp = scheme.femp.find(*gv.template begin<0>());
        auto felv = scheme.femv.find(*gv.template begin<0>());
        auto felw = scheme.femw.find(*gv.template begin<0>());

        // check size of the supplied basis
        if (felu.localBasis().size()!=nu) {
          std::cout << "Basis size mismatch for velocity!" << std::endl;
          exit(1);
        }
        if (felp.localBasis().size()!=np) {
          std::cout << "Basis size mismatch for pressure!" << std::endl;
          exit(1);
        }
       if (felv.localBasis().size()!=nv) {
          std::cout << "Basis size mismatch for transmemberain potential!" << std::endl;
          exit(1);
        }
       if (felw.localBasis().size()!=nw) {
          std::cout << "Basis size mismatch for gating variable!" << std::endl;
          exit(1);
        }

        // find quadrature rules with the given number of points and
        // evaluate quadrature rules
        Dune::GeometryType gt = felu.type();
        std::cout << "Electromechanics LOP on " << gt << std::endl;
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
        int orderbnd=-1; // boundary rule
        Dune::GeometryType gtbnd;
        if (gt.isSimplex()) gtbnd = Dune::GeometryTypes::simplex(dim-1);
        if (gt.isCube()) gtbnd = Dune::GeometryTypes::cube(dim-1);
        for (int order=20; order>=1; order--)
          {
            const Dune::QuadratureRule<RF,dim-1>& rule = Dune::QuadratureRules<RF,dim-1>::rule(gtbnd,order);
            if (rule.size()==mbnd)
              {
                orderbnd = order;
                //std::cout << "order of boundary quadrature with " << mbnd << " points is " << order << std::endl;
                for (int i=0; i<mbnd; i++) {
                  qpbnd[i] = rule[i].position();
                  wbnd[i] = rule[i].weight();
                }
                break;
              }
          }
        if (orderbnd<0) {
          std::cout << "Could not find boundary quadruture rule with that many points!" << std::endl;
          exit(1);
        }

        // evaluate basis functions on reference element for velocity, pressure, potential and gating variable
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
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> zeta(nv);
            felw.localBasis().evaluateFunction(qpvol[k],zeta);
            for (int i=0; i<nv; i++) zetahatvol[i][k] = zeta[i];
          }
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> eta(nw);
            felw.localBasis().evaluateFunction(qpvol[k],eta);
            for (int i=0; i<nw; i++) etahatvol[i][k] = eta[i];
          }
        for (const auto& is : intersections(gv,*gv.template begin<0>()))
          {
            auto fgeo_self = is.geometryInInside();
            int face = is.indexInInside();
            //std::cout << "face " << face << std::endl;
            for (int k=0; k<mbnd; k++) // loop over face quadrature points
              {
                auto qp = fgeo_self.global(qpbnd[k]);
                std::vector<RangeType> phi(nu);
                felu.localBasis().evaluateFunction(qp,phi);
                for (int i=0; i<nu; i++) phihatbnd[face][i][k] = phi[i];
              }
          }
        
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
        // extract local function spaces
        const auto& lfsu = lfstrial.template child<0>();    // velocity node
        const auto& lfsp = lfstrial.template child<1>();  // pressure
        const auto& lfsv = lfstrial.template child<2>();  // transmemberain potential
        const auto& lfsw = lfstrial.template child<3>();  // gating variable

        // extract degrees of freedom
        RF ZT[dim][nu]; // velocity degrees of freedom
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j) ZT[i][j] = x(lfsu_i,j);
          }
        RF YT[np]; // pressure degrees of freedom
        for (size_t j=0; j<np; ++j) YT[j] = x(lfsp,j);
        RF QT[nv]; // potential degrees of freedom
        for (size_t j=0; j<nv; ++j) QT[j] = x(lfsv,j);
        RF TT[nw]; // gating variable degrees of freedom
        for (size_t j=0; j<nw; ++j) TT[j] = x(lfsw,j);
        
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

        // compute pressure at all quadratue points
        RF PT[mvol] = {0.0};
        for (size_t j=0; j<np; ++j) // reduction over j
          for (size_t alpha=0; alpha<mvol; ++alpha)
            PT[alpha] += YT[j]*psihatvol[j][alpha];

        // compute transmemberain potential at all quadratue points
        RF VT[mvol] = {0.0};
        for (size_t j=0; j<nv; ++j) // reduction over j
          for (size_t alpha=0; alpha<mvol; ++alpha)
            VT[alpha] += QT[j]*zetahatvol[j][alpha];

        // compute gating variable at all quadratue points
        RF WT[mvol] = {0.0};
        for (size_t j=0; j<nw; ++j) // reduction over j
          for (size_t alpha=0; alpha<mvol; ++alpha)
            WT[alpha] += TT[j]*etahatvol[j][alpha];

        // compute mu1*J0^n(F^{n+1})(F0^n)^-1(F0^n)^-T
        Dune::FieldMatrix<RF,dim,dim> Inter[mvol];
        {
          RF J0 = 0.0;
          Dune::FieldMatrix<RF,dim,dim> F0, F0invT, KronAt, KronAl;
          for (size_t alpha=0; alpha<mvol; ++alpha)
          {
            //auto v = getv(eg.entity(),qpvol[alpha]);
            Dune::FMatrixHelp::multTransposedMatrix(at,KronAt);
            Dune::FMatrixHelp::multTransposedMatrix(al,KronAl);
            RF Gl = -gil*beta*(VT[alpha]-v_min)/(v_max-v_min+VT[alpha]);
            RF Gt = -git*beta*(VT[alpha]-v_min)/(v_max-v_min+VT[alpha]);
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

        // Now come the residual contributions term by term.
        // These are ultimately reductions over the quadrature points but
        // possibly involve some intermediate results beforehand
        RF RUT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]
        RF RPT[np]  = {0.0};       // residual contribution of this element to pressure test function [j]
        RF RVT[np]  = {0.0};       // residual contribution of this element to potential test function [j]
        RF RWT[np]  = {0.0};       // residual contribution of this element to gating variable test function [j]

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
            RPT[j] += (J[alpha]-1.0)*psihatvol[j][alpha]*factor[alpha];

        // Obtain the Diffusion tensor
        typename IonicModel::Traits::TensorType tensor_i;
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(geo.type()).position(0,0);
        ionicmodel.Di(eg.entity(),localcenter,tensor_i);

        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
        {
          // evaluate gradient of shape functions
          auto& gradzetahat = cache.evaluateJacobian(qpvol[alpha],
                             lfsv.finiteElement().localBasis());
          // transform gradients of shape functions to real element
          const auto S = geo.jacobianInverseTransposed(qpvol[alpha]);
          auto gradzeta = makeJacobianContainer(lfsv);
          for (size_t i=0; i<nv; i++)
            S.mv(gradzetahat[i][0],gradzeta[i][0]);
          // compute gradient of v
          Dune::FieldVector<RF,dim> gradv(0.0), intermediate(0.0);
          for (size_t i=0; i<nv; i++)
            gradv.axpy(x(lfsv,i),gradzeta[i][0]);
          // evaluate J0^{n+1}(F^-1)D(F^-T)GradV
          Dune::FieldMatrix<RF,dim,dim> ret1, ret2, Finv;
          Dune::FMatrixHelp::invertMatrix(F[alpha],Finv);
          Dune::FMatrixHelp::multMatrix(Finv,tensor_i,ret1);
          Dune::FMatrixHelp::multMatrix(ret1,FinvT[alpha],ret2);
          intermediate = Dune::FMatrixHelp::mult(ret2,gradv);

          Dune::FieldVector<RF,ncomp> v(0.0);
				  // for (size_t j=0; j<ncomp; j++)
						   v[0] = VT[alpha];
               v[1] = WT[alpha];

          // evaluate ionic/reaction term
					 Dune::FieldVector<RF,ncomp> rhs(0.0);
           auto xglobal = eg.geometry().global(qpvol[alpha]);
					 ionicmodel.rhs(xglobal,v,rhs);
          
          for (size_t i=0; i<nv; i++)
            residual.accumulate(lfsv,i,(intermediate*gradzeta[i][0]
                             -rhs[0]*zetahatvol[i][alpha])*factor[alpha]);
          for (size_t i=0; i<nw; i++)
            residual.accumulate(lfsw,i,-rhs[1]*etahatvol[i][alpha]*factor[alpha]);
        }

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
                // extract local function spaces
        // const auto& lfsu = lfstest.template child<0>();    // velocity node
        // const auto& lfsp = lfstest.template child<1>();  // pressure
        const auto& lfsv0 = lfstest.template child<2>();  // ionic variables potential
        // const auto& lfsv0 = lfsv.child(0);                // Transmemberain potential
        // const auto& lfsv1 = lfsv.child(1); 
        
       for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
        {
          
           auto xglobal = eg.geometry().global(qpvol[alpha]);
           double Istim = ionicmodel.getStimulusSupport(xglobal);
           if(Istim>1e-5)
           {
             RF factor = wvol[alpha] * eg.geometry().integrationElement(qpvol[alpha]);
             for (size_t i=0; i<nv; i++)
             {
               residual.accumulate(lfsv0,i,-Istim*zetahatvol[i][alpha]*factor);
             }
           }
        }

      }
      private:  
	    const IonicModel& ionicmodel;
      const Scheme& scheme;

    };


// class ElectromechanicsMassLOP
//   : public Dune::PDELab::FullVolumePattern,
//     public Dune::PDELab::LocalOperatorDefaultFlags,
//     public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
// {
  
// public:
//   // pattern assembly flags
//   enum { doPatternVolume = true };

//   // residual assembly flags
//   enum { doAlphaVolume = true };

//   // volume integral depending on test and ansatz functions
//   template<typename EG, typename LFSU, typename X,typename LFSV, typename R>
//   void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
//                      const LFSV& lfsv, R& r) const
//   {
//     // select the two components (assume Galerkin scheme U=V)
//     using namespace Dune::Indices;
//     using LFSU_V = Dune::TypeTree::Child<LFSU,2>;
//     using LFSU_W = Dune::TypeTree::Child<LFSU,3>;
//     using RT_V = typename LFSU_V::Traits::FiniteElementType::
//           Traits::LocalBasisType::Traits::RangeType;
//     using RT_W = typename LFSU_V::Traits::FiniteElementType::
//           Traits::LocalBasisType::Traits::RangeType;
//     using RF = typename LFSU_V::Traits::FiniteElementType::
//           Traits::LocalBasisType::Traits::RangeFieldType;

//     // extract local function spaces      
//     auto lfsuv = lfsu.child(2);
//     auto lfsuw = lfsu.child(3);
//     const unsigned int nv = lfsuv.size();
//     const unsigned int nw = lfsuw.size();

//     std::vector<RT_V> zeta(nv);
//     std::vector<RT_W> eta(nw);

//     // select quadrature rule
//     auto geo = eg.geometry();
//     const int order = 2*lfsuv.finiteElement().localBasis().order();
//     auto rule = Dune::PDELab::quadratureRule(geo,order);

//     // loop over quadrature points
//     for (const auto& ip : rule)
//       {
//         // evaluate transmemberain potential basis functions
//         lfsuv.finiteElement().localBasis().evaluateFunction(ip.position(),zeta);
//         // evaluate gating variable basis functions
//         lfsuw.finiteElement().localBasis().evaluateFunction(ip.position(),eta);

//         // evaluate v and w
//         RF v=0.0, w=0.0;
//         for (std::size_t i=0; i<nv; i++) 
//           v += x(lfsuv,i)*zeta[i];
//         for (std::size_t i=0; i<nw; i++) 
//           w += x(lfsuw,i)*eta[i];   

//         // accumulate residuals
//         RF factor=ip.weight()*geo.integrationElement(ip.position());
//         for (std::size_t i=0; i<nv; i++) 
//           r.accumulate(lfsuv,i,v*zeta[i]*factor);
//         for (std::size_t i=0; i<nw; i++) 
//           r.accumulate(lfsuw,i,w*eta[i]*factor);
        
//       }
//   }

//   //! jacobian contribution of volume term
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename M>
//   void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
//                         const LFSV& lfsv, M& mat) const
//   {
//     // select the two components (assume Galerkin scheme U=V)
//     using namespace Dune::Indices;
//     using LFSU_V = Dune::TypeTree::Child<LFSU,2>;
//     using LFSU_W = Dune::TypeTree::Child<LFSU,3>;
//     using RT_V = typename LFSU_V::Traits::FiniteElementType::
//           Traits::LocalBasisType::Traits::RangeType;
//     using RT_W = typename LFSU_V::Traits::FiniteElementType::
//           Traits::LocalBasisType::Traits::RangeType;
//     using RF = typename LFSU_V::Traits::FiniteElementType::
//           Traits::LocalBasisType::Traits::RangeFieldType;

//     // extract local function spaces      
//     auto lfsuv = lfsu.child(2);
//     auto lfsuw = lfsu.child(3);
//     const unsigned int nv = lfsuv.size();
//     const unsigned int nw = lfsuw.size();

//     std::vector<RT_V> zeta(nv);
//     std::vector<RT_W> eta(nw);

//     // select quadrature rule
//     auto geo = eg.geometry();
//     const int order = 2*lfsuv.finiteElement().localBasis().order();
//     auto rule = Dune::PDELab::quadratureRule(geo,order);

//     // loop over quadrature points
//     for (const auto& ip : rule)
//       {
//         // evaluate transmemberain potential basis functions
//         lfsuv.finiteElement().localBasis().evaluateFunction(ip.position(),zeta);
//         // evaluate gating variable basis functions
//         lfsuw.finiteElement().localBasis().evaluateFunction(ip.position(),eta);

//         // accumulate matrix entries
//         RF factor=ip.weight()*geo.integrationElement(ip.position());
//         for (std::size_t j=0; j<nv; j++)
//           for (std::size_t i=0; i<nv; i++) {
//             mat.accumulate(lfsuv,i,lfsuv,j,
//                            zeta[j]*zeta[i]*factor);
//             mat.accumulate(lfsuw,i,lfsuw,j,
//                            eta[j]*eta[i]*factor);
//           }
//       }
//   }

//   //! apply local jacobian of the volume term -> nonlinear variant
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename R>
//   void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
//                               const X& x, const X& z, const LFSV& lfsv,
//                               R& r) const
//   {
//     alpha_volume(eg,lfsu,z,lfsv,r);
//   }

//   //! apply local jacobian of the volume term -> linear variant
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename R>
//   void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
//                               const X& x, const LFSV& lfsv,
//                               R& r) const
//   {
//     alpha_volume(eg,lfsu,x,lfsv,r);
//   }
// };

template<typename FEM, int ncomp=2> // default number of components in ionic model: 2
class ElectromechanicsMassLOP
  : public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  // types
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  using RF = typename LocalBasis::Traits::RangeFieldType;

  // private data members
  Dune::PDELab::LocalBasisCache<LocalBasis> cache; // a cache for local basis evaluations

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // extract local function spaces      
    auto lfsuv0 = lfsu.template child<2>();
    // auto lfsuv0 = lfsuv.child(0);
    auto lfsuv1 = lfsu.template child<3>();
    const unsigned int nv = lfsuv0.size();

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsuv0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions at first child
        auto& phihat = cache.evaluateFunction(ip.position(),
                         lfsuv0.finiteElement().localBasis());

        // evaluate u0 and 
        RF u0=0.0, u1=0.0;
        for (std::size_t i=0; i<nv; i++) {
          u0 += x(lfsuv0,i)*phihat[i];
          u1 += x(lfsuv1,i)*phihat[i];
        }   

        // accumulate residuals
        RF factor=ip.weight()*geo.integrationElement(ip.position());
        for (std::size_t i=0; i<nv; i++) 
        {
          r.accumulate(lfsuv0,i,u0*phihat[i]*factor);
          r.accumulate(lfsuv1,i,u1*phihat[i]*factor);
        }
        
      }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
                        const LFSV& lfsv, M& mat) const
  {
    // extract local function spaces      
    auto lfsuv0 = lfsu.template child<2>();
    // auto lfsuv0 = lfsuv.child(0);
    auto lfsuv1 = lfsu.template child<3>();
    const unsigned int nv = lfsuv0.size();

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsuv0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
         // evaluate basis functions at first child
        auto& phihat = cache.evaluateFunction(ip.position(),
                         lfsuv0.finiteElement().localBasis());

        // accumulate matrix entries
        RF factor=ip.weight()*geo.integrationElement(ip.position());
        for (std::size_t j=0; j<nv; j++)
          for (std::size_t i=0; i<nv; i++) {
            mat.accumulate(lfsuv0,i,lfsuv0,j,
                           phihat[j]*phihat[i]*factor);
            mat.accumulate(lfsuv1,i,lfsuv1,j,
                           phihat[j]*phihat[i]*factor);
          }
      }
  }

  //! apply local jacobian of the volume term -> nonlinear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,x,lfsv,r);
  }
};
