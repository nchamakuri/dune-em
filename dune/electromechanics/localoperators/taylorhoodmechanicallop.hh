// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODMECHANICAL_HH
#define DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODMECHANICAL_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

namespace Dune {
  namespace PDELab {
    namespace Elasticity {
    

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** \brief A local operator for the nonlinear elasticity model used in electro-mechanical coupling models.

        This class implements a local operator for conforming finite element 
        discretizations of the active strain formulation of Neo-Hookean materials with
        TaylorHood basis.

        Find $\boldsymbol{d}, p$ in suitable admissible displacement and pressure spaces such that

        $$
        \begin{gathered}
        \int_{\Omega_{o}}\left(\mu_{1} J_{o}(\mathbf{I}+\nabla \boldsymbol{d}) \mathbf{F}_{o}^{-1} \mathbf{F}_{o}^{-T}: 
          \nabla \boldsymbol{\varphi}-p (\mathbf{I}+\nabla \boldsymbol{d})^{-T}: \nabla \boldsymbol{\varphi}\right)=0 \\
        \int_{\Omega_{o}}(J-1) q=0
        \end{gathered}
        $$

        for all test functions $\varphi, q$. Concerning boundary conditions, specifying a displacement field 
        on a portion $\Gamma_{D}$ of $\partial \Omega_{o}$, which for simplicity has been taken as homogeneous
         Dirichlet boundary data

        $$
        \boldsymbol{d}=\mathbf{0} \quad \text { on } \Gamma_{D} \subset \partial \Omega_{o}
        $$

        and homogeneous Neumann conditions on $\partial \Omega_{o} \backslash \Gamma_{D}$.

    */

        template<typename Params, typename GetV = void, typename P0GFS=void, typename TD=void>
    class FiniteElasticityLOP :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags
    {
    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      //enum { doLambdaVolume = true };
      //enum { doLambdaBoundary = true };

        typedef typename Dune::PDELab::LocalFunctionSpace<P0GFS, Dune::PDELab::TrialSpaceTag> P0LFSU;
        			
      // Constructor for Electrical-Mechanical coupled problem
        FiniteElasticityLOP  (const Params& params_, const GetV& getv_,  const P0GFS& p0gfs_, const TD& sigma_, int superintegration_order_ = 0):params(params_),
                                                                                                                                                 getv(getv_), active_strain(true), p0gfs(p0gfs_), sigma(sigma_), superintegration_order(superintegration_order_)
      {}

      // Constructor for Mechanical problem
      FiniteElasticityLOP (const Params& params_, int superintegration_order_ = 0): params(params_), 
                        active_strain(false), superintegration_order(superintegration_order_)
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // define types
        using namespace Indices;
        using LFSU_V_PFS = TypeTree::Child<LFSU,_0>;
        using LFSU_V = TypeTree::Child<LFSU_V_PFS,_0>;
        using LFSU_P = TypeTree::Child<LFSU,_1>;
        using DF = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType;
        using RF = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RT_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using JacobianType_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using RT_P = typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // extract local function spaces
        const auto& lfsu_v_pfs = child(lfsu,_0);
        const unsigned int nu = lfsu_v_pfs.child(0).size();
        const auto& lfsu_p = child(lfsu,_1);
        const unsigned int np = lfsu_p.size();

        // dimensions
        const int dim = EG::Geometry::mydimension;

	// get geometry
        auto geo = eg.geometry();

        // determine quadrature order
        const int v_order = lfsu_v_pfs.child(0).finiteElement().localBasis().order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        typename EG::Geometry::JacobianInverseTransposed jac;
        std::vector<Dune::FieldVector<RF,dim> > gradphi(nu);
        std::vector<RT_P> psi(np);
        Dune::FieldVector<RF,dim> vu(0.0);
        std::vector<RT_V> phi(nu);
        Dune::FieldMatrix<RF,dim,dim> jacu(0.0), F(0.0), FinvT(0.0);
        RF J;
        Dune::FieldMatrix<RF,dim,dim> F0, F0invT, KronAt, KronAl;
        Dune::FieldMatrix<RF,dim,dim> inter(0.0);

        // Extract tensor values
        typedef Dune::PDELab::LFSIndexCache<P0LFSU> LFSU0Cache;
        typedef typename TD::template ConstLocalView<LFSU0Cache> U0View;
        P0LFSU p0lfsu(p0gfs);
        LFSU0Cache lfsu0Cache(p0lfsu);
        U0View xview(sigma);
        p0lfsu.bind(eg.entity());
        // local coefficients
        LocalVector<typename TD::ElementType, TrialSpaceTag> sigma_data(p0lfsu.maxSize());
        lfsu0Cache.update();
        xview.bind(lfsu0Cache);
        xview.read(sigma_data);
        xview.unbind();
        // store it in a FieldVector
        Dune::FieldVector<DF,dim> y(0.0); 
        for (size_t i=0; i<dim; i++)
				  y[i] = sigma_data(p0lfsu.child(i),0);
        
        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> js(nu);
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(ip.position(),js);

            // transform gradient to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_t i=0; i<nu; i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // evaluate pressure basis functions
            lfsu_p.finiteElement().localBasis().evaluateFunction(ip.position(),psi);

            // evaluate displacement basis functions
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateFunction(ip.position(),phi);
            
            // compute displacement
            for(int d=0; d<dim; ++d)
              {
                vu[d] = 0.0;
                const auto& lfsu_v = lfsu_v_pfs.child(d);
                for (size_t i=0; i<lfsu_v.size(); i++)
                  vu[d] += x(lfsu_v,i) * phi[i];
              }

            // Compute displacement jacobian
            for(int d=0; d<dim; ++d){
              jacu[d] = 0.0;
              const auto& lfsu_v = lfsu_v_pfs.child(d);
              for (size_t i=0; i<lfsu_v.size(); i++)
                jacu[d].axpy(x(lfsu_v,i),gradphi[i]);
            }

            // compute pressure
            RT_P func_p(0.0);
            for (size_t i=0; i<lfsu_p.size(); i++)
              func_p += x(lfsu_p,i) * psi[i];

            // compute deformation gradient
            for (size_t i=0; i<dim; ++i)
              for (size_t j=0; j<dim; ++j)
                F[i][j] = 1.0*(i==j)+jacu[i][j];

            // compute the F^(-T) and J = det(F)
            J = Dune::FMatrixHelp::invertMatrix_retTransposed(F,FinvT);

            // compute J0*F0^-1*F0^-T for active strain model
            Dune::FieldMatrix<RF,dim,dim> J0F0invF0invT(0.0);
            for (size_t i=0; i<dim; ++i) J0F0invF0invT[i][i] = 1.0;
            if(active_strain) // Incase of electromechanical problem
            {  
                  auto v = getv(eg.entity(),ip.position());
		  //const auto at = params.at(eg.entity(),ip.position());
		  //const auto al = params.al(eg.entity(),ip.position());
                  Dune::FMatrixHelp::multTransposedMatrix(params.at,KronAt);
                  Dune::FMatrixHelp::multTransposedMatrix(params.al,KronAl);

                  RF Gl = -params.gil*params.beta*(v-params.v_min)/(params.v_max-params.v_min+v);
                  RF Gt = -params.git*params.beta*(v-params.v_min)/(params.v_max-params.v_min+v);

                  for (size_t i=0; i<dim; ++i)
                    for (size_t j=0; j<dim; ++j)
                        F0[i][j] = 1.0*(i==j)+Gl*KronAl[i][j]+Gt*KronAt[i][j];

                  RF J0 = Dune::FMatrixHelp::invertMatrix_retTransposed(F0,F0invT);

                  Dune::FMatrixHelp::multTransposedMatrix(F0invT,J0F0invF0invT);
                  J0F0invF0invT *= J0;
            }
             inter = F*J0F0invF0invT;

            // Viscosity and density
            // const auto mu = params.mu(eg.entity(),ip.position());
            // const auto rho = params.rho(eg.entity(),ip.position());

            // geometric weight
            const auto factor = ip.weight() * geo.integrationElement(ip.position());

            for(int d=0; d<dim; ++d)
            {
              const auto& lfsu_v = lfsu_v_pfs.child(d);
              for (size_t i=0; i<nu; i++)
              {
                // integrate mu * J0*F*F0^-1*F0^-T : e_d*grad phi_i
                r.accumulate(lfsu_v,i, params.mu * (inter[d] * gradphi[i]) * factor);

                // integrate -p*FinvT : e_d*grad phi_i
                r.accumulate(lfsu_v,i,- func_p * (FinvT[d]* gradphi[i]) * factor);
              }
            }
            // integrate (J-1) * psi_i
            for (size_t i=0; i<np; i++)
                r.accumulate(lfsu_p,i, -(J - 1.0) * psi[i] * factor);

          }
      }


      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // define types
        using namespace Indices;
        using LFSU_V_PFS = TypeTree::Child<LFSU,_0>;
        using LFSU_V = TypeTree::Child<LFSU_V_PFS,_0>;
        using LFSU_P = TypeTree::Child<LFSU,_1>;
        using RF = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using RT_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;
        using JacobianType_V = typename LFSU_V::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType;
        using RT_P = typename LFSU_P::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType;

        // extract local function spaces
        const auto& lfsu_v_pfs = child(lfsu,_0);
        const unsigned int vsize = lfsu_v_pfs.child(0).size();
        const auto& lfsu_p = child(lfsu,_1);
        const unsigned int psize = lfsu_p.size();

        // dimensions
        const int dim = EG::Geometry::mydimension;
        
	// get geometry
        auto geo = eg.geometry();

        // determine quadrature order
        const int v_order = lfsu_v_pfs.child(0).finiteElement().localBasis().order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*v_order - 1 + jac_order + det_jac_order + superintegration_order;

        // Initialize vectors outside for loop
        typename EG::Geometry::JacobianInverseTransposed jac;
        std::vector<JacobianType_V> js(vsize);
        std::vector<Dune::FieldVector<RF,dim> > gradphi(vsize);
        std::vector<RT_P> psi(psize);
        std::vector<RT_V> phi(vsize);
        Dune::FieldVector<RF,dim> vu(0.0);
        
        Dune::FieldMatrix<RF,dim,dim> jacu(0.0), F(0.0), Finv(0.0);
        Dune::FieldMatrix<RF,dim,dim> F0, F0invT, KronAt, KronAl;
        Dune::FieldMatrix<RF,dim,dim> inter(0.0);

        // loop over quadrature points
        for (const auto& ip : quadratureRule(geo,qorder))
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            lfsu_v_pfs.child(0).finiteElement().localBasis().evaluateJacobian(ip.position(),js);

            // transform gradient to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_t i=0; i<vsize; i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // evaluate basis functions
            lfsu_p.finiteElement().localBasis().evaluateFunction(ip.position(),psi);

            // Compute displacement jacobian
            for(int d=0; d<dim; ++d){
              jacu[d] = 0.0;
              const auto& lfsu_v = lfsu_v_pfs.child(d);
              for (size_t i=0; i<lfsu_v.size(); i++)
                jacu[d].axpy(x(lfsu_v,i),gradphi[i]);
            }

            // compute pressure
            RT_P func_p(0.0);
            for (size_t i=0; i<lfsu_p.size(); i++)
              func_p += x(lfsu_p,i) * psi[i];

            // compute deformation gradient
            for (size_t i=0; i<dim; ++i)
              for (size_t j=0; j<dim; ++j)
                F[i][j] = 1.0*(i==j)+jacu[i][j];

            // compute the F^(-T) and J = det(F)
            RF J = Dune::FMatrixHelp::invertMatrix(F,Finv);

            // Viscosity and density # 
            //const auto mu = params.mu(eg.entity(),ip.position());
            //const auto rho = params.rho(eg.entity(),ip.position());

            // compute J0*F0^-1*F0^-T for active strain model
            Dune::FieldMatrix<RF,dim,dim> J0F0invF0invT(0.0);
            for (size_t i=0; i<dim; ++i) J0F0invF0invT[i][i] = 1.0;
            if(active_strain) // incase of electromechanical problem 
            {  
                  auto v = getv(eg.entity(),ip.position());
		  //const auto at = params.at(eg.entity(),ip.position());
		  //const auto al = params.al(eg.entity(),ip.position());
                  Dune::FMatrixHelp::multTransposedMatrix(params.at,KronAt);
                  Dune::FMatrixHelp::multTransposedMatrix(params.al,KronAl);

                  RF Gl = -params.gil*params.beta*(v-params.v_min)/(params.v_max-params.v_min+v);
                  RF Gt = -params.git*params.beta*(v-params.v_min)/(params.v_max-params.v_min+v);

                  for (size_t i=0; i<dim; ++i)
                    for (size_t j=0; j<dim; ++j)
                        F0[i][j] = 1.0*(i==j)+Gl*KronAl[i][j]+Gt*KronAt[i][j];

                  RF J0 = Dune::FMatrixHelp::invertMatrix_retTransposed(F0,F0invT);
                  Dune::FMatrixHelp::multTransposedMatrix(F0invT,J0F0invF0invT);
                  J0F0invF0invT *= J0;
            }
             inter = J0F0invF0invT;

            const auto factor = ip.weight() * geo.integrationElement(ip.position());

            for(int d=0; d<dim; ++d)
            {

              const auto& lfsv_v = lfsu_v_pfs.child(d);
              const auto& lfsu_v = lfsv_v;

              for (size_t i=0; i<lfsv_v.size(); i++){

                // integrate mu1* J0 * Grad phi_u_j e_d * F0^-1 *F0^-T : grad phi_v_i e_d
                for (size_t j=0; j<lfsv_v.size(); j++){
                  for (int jd=0; jd<dim; ++jd)
                     for (int id=0; id<dim; ++id )
                        mat.accumulate(lfsv_v,i,lfsu_v,j, params.mu * inter[jd][id]*(gradphi[i][id] * gradphi[j][jd]) * factor);
                }
                  // integrate p*(F^-1*phi_u_j e_dd )^T : F^-1*phi_v_i e_d
                for (size_t j=0; j<lfsv_v.size(); j++){
                  for(int dd=0; dd<dim; ++dd){
                    const auto& lfsu_v = lfsu_v_pfs.child(dd);
                    for (int jd=0; jd<dim; ++jd)
                      for (int id=0; id<dim; ++id )
                        mat.accumulate(lfsv_v,i,lfsu_v,j, func_p * Finv[id][dd]*gradphi[j][jd] * Finv[jd][d]*gradphi[i][id] * factor);
                    }
                }

                // integrate -psi_j * F^-T : phi_v_i e_d
                for (size_t j=0; j<psize; j++)
                  for (int id=0; id<dim; ++id )
                    mat.accumulate(lfsv_v,i,lfsu_p,j, - (Finv[id][d] * gradphi[i][d] * psi[j]) * factor);
              }

              // integrate J*F^-T : phi_u_j e_d* psi_i
              for (size_t i=0; i<psize; i++){
                for (size_t j=0; j<lfsu_v.size(); j++)
                  for (int jd=0; jd<dim; ++jd)
                    mat.accumulate(lfsu_p,i,lfsu_v,j,  -1.0 * (Finv[jd][d] * gradphi[j][jd] * psi[i]) * factor);
              }
            } // d
          } // ip
      }

      //! apply local jacobian of the volume term -> nonlinear variant
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv, R& r) const
      {
          alpha_volume(eg,lfsu,z,lfsv,r);
      }

      // dummy jacobian_boundary for local operator to assemble only in the overlap with other subdomains
      template<typename IG, typename LFSU, typename X, typename LFSV, typename LocalMatrix>
      void jacobian_boundary(const IG& ig, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                                   LocalMatrix& mat) const
      {
      }


      const bool active_strain;
    private:
      const Params& params;
      const int superintegration_order;
      const GetV& getv;
        const TD& sigma;
        const P0GFS& p0gfs;
    };

    } // namespace Elasticity
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODMECHANICAL_HH
