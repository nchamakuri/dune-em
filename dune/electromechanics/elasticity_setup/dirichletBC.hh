#ifndef DIRICHLET_BC_HH
#define DIRICHLET_BC_HH

namespace Dune{
    namespace Heart{
	namespace Elasticity{

        //! Define Scalar Dirichlet Boundary Conditions for displacement and pressure
        template<typename GV, typename MODEL, typename RF>
            class Scalar_BC :
                public Dune::PDELab::AnalyticGridFunctionBase<
                Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                Scalar_BC<GV,MODEL,RF> >,
                public Dune::PDELab::InstationaryFunctionDefaults
        {
            public:
		enum{dim = GV::dimension};
                typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
                typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Scalar_BC<GV,MODEL,RF>> BaseT;

                typedef typename Traits::DomainType DomainType;
                typedef typename Traits::RangeType RangeType;

                // Constructor (if no i_dim specifyed then it is pressure BC)
                Scalar_BC(const GV & gv, MODEL& model_, int i_dim_=-1) :
                    BaseT(gv), model(model_), i_dim(i_dim_){ }

                template<typename I>
                    bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::mydimension> & x) const
                    {
                        Dune::FieldVector<double,dim> xg = ig.geometry().global( x );
                        return model.isDirichlet(xg, i_dim);
                    }

                inline void evaluateGlobal(const DomainType & x, RangeType & u) const
                {
                    u = model.evaluateDirichlet(x,i_dim);
                } // end inline function evaluateGlobal

                void setDof(int degree_of_freedom){
                    dof = degree_of_freedom;
                }

            private:
                int dof;
                MODEL& model;
                int i_dim;
        };

       }
    }
}

#endif
