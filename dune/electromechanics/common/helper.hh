#ifndef DUNE_COMMON_HELPER_HH
#define DUNE_COMMON_HELPER_HH


template<typename GV>
typename GV::ctype smallest_edge(const GV &gv) {
    typedef typename GV::ctype DF;
    typedef typename GV::template Codim<0>::
        template Partition<Dune::Interior_Partition>::Iterator Iterator;
    typedef typename GV::template Codim<0>::Geometry Geometry;
    static const std::size_t dim = GV::dimension;
    typedef Dune::FieldVector<DF, GV::dimensionworld> DomainW;
    
    typedef Dune::ReferenceElements<DF, dim> Refelems;
    typedef Dune::ReferenceElement<DF, dim> Refelem;
    
    DF smallest = std::numeric_limits<DF>::infinity();
    
    const Iterator &end = gv.template end<0, Dune::Interior_Partition>();
    for(Iterator it = gv.template begin<0, Dune::Interior_Partition>();
        it != end; ++it)
        {
            auto refelem = Refelems::general(it->type());
            const Geometry &geo = it->geometry();
            for(int edge = 0; edge < refelem.size(dim-1); ++edge) {
                const DomainW &cornerpos0 =
                    geo.corner(refelem.subEntity(edge, dim-1, 0, dim));
                const DomainW &cornerpos1 =
                    geo.corner(refelem.subEntity(edge, dim-1, 1, dim));
                smallest = std::min(smallest, (cornerpos0-cornerpos1).two_norm2());
            }
        }
    
    return std::sqrt(gv.comm().min(smallest));
}

// return only internal nodes which are unique
template<typename GV>
int compute_internalnodes(const GV& gv)
{ 
    const int dim = GV::dimension;
    int gsize = 0; // gv.ghostSize(dim);
    int osize = 0; //gv.overlapSize(dim);
    int isize = 0; // internal node size
    typedef typename GV::template Codim<dim>::Iterator Iterator; 
    Iterator it = gv.template begin<dim>();
    const Iterator &endIt = gv.template end<dim>();
    // check the intersections
    for (; it != endIt; ++it) 
        {
            if (it->partitionType() == Dune::GhostEntity)
                gsize++;
            else if(it->partitionType() == Dune::OverlapEntity)
                osize++;
            else
                isize++;
        }
 
    gsize =gv.comm().sum(gsize);
    osize =gv.comm().sum(osize);
    isize =gv.comm().sum(isize);
  
    // if(gv.comm().rank()==0)
    //     std::cout<<" internal_nodes "<< isize << " ghost nodes "<<gsize<<" overlap nodes "<<osize<<std::endl;
    return isize;

}

// return only internal elements which are unique
template<typename GV>
int compute_internalelements(const GV& gv)
{ 
    const int dim = GV::dimension;
    int gsize = 0; // gv.ghostSize(dim);
    int osize = 0; //gv.overlapSize(dim);
    int isize = 0; // internal node size
    typedef typename GV::template Codim<0>::Iterator Iterator; 
    Iterator it = gv.template begin<0>();
    const Iterator &endIt = gv.template end<0>();
    // check the intersections
    for (; it != endIt; ++it) 
        {
            if (it->partitionType() == Dune::GhostEntity)
                gsize++;
            else if(it->partitionType() == Dune::OverlapEntity)
                osize++;
            else
                isize++;
        }
 
    gsize =gv.comm().sum(gsize);
    osize =gv.comm().sum(osize);
    isize =gv.comm().sum(isize);
  
    // if(gv.comm().rank()==0)
    //     std::cout<<" internal_elements "<< isize << " ghost ele "<<gsize<<" overlap ele "<<osize<<std::endl;
    return isize;

}


/*! set ID for specified point
 */ 
template< typename GV, typename U, int ncomp>
std::array<double,ncomp> getSolutionAtPoint(const GV& gv, const U& u)
{
    typedef double Real;
    typedef typename GV::ctype ct;
    const int dim = GV::dimension;
    int index_id = -1, rank_id = -1;
    Dune::FieldVector<ct,dim> sol_vertex;
    int visited = 0;
    std::array<Real,ncomp> sol_data;
    for(int i=0; i<ncomp; i++)
        sol_data[i] = -1000.0;
  
  
    Dune::FieldVector<ct,dim> given_point(3.0);
    //given_point[0] = 0.422834;  given_point[1] = -0.09347700;   given_point[2] = 1.1084599;
    given_point[0] = 1.65302;  given_point[1] = -0.25779;   given_point[2] = 0.276481;
    //given_point[0] = 1.9236; given_point[1] = 0.108066; given_point[2] = -0.174734;
  
    typedef typename GV::Traits::template Codim<dim>::Iterator VIterator;
    VIterator vendit = gv.template end<dim>();

    using EntitySet = Dune::PDELab::NonOverlappingEntitySet<GV>;
    auto entity_set = EntitySet(gv);      

  	 
    Real small_dist = 0.001;
    //for (VIterator it = gv.template begin<dim>(); it!=vendit; ++it)
    for (const auto& vertex : vertices(entity_set))
        { 
            Dune::FieldVector<ct,dim> globalX =   vertex.geometry().center();//corner(0);
            Real  Dist  = (given_point- globalX).two_norm();
            if(Dist<small_dist)
                {
                    visited = 1;
                    small_dist = std::min(small_dist,Dist);
                    index_id = entity_set.indexSet().index(vertex);//indexSet.index(*it);
                    rank_id =gv.comm().rank();
                    for(int i=0; i<ncomp; i++)
                        sol_data[i] = Dune::PDELab::Backend::native(u)[index_id][i];
                    //gv.comm().broadcast(&sol_data,ncomp,rank_id);
	  
                    //std::cout<<"point "<<globalX[0]<<" "<<globalX[1]<<" "<<globalX[2]<<" "<<rank_id<<" "<<index_id<<std::endl;
                }
        }
	  
    visited = gv.comm().max(visited);
    Real smallest_dist = gv.comm().min(small_dist);
    for(int i=0; i<ncomp; i++)
        sol_data[i] = gv.comm().max(sol_data[i]);
    // make sure that atleast one point surrounding the specified distance
    if(!visited)
        {
            if(gv.comm().rank()==0)
                std::cout<<" WARNING: Could not find the point (to extract solution) which is near to the specified location: smallest distance "<<smallest_dist<<" and visited "<<visited<<std::endl;
            exit(1);
        }

    return sol_data;
}

template< typename GV, typename U, typename V, int ncomp, int points=2>
std::array<double,ncomp*points> getSolutionAtSeveralPoints(const GV& gv,const U& u, const V& v)
{
    typedef double Real;
    typedef typename GV::ctype ct;
    const int dim = GV::dimension;
    int index_id = -1, rank_id = -1;
    Dune::FieldVector<ct,dim> sol_vertex;
    std::array<int,points> visited = {0,0,0};
    std::array<Real,ncomp*points> sol_data;
    for(int i=0; i<ncomp*points; i++)
        sol_data[i] = -1000.0;
  
  
    std::array<Dune::FieldVector<ct,dim>,points> given_point;
  
#if TBunnyC
    given_point[0][0] = 1.4574;  given_point[0][1] = -0.3846;  given_point[0][2] = 0.2531;
    given_point[1][0] = 0.8945;  given_point[1][1] = -0.7356;  given_point[1][2] = 0.2600;
    given_point[2][0] = 0.3554;  given_point[2][1] = -0.8780;  given_point[2][2] = 0.2657;
#else
    given_point[0][0] = 1.0;  given_point[0][1] = 1.0;  given_point[0][2] = 1.0;
    given_point[1][0] = 0.0;  given_point[1][1] = 0.0;  given_point[1][2] = 0.0;
    given_point[2][0] = 0.5;  given_point[2][1] = 0.5;  given_point[2][2] = 0.5;

#endif
    typedef typename GV::Traits::template Codim<dim>::Iterator VIterator;
    VIterator vendit = gv.template end<dim>();

    using EntitySet = Dune::PDELab::NonOverlappingEntitySet<GV>;
    auto entity_set = EntitySet(gv);      

  	 
    std::array<double, points> small_dist={0.0001,0.001,0.001};
  
    for (const auto& vertex : vertices(entity_set))
        { 
            Dune::FieldVector<ct,dim> globalX =   vertex.geometry().center();
            index_id = entity_set.indexSet().index(vertex);;
            for(int p=0; p<points; p++)
                {
                    Real  Dist  = (given_point[p]- globalX).two_norm();
                    if(Dist<small_dist[p])
                        {
                            visited[p] = 1;
                            small_dist[p] = std::min(small_dist[p],Dist);
	      
                            for(int i=0; i<(ncomp-1); i++)
                                sol_data[p*ncomp+i] = Dune::PDELab::Backend::native(v)[index_id][i];
                            sol_data[p*ncomp+(ncomp-1)] = Dune::PDELab::Backend::native(u)[index_id][0];
                        }
                }
        }

    for(int p=0; p<points; p++)
        visited[p] = gv.comm().max(visited[p]);

    Real smallest_dist = gv.comm().min(small_dist[0]);
    for(int i=0; i<ncomp*points; i++)
        sol_data[i] = gv.comm().max(sol_data[i]);
    // make sure that atleast one point surrounding the specified distance
    if(!visited[0] || !visited[1] || !visited[2] )
        {
            if(gv.comm().rank()==0)
                std::cout<<" WARNING: Could not find the point (to extract solution) which is near to the specified location: smallest distance "<<smallest_dist<<" and visited "<<visited [0]<<std::endl;
            exit(1);
        }

    return sol_data;
}

template< typename GV, typename V, int ncomp, int points=2>
std::array<double,ncomp*points> getSolutionAtSeveralPoints(const GV& gv, const V& v)
{
    typedef double Real;
    typedef typename GV::ctype ct;
    const int dim = GV::dimension;
    int index_id = -1, rank_id = -1;
    Dune::FieldVector<ct,dim> sol_vertex;
    std::array<int,points> visited = {0,0,0};
    std::array<Real,ncomp*points> sol_data;
    for(int i=0; i<ncomp*points; i++)
        sol_data[i] = -1000.0;
  
  
    std::array<Dune::FieldVector<ct,dim>,points> given_point;
  
#if TBunnyC
    given_point[0][0] = 1.4574;  given_point[0][1] = -0.3846;  given_point[0][2] = 0.2531;
    given_point[1][0] = 0.8945;  given_point[1][1] = -0.7356;  given_point[1][2] = 0.2600;
    given_point[2][0] = 0.3554;  given_point[2][1] = -0.8780;  given_point[2][2] = 0.2657;
#else
    given_point[0][0] = 1.0;  given_point[0][1] = 1.0;  given_point[0][2] = 1.0;
    given_point[1][0] = 0.0;  given_point[1][1] = 0.0;  given_point[1][2] = 0.0;
    given_point[2][0] = 0.5;  given_point[2][1] = 0.5;  given_point[2][2] = 0.5;

#endif
    typedef typename GV::Traits::template Codim<dim>::Iterator VIterator;
    VIterator vendit = gv.template end<dim>();

    using EntitySet = Dune::PDELab::NonOverlappingEntitySet<GV>;
    auto entity_set = EntitySet(gv);      

  	 
    std::array<double, points> small_dist={0.0001,0.001,0.001};
  
    for (const auto& vertex : vertices(entity_set))
        { 
            Dune::FieldVector<ct,dim> globalX =   vertex.geometry().center();
            index_id = entity_set.indexSet().index(vertex);;
            for(int p=0; p<points; p++)
                {
                    Real  Dist  = (given_point[p]- globalX).two_norm();
                    if(Dist<small_dist[p])
                        {
                            visited[p] = 1;
                            small_dist[p] = std::min(small_dist[p],Dist);
	      
                            for(int i=0; i<(ncomp-1); i++)
                                sol_data[p*ncomp+i] = Dune::PDELab::Backend::native(v)[index_id][i];
	      
                        }
                }
        }

    for(int p=0; p<points; p++)
        visited[p] = gv.comm().max(visited[p]);

    Real smallest_dist = gv.comm().min(small_dist[0]);
    for(int i=0; i<ncomp*points; i++)
        sol_data[i] = gv.comm().max(sol_data[i]);
    // make sure that atleast one point surrounding the specified distance
    if(!visited[0] || !visited[1] || !visited[2] )
        {
            if(gv.comm().rank()==0)
                std::cout<<" WARNING: Could not find the point (to extract solution) which is near to the specified location: smallest distance "<<smallest_dist<<" and visited "<<visited [0]<<std::endl;
            exit(1);
        }

    return sol_data;
}

// boundary grid function selecting boundary conditions 
template<typename GV, int ncomp>
class NeumannBC
   : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
												   BoundaryGridFunctionTraits<
													  GV,Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type,ncomp,
													  Dune::FieldVector<
														 Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type,ncomp> >,
												   NeumannBC<GV,ncomp> >
{
   const GV& gv;

public:
   typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions BC;
   typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,
													Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type,ncomp,
													Dune::FieldVector<Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type,ncomp> > Traits;
   typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,NeumannBC<GV,ncomp> > BaseT;

   NeumannBC (const GV& gv_) : gv(gv_) {
   std::cout << " NeumannBc: Grid size: elements =" << gv.size(0)<<" on rank " << gv.comm().rank()<<" size "<<gv.comm().size()<<std::endl;
   }

   template<typename I>
   inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
						 const typename Traits::DomainType& x,
						 typename Traits::RangeType& y) const
   {  
	  y = BC::Neumann;
   }

  /** Predicate identifying Dirichlet boundaries for velocity. */
      template<typename I>
      bool isDirichlet(const I & intersection,
                       const Dune::FieldVector<typename I::ctype, I::mydimension> & coord) const
      { return false; }

   //! get a reference to the GridView
   inline const GV& getGridView ()
   {
	  return gv;
   }
};

#endif
