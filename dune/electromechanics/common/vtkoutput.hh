#ifndef DUNE_USER_VTKOUT_HH
#define DUNE_USER_VTKOUT_HH

#include<dune/pdelab/gridfunctionspace/vtk.hh>

// VTK output using recursie functions
template<typename GV, typename BidHelper>
class VTKOutput : public Output
{
    //enum {dim = };
    
    typedef double real;
    typedef Dune::VTKSequenceWriter<GV> VTKWriter;
  
public:
  
  VTKOutput(const GV& gv_, const BidHelper& bidhelper_, real intervall = -1, //ConfigParser::get("VTKOutput.vtk_interval",-1),
	    bool force = false //ConfigParser::get("VTKOutput.vtk_force_interval",false) ):
	    ): gv(gv_), bidhelper(bidhelper_), vtkwriter(gv_,"bidomain",this->path+ "/vtk","",Dune::VTK::conforming),vtkwriter_vm(gv_,"bidomain_vm",this->path+ "/vtk","",Dune::VTK::conforming)
    {
      //std::cout<<" VTKOut: path in constructor "<<this->path<<std::endl;
    }

protected:
    // write J.th component, use dummy bool argument. so only partial specialization is required, which is allowed by the standard...
    template<bool dummy, unsigned J>
    struct Helper
    {
	template<typename V>
	static void write_vtk(VTKWriter& writer, const BidHelper& bidhelper, const V& v, real time)
        {
          // extract GFS from PGFS
          typedef Dune::PDELab::GridFunctionSubSpace<typename BidHelper::VGFS,Dune::TypeTree::StaticTreePath<J>> GFS;
          GFS gfs(bidhelper.vgfs);
          typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
          
          DGF dgf(gfs, v);

          // add data to writer
          writer.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,bidhelper.ionicmodel.component_name(J)));

          // recurse to next component
          Helper<true,J-1>::write_vtk(writer, bidhelper,v, time);

          // here the dgf and gfs go out of scope and will be deleted, hence its lifetime exceeds the writer.write() and writer.clear() calls
        }
    };


    // partial specialization for 0.th component, dont recurse!
    template<bool dummy>
    struct Helper<dummy,0>
    {
        template<typename V>
        static void write_vtk(VTKWriter& writer, const BidHelper& bidhelper, const V& v,real time)
        {
          // extract GFS from PGFS
          typedef Dune::PDELab::GridFunctionSubSpace<typename BidHelper::VGFS,Dune::TypeTree::StaticTreePath<0>> GFS;
	  //typename BidHelper::VGFS vgfs = bidhelper.vgfs;
           GFS gfs(bidhelper.vgfs);

          // generate DiscreteGridFunction
          typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
          DGF dgf(gfs, v);

          // add data to writer
          writer.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,bidhelper.ionicmodel.component_name(0)));
          writer.write(time,Dune::VTK::appendedraw);
		  writer.clear();
        }
    };

public:
    template<typename V>
    void write(const V& v, real t0)
   {
     // extract the solution here since otherwise it will be extractet for each component recursion wich might become expensive

      // recursive writing of all components. this can be used for different PGFS components
       Helper<true,BidHelper::ncomp-1>::write_vtk(vtkwriter, bidhelper, v, t0);
   }

  template<typename U, typename V>
  void write(const U& u, const V& v, real t0)
  {
    // extract the solution here since otherwise it will be extractet for each component recursion wich might become expensive
    
    // recursive writing of all components. this can be used for different PGFS components
    typedef Dune::PDELab::DiscreteGridFunction<typename BidHelper::UGFS,U> DGF;
    DGF dgfu(bidhelper.ugfs,u);
   
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfu,"u"));
    Helper<true,BidHelper::ncomp-1>::write_vtk(vtkwriter, bidhelper, v, t0);
    
  }

  template<typename V>
    void write_vm(const V& v,real t0)
    {
    	typedef Dune::PDELab::DiscreteGridFunction<typename BidHelper::UGFS,V> DGF;
    	DGF dgfu(bidhelper.ugfs,v);
   
    	vtkwriter_vm.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfu,"v"));
    	vtkwriter_vm.write(t0,Dune::VTK::appendedraw); 
    	vtkwriter_vm.clear();
    }

   virtual void output(real t0)
   {
   }
    virtual void output_dns(real t0)
   {
   }

protected:

    const GV& gv;
    const BidHelper& bidhelper;
    // make them static so that they can be deleted when signal received
    Dune::VTKSequenceWriter<GV> vtkwriter, vtkwriter_vm;
};


template<typename GV, typename BidHelper, typename FiniteElasticity>
class VTKOutput_mechano : public Output
{
    //enum {dim = };
    
    typedef double real;
    typedef Dune::VTKSequenceWriter<GV> VTKWriter;
    using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
    using Path1 = Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>;
public:
  
    VTKOutput_mechano(const GV& gv_, const BidHelper& bidhelper_, const FiniteElasticity& fiel_,real intervall = ConfigParser::get("VTKOutput.vtk_interval",-1),
	      bool force = ConfigParser::get("VTKOutput.vtk_force_interval",false) ):
      gv(gv_), bidhelper(bidhelper_), fiel(fiel_), vtkwriter(gv_,"electrophysiology",this->path+ "/vtk","",Dune::VTK::conforming),
      vtkwriter_vm(gv_,"bidomain_vm",this->path+ "/vtk","",Dune::VTK::conforming)
    {
      //std::cout<<" VTKout: path in constructor "<<this->path<<std::endl;
    }

protected:
    // write J.th component, use dummy bool argument. so only partial specialization is required, which is allowed by the standard...
    template<bool dummy, unsigned J>
    struct Helper
    {
	template<typename V>
	static void write_vtk(VTKWriter& writer, const BidHelper& bidhelper, const V& v, real time)
        {
          // extract GFS from PGFS
          typedef Dune::PDELab::GridFunctionSubSpace<typename BidHelper::VGFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<J>>> GFS;
          GFS gfs(bidhelper.vgfs);
          typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
          
          DGF dgf(gfs, v);

          // add data to writer
          writer.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,bidhelper.ionicmodel.component_name(J)));

          // recurse to next component
          Helper<true,J-1>::write_vtk(writer, bidhelper,v, time);

          // here the dgf and gfs go out of scope and will be deleted, hence its lifetime exceeds the writer.write() and writer.clear() calls
        }
    };


    // partial specialization for 0.th component, dont recurse!
    template<bool dummy>
    struct Helper<dummy,0>
    {
        template<typename V>
        static void write_vtk(VTKWriter& writer, const BidHelper& bidhelper, const V& v,real time)
        {
          // extract GFS from PGFS
          typedef Dune::PDELab::GridFunctionSubSpace<typename BidHelper::VGFS,Path0> GFS;
	  //typename BidHelper::VGFS vgfs = bidhelper.vgfs;
           GFS gfs(bidhelper.vgfs);

          // generate DiscreteGridFunction
          typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
          DGF dgf(gfs, v);

          // add data to writer
          writer.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,bidhelper.ionicmodel.component_name(0)));
          writer.write(time,Dune::VTK::appendedraw);
		  writer.clear();
        }
    };

public:
    template<typename V>
    void write(const V& v, real t0)
   {
     // extract the solution here since otherwise it will be extractet for each component recursion wich might become expensive

      // recursive writing of all components. this can be used for different PGFS components
       Helper<true,BidHelper::ncomp-1>::write_vtk(vtkwriter, bidhelper, v, t0);
   }

  template<typename U, typename V>
  void write(const U& u, const V& v, real t0)
  {
    // extract the solution here since otherwise it will be extractet for each component recursion wich might become expensive
    
    // recursive writing of all components. this can be used for different PGFS components
    typedef Dune::PDELab::DiscreteGridFunction<typename BidHelper::UGFS,U> DGF;
    DGF dgfu(bidhelper.ugfs,u);
   
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfu,"u"));
    Helper<true,BidHelper::ncomp-1>::write_vtk(vtkwriter, bidhelper, v, t0);
    
  }

    template<typename UP, typename V>
  void write_all(const UP& up, const V& v, real t0)
  {
    typedef  Dune::PDELab::GridFunctionSubSpace<typename FiniteElasticity::GFSUP,Path0> DisplacementSubGFS;
    DisplacementSubGFS dispsubgfs(fiel.gfsup);
    typedef Dune::PDELab::VectorDiscreteGridFunction<DisplacementSubGFS,UP> DGFSU;
    DGFSU dgfsu(dispsubgfs,up);
    typedef Dune::PDELab::GridFunctionSubSpace<typename FiniteElasticity:: GFSUP,Path1> PressureSubGFS;
    PressureSubGFS psubgfs(fiel.gfsup);
    typedef Dune::PDELab::DiscreteGridFunction<PressureSubGFS,UP> PressureDGFS;
    PressureDGFS pdgfs(psubgfs,up);

    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFSU> >(dgfsu,"displacement"));
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PressureDGFS> >(pdgfs,"pressure"));
    Helper<true,BidHelper::ncomp-1>::write_vtk(vtkwriter, bidhelper, v, t0);
    
  }

  template<typename V>
    void write_vm(const V& v,real t0)
    {
    	typedef Dune::PDELab::DiscreteGridFunction<typename BidHelper::UGFS,V> DGF;
    	DGF dgfu(bidhelper.ugfs,v);
   
    	vtkwriter_vm.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfu,"v"));
    	vtkwriter_vm.write(t0,Dune::VTK::appendedraw); 
    	vtkwriter_vm.clear();
    }

   virtual void output(real t0)
   {
   }
    virtual void output_dns(real t0)
   {
   }

protected:

    const GV& gv;
    const BidHelper& bidhelper;
    const FiniteElasticity& fiel;
    // make them static so that they can be deleted when signal received
    Dune::VTKSequenceWriter<GV> vtkwriter, vtkwriter_vm;
};


template<typename GV, typename CancerModel,typename VGFS>
class VTKOutput_cancer : public Output
{
    typedef double real;
    typedef Dune::VTKSequenceWriter<GV> VTKWriter;

public:
  VTKOutput_cancer(const GV& gv_, const CancerModel& cancermodel_, real intervall = -1,  bool force = false): 
        gv(gv_), cancermodel(cancermodel_),  vtkwriter(gv_,"cancer",this->path+ "/vtk","",Dune::VTK::conforming)
    {
      std::cout<<" VTKOut: path in constructor "<<this->path<<std::endl;
    }

    template<typename V>
    void write(const V& v, const VGFS& vgfs, real t0)
   {
       Dune::PDELab::addSolutionToVTKWriter(vtkwriter,vgfs,v);
       vtkwriter.write(t0,Dune::VTK::appendedraw);
       vtkwriter.clear();
   }
   virtual void output(real t0)
   {
   }
    virtual void output_dns(real t0)
   {
   }
protected:
    const GV& gv;
    const CancerModel& cancermodel;
    Dune::VTKSequenceWriter<GV> vtkwriter;
};


// VTK output using recursie functions
template<typename Domain>
class VTKOutput_intra : public Output
{
    //enum {dim = };
    
    typedef double real;
    typedef Dune::VTKSequenceWriter<typename Domain::GV> VTKWriter;
public:
  
    VTKOutput_intra(const Domain& domain_, 
	      real intervall = ConfigParser::get("VTKOutput.vtk_interval",-1),
         bool force = ConfigParser::get("VTKOutput.vtk_force_interval",false) ):
	   domain(domain_),
	   vtkwriter(domain.gv,domain.name(),this->path+ "/vtk","",Dune::VTK::conforming),
	   vtkwriter_dns(domain.gv,domain.name_dns(),this->path+ "/vtk","",Dune::VTK::conforming)
    {
    }

protected:
    // write J.th component, use dummy bool argument. so only partial specialization is required, which is allowed by the standard...
    template<bool dummy, unsigned J>
    struct Helper
    {
	template<typename V>
	static void write_vtk(VTKWriter& writer, const Domain& domain, const V& v, real time)
        {
          // extract GFS from PGFS
          typedef Dune::PDELab::GridFunctionSubSpace<typename Domain::VGFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<J>>> GFS;
          GFS gfs(domain.vgfs);
          typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
          
          DGF dgf(gfs, v);

          // add data to writer
          writer.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,domain.ionicmodel.component_name(J)));

          // recurse to next component
          Helper<true,J-1>::write_vtk(writer, domain, v, time);

          // here the dgf and gfs go out of scope and will be deleted, hence its lifetime exceeds the writer.write() and writer.clear() calls
        }
    };


    // partial specialization for 0.th component, dont recurse!
    template<bool dummy>
    struct Helper<dummy,0>
    {
        template<typename V>
        static void write_vtk(VTKWriter& writer, const Domain& domain, const V& v,real time)
        {
          // extract GFS from PGFS
          typedef Dune::PDELab::GridFunctionSubSpace<typename Domain::VGFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>> GFS;
          GFS gfs(domain.vgfs);

          // generate DiscreteGridFunction
          typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
          DGF dgf(gfs, v);

          // add data to writer
          writer.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,domain.ionicmodel.component_name(0)));
          writer.write(time,Dune::VTK::appendedraw);
          writer.clear();
        }
    };

public:
   virtual void output(real t0)
   {
     // extract the solution here since otherwise it will be extractet for each component recursion wich might become expensive
//      typename CalciumDynamics::ERDomain::V erv = dynamics.erDomain().solution(t0);
//      typename CalciumDynamics::CYDomain::V cyv = dynamics.cyDomain().solution(t0);
      typename Domain::V v = domain.v;

      // recursive writing of all components. this can be used for different PGFS components
	  Helper<true,Domain::ncomp-1>::write_vtk(vtkwriter, domain, v, t0);
   }

   virtual void output_dns(real t0)
   {
	  // extract the solution here since otherwise it will be extractet for each component recursion wich might become expensive
	  //      typename CalciumDynamics::ERDomain::V erv = dynamics.erDomain().solution(t0);
	  //      typename CalciumDynamics::CYDomain::V cyv = dynamics.cyDomain().solution(t0);
      typename Domain::V v = domain.v;

      // recursive writing of all components. this can be used for different PGFS components
      Helper<true,Domain::ncomp-1>::write_vtk(vtkwriter_dns, domain, v, t0);
   }
protected:

//    const CalciumDynamics& dynamics;
   const Domain& domain;
   // make them static so that they can be deleted when signal received
   Dune::VTKSequenceWriter<typename Domain::GV> vtkwriter;
   Dune::VTKSequenceWriter<typename Domain::GV> vtkwriter_dns;
};

// VTK output using recursie functions
template<typename Domain>
class VTKOutput_extra : public Output
{
    //enum {dim = };
    
    typedef double real;
    typedef Dune::VTKSequenceWriter<typename Domain::GV> VTKWriter;
public:
  
    VTKOutput_extra(const Domain& domain_, 
	      real intervall = ConfigParser::get("VTKOutput.vtk_interval",-1),
         bool force = ConfigParser::get("VTKOutput.vtk_force_interval",false) ):
	domain(domain_),
	vtkwriter(domain.gv,domain.name(),this->path+ "/vtk","",Dune::VTK::conforming),
	vtkwriter_dns(domain.gv,domain.name_dns(),this->path+ "/vtk","",Dune::VTK::conforming)
    {
    }

public:
   virtual void output(real time)
   { 
	  typename Domain::U u = domain.u;
	  typedef Dune::PDELab::DiscreteGridFunction<typename Domain::UGFS,typename Domain::U> DGF;
	  DGF dgfu(domain.ugfs,u);
   
	  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfu,"u"));
	  vtkwriter.write(time,Dune::VTK::appendedraw); 
	  vtkwriter.clear();
   }
 
   virtual void output_dns(real time)
   { 
	  typename Domain::U u = domain.u;
	  typedef Dune::PDELab::DiscreteGridFunction<typename Domain::UGFS,typename Domain::U> DGF;
	  DGF dgfu(domain.ugfs,u);
   
	  vtkwriter_dns.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfu,"u"));
	  vtkwriter_dns.write(time,Dune::VTK::appendedraw); 
	  vtkwriter_dns.clear();
   }
protected:

//    const CalciumDynamics& dynamics;
   const Domain& domain;
   // make them static so that they can be deleted when signal received
   Dune::VTKSequenceWriter<typename Domain::GV> vtkwriter;
   Dune::VTKSequenceWriter<typename Domain::GV> vtkwriter_dns;
};


// VTK output writer for 1 component
// VTK output writer
template<class GV, class GFS, class V>
void vtkout_cell (const GV& gv, const GFS& gfs, const V& v, const char* name)
{
   typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
   DGF dgfv(gfs,v);
   Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
   vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfv,"vm_sol"));
   vtkwriter.write(name,Dune::VTK::appendedraw); 
}

// VTK output writer
template<class GV, class GFS, class V>
void vtkout (const GV& gv, const GFS& gfs, const V& v, const char* name)
{
   typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
   DGF dgfv(gfs,v);
   
   Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
   vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfv,"v"));
   vtkwriter.write(name,Dune::VTK::appendedraw); 
}

template<class GV, class GFS, class V>
void vtkout_1comp(const GV& gv, const GFS& gfs, const V& v, const char* name)
{
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>> sol1_SUB;
    sol1_SUB sol1_sub(gfs);
    typedef Dune::PDELab::DiscreteGridFunction<sol1_SUB,V> sol1_DGF;
    sol1_DGF sol1_dgf(sol1_sub,v);
    
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol1_DGF> >(sol1_dgf,"v"));
    vtkwriter.write(name,Dune::VTK::appendedraw);
    vtkwriter.clear();
}

template<class GV, class GFS, class V>
void vtkout_2comp(const GV& gv, const GFS& gfs, const V& v, const char* name)
{
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>> sol1_SUB;
	typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>> sol2_SUB;
    sol1_SUB sol1_sub(gfs);
	sol2_SUB sol2_sub(gfs);
    typedef Dune::PDELab::DiscreteGridFunction<sol1_SUB,V> sol1_DGF;
	typedef Dune::PDELab::DiscreteGridFunction<sol2_SUB,V> sol2_DGF;
    sol1_DGF sol1_dgf(sol1_sub,v);
	sol2_DGF sol2_dgf(sol2_sub,v);
    
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol1_DGF> >(sol1_dgf,"v"));
	vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol2_DGF> >(sol2_dgf,"w"));
    vtkwriter.write(name,Dune::VTK::appendedraw);
    vtkwriter.clear();
}



template<class GV, class GFS, class V>
void vtkout_tensors (const GV& gv, const GFS& gfs, const V& v, const char* name)
{
   typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> DGF;
   DGF dgfv(gfs,v);
   //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1);
   Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
   vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgfv,"tensors"));
   vtkwriter.write(name,Dune::VTK::appendedraw); 
}

        
// VTK output writer
template<class GV, class GFS, class V>
void vtkout (const GV& gv, const GFS& gfs, const V& u, const char* name,const char* pathname)
{
#ifdef Q1 
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
#else
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
#endif
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>> sol1_SUB;
    sol1_SUB sol1_sub(gfs);
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>> sol2_SUB;
    sol2_SUB sol2_sub(gfs); 
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<2>>> sol3_SUB;
    sol3_SUB sol3_sub(gfs); 
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<3>>>  sol4_SUB;
    sol4_SUB sol4_sub(gfs);
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<4>>>  sol5_SUB;
    sol5_SUB sol5_sub(gfs); 
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<5>>> sol6_SUB;
    sol6_SUB sol6_sub(gfs); 
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::HybridTreePath<Dune::index_constant<6>>> sol7_SUB;
    sol7_SUB sol7_sub(gfs);

    typedef Dune::PDELab::DiscreteGridFunction<sol1_SUB,V> sol1_DGF;
    sol1_DGF sol1_dgf(sol1_sub,u);
    typedef Dune::PDELab::DiscreteGridFunction<sol2_SUB,V> sol2_DGF;
    sol2_DGF sol2_dgf(sol2_sub,u); 
    typedef Dune::PDELab::DiscreteGridFunction<sol3_SUB,V> sol3_DGF;
    sol3_DGF sol3_dgf(sol3_sub,u);
    typedef Dune::PDELab::DiscreteGridFunction<sol4_SUB,V> sol4_DGF;
    sol4_DGF sol4_dgf(sol4_sub,u);
    typedef Dune::PDELab::DiscreteGridFunction<sol5_SUB,V> sol5_DGF;
    sol5_DGF sol5_dgf(sol5_sub,u); 
    typedef Dune::PDELab::DiscreteGridFunction<sol6_SUB,V> sol6_DGF;
    sol6_DGF sol6_dgf(sol6_sub,u);
    typedef Dune::PDELab::DiscreteGridFunction<sol7_SUB,V> sol7_DGF;
    sol7_DGF sol7_dgf(sol7_sub,u);

    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol1_DGF>>(sol1_dgf,"sol1"));
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol2_DGF>>(sol2_dgf,"sol2"));
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol3_DGF>>(sol3_dgf,"sol3")); 
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol4_DGF>>(sol4_dgf,"sol4"));
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol5_DGF>>(sol5_dgf,"sol5"));
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol6_DGF>>(sol6_dgf,"sol6")); 
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<sol7_DGF>>(sol7_dgf,"sol7"));
    vtkwriter.pwrite(name,pathname,"",Dune::VTK::appendedraw); 
}

#endif

