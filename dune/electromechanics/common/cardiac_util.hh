#ifndef CARDIAC_UTIL_HH
#define CARDIAC_UTIL_HH

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

namespace Dune
{
    namespace PDELab
    {
	namespace BiDomain
	{
	    template<typename GV, typename RF, int ncomp>
	    struct BiDomainParameterTraits
	    {
		//! \brief the grid view
		typedef GV GridViewType;
		
		//! \brief Enum for domain dimension
		enum { 
		    //! \brief dimension of the domain
		    dimDomain = GV::dimension
		}; 
		   
		//! \brief Export type for domain field
		typedef typename GV::Grid::ctype DomainFieldType;
    
		//! \brief domain type
		typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;
    
		//! \brief domain type
		typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;
    
		//! \brief Export type for range field
		typedef RF RangeFieldType;
    
		//! \brief range type
		typedef Dune::FieldVector<RF,ncomp> RangeType;

		//! \brief Jacobian range type (ncomp x ncomp)
		typedef Dune::FieldMatrix<RF,ncomp,ncomp> JacobianType;
    
		//! \brief permeability tensor type
		typedef Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain> TensorType;
		
		//! grid types
		typedef typename GV::Traits::template Codim<0>::Entity ElementType;
		typedef typename GV::Intersection IntersectionType;
	    };
	    
	    /** \brief Class to define the boundary condition types
	     */
	    struct ElectrodeBoundary
	    {
		enum Type { Anode = 1, None = 0, Cathode = -1}; 
	    };
	    
	}
    } 
}

template <class GV, class RF, int ncomp>
class BidParameters
{ 
    enum {dim=GV::dimension, dimworld = GV::dimensionworld};  
    typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,RF,ncomp> Traits;
public:
  BidParameters(): gil(3.0E-3),git(3.1525E-4),gel(2.0E-3),get(1.3514E-3)
   {     
       //gil= 1E-3; git = 1E-3; gel = 3E-3; get = 3E-3;// works without anisotrophy for LR3v model
       //gil= 1.0E-3; git = 1.0E-3; gel = 1E-3; get = 1E-3;// works without anisotrophy for LR3v model
	  //gil = 1.6E-2; gel=  2.0E-3;git = 2.1E-3; get = 1.0E-3;
       
	  //git = 3.6E-4; get = 1.9E-3;
	  //gil= 1.2E-2; git = 1.24E-3; gel = 8E-3; get = 5.4E-3;
	 
         gil = 3.28e-2; git = 6.99e-3; // MSR model  Am = 50
        double ang = (22.0/7.0 )/4.0 *0.0;
	    
	for (std::size_t i=0; i<dim; i++)
	    for (std::size_t j=0; j<dim; j++)
		if (i==j)
		    { 
			if(i==0)
			    sigma_i[i][j] = gil * cos(ang) * cos(ang) + git * sin(ang) * sin(ang);
			else
			    sigma_i[i][j] = git * cos(ang) * cos(ang) + gil * sin(ang) * sin(ang);
		    }
		else
		    { 
			if(i<j)
			    sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
			else
			    sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
		    }
 
	for (int i=0; i<dim; i++)
	    for (int j=0; j<dim; j++)
		if (i==j)
		    { 
			if(i==0)
			    sigma_e[i][j] = gel * cos(ang) * cos(ang) + get * sin(ang) * sin(ang);
			else
			    sigma_e[i][j] = get * cos(ang) * cos(ang) + gel * sin(ang) * sin(ang);
		    }
		else
		    { 
			if(i<j)
			    sigma_e[i][j] = (gel-get) * cos(ang) * sin(ang);
			else
			    sigma_e[i][j] = (gel-get) * cos(ang) * sin(ang);
		    }
    // sigma_i *= (1.0/8.0);
    // sigma_e *= 0.2;

    }

public: //! \brief permeability tensor type
    typename Traits::TensorType sigma_i, sigma_e;
    double gil,git,gel,get;
};


template <class GV, class RF,int ncomp=1>
class ConductivityTensors_i
{ 

public:
    enum {dim=GV::dimension, dimworld = GV::dimensionworld};  
    typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,RF,ncomp> Traits;

    ConductivityTensors_i(): gil(3.0E-3),git(3.1525E-4),gel(2.0E-3),get(1.3514E-3)  //gil(3.0E-3),git(3.1525E-4)
    {     
	//gil= 4E-4; git = 4E-4; gel = 1E-3; get = 1E-3;// works without anisotrophy
	gil= 1E-3; git = 1E-3; gel = 1E-3; get = 1E-3;// works without anisotrophy
	//gil = 1.6E-2; gel=  2.0E-3;git = 2.1E-3; get = 1.0E-3;
	//gil = 3.0E-3; gel=  3.0E-3; //gil = 8.8E-3;
	//git = 3.6E-4; get = 1.9E-3;
        //gil= 1.0; git = 1.0; gel = 1.0; get = 1.0; // works without anisotrophy   for LR3V model
    }

 //! diffusion tensor
    void Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::TensorType& sigma_i) const
   {
	  sigma_i = 0.0;
	  double ang = 0.0; double factor = 1.0;
      for (std::size_t i=0; i<dim; i++)
          for (std::size_t j=0; j<dim; j++)
              if (i==j)
                  { 
                      if(i==0)
                          sigma_i[i][j] = factor*gil * cos(ang) * cos(ang) + git * sin(ang) * sin(ang);
                      else
                          sigma_i[i][j] = factor*git * cos(ang) * cos(ang) + gil * sin(ang) * sin(ang);
                  }
              else
                  { 
                      if(i<j)
                          sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
                      else
                          sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
                  }
#if PACI20Model
      sigma_i *= 1000;//0.2;
#endif      
   }
    
 //! diffusion tensor
   void Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x,  const typename Traits::DomainType& sigma,
			typename Traits::TensorType& res) const
   {
	  res = 0.0;
	   for (int i=0; i<dim; i++)
		 for (int j=0; j<dim; j++)
			if(i==j)
			   res[i][j] = gil*sigma[i]*sigma[j] + (1.0-sigma[i]*sigma[j])*git;
			else
			   res[i][j] = gil*sigma[i]*sigma[j] + (0.0-sigma[i]*sigma[j])*git;

	   //this is only valid for MS model
#if MSRModel	   
#ifdef YASP2D
           res *= 1.0/8.0*4.0*1.5;
#else
           res *= 1.0/8.0;//(1.0/15.0);
#endif
#endif
   }

   //! diffusion tensor; extracellular for compatibility
   void De (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
			typename Traits::TensorType& sigma) const
   {
	  sigma = 0.0;
   } 

public: //! \brief permeability tensor type
    typename Traits::TensorType sigma_i, sigma_e;
    double gil,git,gel,get;
};


template <class GV, class RF,int ncomp=1>
class ConductivityTensors_e
{ 
public:
    enum {dim=GV::dimension, dimworld = GV::dimensionworld};  
    typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,RF,ncomp> Traits;

    ConductivityTensors_e(): gil(3.0E-3),git(3.1525E-4),gel(2.0E-3),get(1.3514E-3)
    {     
	//gil= 4E-4; git = 4E-4; gel = 1E-3; get = 1E-3;// works without anisotrophy
	gil= 1E-3; git = 1E-3; gel = 1E-3; get = 1E-3;// works without anisotrophy
	//gil = 1.6E-2; gel=  2.0E-3;git = 2.1E-3; get = 1.0E-3;
	//gil = 3.0E-3; gel=  3.0E-3; //gil = 8.8E-3;
	//git = 3.6E-4; get = 1.9E-3;
    // gil= 1.0; git = 1.0; gel = 1.0; get = 1.0;    
        sigma_e = 0.0;
        sigma_e[0][0] = gel;
        sigma_e[1][1] = get;
        if(dim==3)
            sigma_e[2][2] = get;
	   
          //this is only valid for MS model
#if MSRModel	
#ifdef YASP2D
           sigma_e *= 0.4;//0.2;
#else
           sigma_e *= 0.2;
#endif
#endif
    }
    //! diffusion tensor
    void Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
             typename Traits::TensorType& sigma_i) const
    {
        sigma_i = 0.0;
        double ang = 0.0;double factor = 1.0;
        for (std::size_t i=0; i<dim; i++)
            for (std::size_t j=0; j<dim; j++)
                if (i==j)
                    { 
                        if(i==0)
                            sigma_i[i][j] = factor*gil * cos(ang) * cos(ang) + git * sin(ang) * sin(ang);
                        else
                            sigma_i[i][j] = factor*git * cos(ang) * cos(ang) + gil * sin(ang) * sin(ang);
                    }
                else
                    { 
                        if(i<j)
                            sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
                        else
                            sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
		    }
#if PACI20Model
	sigma_i *= 1000;//0.2;
#endif 	
    }

 //! diffusion tensor
   void Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x,  const typename Traits::DomainType& sigma,
	   typename Traits::TensorType& res) const
   { 
	  // typename Traits::DomainType xg = e.geometry().global(x);
	  // if (xg[0]>=0 && xg[0]<=2 )
	  // 	 {
	  // 		if(xg[1]>=0 && xg[1]<=2)
	  // 		   res = sigma_i;
	  // 	 } 
	  // else
	  // 	 res = 0.0;
	  res = 0.0;
	  bool sigma_proceed = false;
	  for (int i=0; i<dim; i++)
        if(std::fabs(sigma[i])>0.0) // if one of them is fine and then proceed
			sigma_proceed = true;
	  if(sigma_proceed)
	  for (int i=0; i<dim; i++)
		 for (int j=0; j<dim; j++)
			if(i==j)
			   res[i][j] = gil*sigma[i]*sigma[j] + (1.0-sigma[i]*sigma[j])*git;
			else
			   res[i][j] = gil*sigma[i]*sigma[j] + (0.0-sigma[i]*sigma[j])*git;
	  
  //this is only valid for MS model
#if MSRModel
#ifdef YASP2D
          res *= 1.0/8.0*4.0*1.5;//(1.0/15.0);
#else
          res *= 1.0/8.0;//(1.0/15.0);
#endif
#endif
   }

   //! diffusion tensor; extracellular for compatibility
   void De (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
			typename Traits::TensorType& sigma) const
   {
	  sigma = sigma_e;
   }

public: //! \brief permeability tensor type
    typename Traits::TensorType sigma_e;
    double gil,git,gel,get;
};

template <class GV, class RF,int ncomp=1>
class ConductivityTensors_i_e
{ 

public:
    enum {dim=GV::dimension, dimworld = GV::dimensionworld};  
    typedef Dune::PDELab::BiDomain::BiDomainParameterTraits<GV,RF,ncomp> Traits;

    ConductivityTensors_i_e(): gil(3.0E-3),git(3.1525E-4)
    {     
	//gil= 4E-4; git = 4E-4; gel = 1E-3; get = 1E-3;// works without anisotrophy
	gil= 1E-3; git = 1E-3; // works without anisotrophy
	//gil = 1.6E-2; gel=  2.0E-3;git = 2.1E-3; get = 1.0E-3;
	//gil = 3.0E-3; gel=  3.0E-3; //gil = 8.8E-3;
	//git = 3.6E-4; get = 1.9E-3;
        //gil= 1.0; git = 1.0; 
    }

        //! diffusion tensor
    void Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::TensorType& sigma_i) const
   {
	  sigma_i = 0.0;
      double ang = 0.0;
      for (std::size_t i=0; i<dim; i++)
          for (std::size_t j=0; j<dim; j++)
              if (i==j)
                  { 
                      if(i==0)
                          sigma_i[i][j] = gil * cos(ang) * cos(ang) + git * sin(ang) * sin(ang);
                      else
                          sigma_i[i][j] = git * cos(ang) * cos(ang) + gil * sin(ang) * sin(ang);
                  }
              else
                  { 
                      if(i<j)
                          sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
                      else
                          sigma_i[i][j] = (gil-git) * cos(ang) * sin(ang);
                  }
#if PACI20Model
      sigma_i *= 1000;//0.2;
#endif 
   }
    
    
 //! diffusion tensor
   void Di (const typename Traits::ElementType& e, const typename Traits::DomainType& x,  const typename Traits::DomainType& sigma,
			typename Traits::TensorType& res) const
   {
	  res = 0.0;
	  bool sigma_proceed = false;
	  for (int i=0; i<dim; i++)
		 if(std::fabs(sigma[i])>0.0)
			sigma_proceed = true;
	  if(sigma_proceed)
		 for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			   if(i==j)
				  res[i][j] = gil*sigma[i]*sigma[j] + (1.0-sigma[i]*sigma[j])*git;
			   else
				  res[i][j] = gil*sigma[i]*sigma[j] + (0.0-sigma[i]*sigma[j])*git;
#if MSRModel
#ifdef YASP2D
          res *= 1.0/8.0*4.0*1.5;//(1.0/15.0);
#else
          res *= 1.0/8.0;//(1.0/15.0);
#endif
#endif

   }

   //! diffusion tensor; extracellular for compatibility
   void De (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
			typename Traits::TensorType& sigma) const
   {
	  sigma = 0.0;
   } 

public: //! \brief permeability tensor type
    typename Traits::TensorType sigma_e;
    double gil,git;
};


template<int dim, class iTupel>
class YaspPartition : public Dune::YLoadBalance<dim>
{
private:
    const iTupel yasppartitions;

public:
    //constructor:
    YaspPartition(const  iTupel& yasppartitions_ )
	: yasppartitions(yasppartitions_ )
    {
    }

    void loadbalance (const iTupel& size, int P, iTupel& dims) const
    {
	dims = yasppartitions;
    }
};

// From Martin's code
/* little extension of std::exception to be used at some point...*/
class BidException : public std::exception
{

public:
    BidException(std::string str) : msg(str)
    {}

    template<typename T>
    BidException(std::string str, const T& info) : msg(str)
    {
	std::stringstream ss;
	ss<<msg<<info;
	msg=ss.str();
    }

    const char* what() const throw()
    {
	if(Dune::MPIHelper::isFake)
	    return msg.c_str();
	std::stringstream ss;
	ss << "Rank #" << Dune::MPIHelper::getCommunication().rank()<<": "<<msg;
	return ss.str().c_str();
    }
    
    virtual ~BidException() throw(){}
    
protected:
    std::string msg;
};



class ConfigParser
{
public:
    static void parse(std::string filename)
    {
	//if(instance == 0) throw BidException("ConfigSet already initialized!");
	Dune::ParameterTreeParser::readINITree(filename,instance);
	verbose = get<bool>("Verbosity.config_parser",false);
    }

    template<class Type>
    static Type get(std::string str)
    {
	// if(instance == 0) throw BidException("ConfigSet is not initialized properly!");
	Type t = instance[str];
	if(verbose && Dune::MPIHelper::getCommunication().rank()==0)
	    std::cout << "Read configuration entry: " << std::setw(15) << std::setfill(' ') << std::left << str << " = " << std::setw(15) << std::setfill(' ') << std::right << t << std::endl;
	return t;
    }

    template<class Type>
    static Type get(std::string str, const Type& _default)
    {
	//if(instance == 0) throw BidException("ConfigSet is not initialized properly!");
	Type t = instance.get(str,_default);
	if(verbose && Dune::MPIHelper::getCommunication().rank()==0)
	    std::cout << "Read configuration entry: " << std::setw(15) << std::setfill(' ') << std::left << str << " = " << std::setw(15) << std::setfill(' ') << std::right << t << std::endl;
	return t;
    
    }
    
protected:
    static Dune::ParameterTree instance;
    static bool verbose;
};

//instantiate static members ; Otherwise liking error
Dune::ParameterTree ConfigParser::instance;
bool ConfigParser::verbose = false;








/*! \brief Configuration Setup
 */
class BidConfigSet
{
public:
   static void parse(std::string filename)
   {
      if(instance != NULL) throw BidException("ConfigSet already initialized!");
      instance = new boost::property_tree::ptree();
      boost::property_tree::ini_parser::read_ini(filename, *instance);
      verbose = get<bool>("Verbosity.config_parser",false);
   }

    static void write(std::string filename)
    {
	if(instance == NULL) throw BidException("ConfigSet not initialized!");
       boost::property_tree::ini_parser::write_ini(filename, *instance);
    }
    
   static boost::property_tree::ptree& getInstance()
   {
      if(instance == 0) throw BidException("ConfigSet is not initialized properly!");
      return *instance;
   }


   template<class Type>
   static Type get(std::string str)
   {
      if(instance == 0) throw BidException("ConfigSet is not initialized properly!");
      Type t = instance->get<Type>(str);
      if(verbose && Dune::MPIHelper::getCommunication().rank() ==0)
         std::cout << "Read configuration entry: " << std::setw(15) << std::setfill(' ') << std::left << str << " = " << std::setw(15) << std::setfill(' ') << std::right << t << std::endl;
      return t;
   }
   template<class Type>
   static Type get(std::string str, const Type& _default)
   {
      try
      {
         return BidConfigSet::get<Type>(str);
      }
      catch(...)
      {
         if(verbose)
            std::cout << "Can't read configuration entry \"" << str << "\" using default value '" << _default << "'" <<std::endl;
         return _default;
      }
   }

   template<class Type>
   static void set(std::string key, const Type& val)
   {
      if(instance == 0) instance = new boost::property_tree::ptree();
      instance->add(key,val);
   }

   static boost::property_tree::ptree& get_child(std::string str)
   {
      if(instance == 0) throw BidException("ConfigSet is not initialized properly!");
      return instance->get_child(str);
   }

protected:
   static boost::property_tree::ptree *instance;
   static bool verbose;
};

/*! instantiate static members
 */
bool BidConfigSet::verbose = true;
boost::property_tree::ptree* BidConfigSet::instance = 0;




#endif
