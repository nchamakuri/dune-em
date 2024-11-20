#include <boost/filesystem.hpp>
#include <boost/bind.hpp> 


class OutputBase
{
protected:
   std::string path;

   OutputBase()
   {
      path = getPath();
      
      // write the parameter set
      // BidConfigSet::write(path + "/parameters.txt");
      //std::cout<<" Created Path is "<<path<<std::endl;
   }

public:
   /* NOTE: in the first place, this method has to be called on all processes since it requires global communication! */
   static std::string getPath()
   {
      static std::string path = "";										// initialize to empty

      if(path == "") 														// path not yet set
      {
         path = ConfigParser::get<std::string>("Output.path");

         std::string temp = path;
         unsigned suffix  = 0;

         while(boost::filesystem::exists(path))							// loop until i found a new suffix
         {
            std::stringstream ss;
            ss << temp << ++suffix; 									// add suffix to path
            path = ss.str();
         }

         Dune::MPIHelper::getCommunication().barrier();		// barrier until all processes have found the right prefix

         if(Dune::MPIHelper::getCommunication().rank()==0) 	// create folder on root processor
         {
            boost::filesystem::create_directory(path);
			boost::filesystem::create_directory(path + "/vtk");
            boost::filesystem::create_directory(path + "/xdmf");
            BidConfigSet::write(path + "/parameters.txt");
            
            if(temp != path)
               std::cout << "Path \"" << temp << "\" already exists, creating " << path << "!" << std::endl;
            else
               std::cout << "Path \"" << temp << "\" does not exist. it will now be created!" << std::endl;
         }


         while(!boost::filesystem::exists(path))
         {
            //idle until every rank can see the created directory
         }
	 
         Dune::MPIHelper::getCommunication().barrier();		// barrier all processes until output path has been created
	 
	 
      }
      
      return path;
   }

   static void writeCPUTime(double t)
   {
      if(Dune::MPIHelper::getCommunication().rank()!=0) return;

      std::ofstream extra((OutputBase::getPath() + "/cputime").c_str(), std::ios::out);
      extra<<t;
   }
};

/*
 * intermediate template layer, just store the constant references of
 * interesting objects and provide timing variables
 */
class Output : public OutputBase
{
protected:

   typedef double real;

   real out_intervall;
   real next_time_step;
   bool flg_force_time_stepping;

    //template<typename Dynamics>
   Output(real _out_intervall = -1., bool _force_time_stepping = false):
      out_intervall(_out_intervall),
      next_time_step(_out_intervall),
      flg_force_time_stepping(_force_time_stepping)
   {
      //TODO: slots are bound to object lifetime, at object destructions slots should be disconnected
      //      this means if the outputs lifetime is shorter than the lifetime of the Dynamics object it can lead to errors
      //dyn.connectTimeStepListener(boost::bind(&Output::time_step_slot, this, _1));
      //dyn.connectTauRestriction(boost::bind(&Output::propose_tau_slot, this,_1));
     //std::cout<<" Output: path in constructor "<<this->path<<std::endl;
   }

public:
   real propose_tau_slot(real t0) const
   {
      if(flg_force_time_stepping)
      {
         // return remaining time until next output time
         if (next_time_step == t0)
            return next_time_step + out_intervall- t0;
         else
            return next_time_step - t0;
      }
      else
         return 100000;// writer does not care about time step
   }

   void time_step_slot(real t0)
   {
      if(t0 >= next_time_step || t0 == 0.)
      {
         next_time_step = t0 + out_intervall;
         output(t0);
      }
   }

   virtual void output(real t0) = 0; // abstract base method for type dependent output
   virtual void output_dns(real t0) = 0; // abstract base method for type dependent output
      
   virtual ~Output()
   {
      //TODO: implement slot disconnections
   }
};

