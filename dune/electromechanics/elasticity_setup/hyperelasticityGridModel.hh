#ifndef CREATE_STRUCTUREMODEL__HH
#define CREATE_STRUCTUREMODEL__HH

#include <dune/composites/Setup/Cijkl.hh>
#include <dune/composites/Setup/Region.hh>

#include <dune/composites/Driver/Solvers/solver_info.hh>
#include <dune/composites/Driver/Solvers/solvers.hh>

using namespace std;

namespace Dune{
    namespace Composites{

        //!  structuredGridModel class.
        /*!
          structuredGridModel provides a base class for 3D structure grid models built on
          using Dune::YaspGrid and Dune::GeometryGrid. The class provides all default functions for definition of a user defined model.
          */
        template <typename YGRID>
        class baseGridModel{

            public:
                //Types
                typedef typename YGRID::LeafGridView YGV;
                typedef Dune::FieldVector<double,3> FieldVec;
                typedef Dune::FieldVector<double,6> Tensor2;
                typedef Dune::FieldMatrix<double,6,6> Tensor4;

                std::string UG_input;
                //! Size of overlap between subdomains for domain decomposition
                int overlap;
		double E, nu, mu, lambda;

                std::array<std::vector<double>,3>  coords;

                //! Array which stores periodicity of grid in each of the 3 dimensions, default value are set to all false.
                std::bitset<3> periodic;
                Dune::MPIHelper& helper;
                std::vector<Tensor4> baseTensor;

                //! Number of Different types of materials defined for applications
                int numMaterials;
                std::vector<Material> myMaterials; //! std::vector containing information of each material defined
                std::vector<Region> R; //! std::vector containing information of each Region

                std::vector<int> elemIndx2PG; //! std::vector containing information of each material defined
                std::vector<double> rot123; //! Vector defining rotation of ply in material coordinates.
                std::vector<Tensor4> C; //! Vector storing elastic tensor within each element

                bool verbosity = false; //! Variable controls level of output to interface
                double tolerance; //! Tolerance of iterative solver
                int solver_type; //! Variable which controls solver type (UMFPack for sequential and GENEO for parallel)
                int solver_tolerance; //! Tolerance of Solver (ignored for a direct solver)
                int Krylov_Verb; //! Variable Controlling Verbosity of Krylov Solver
                int Krylov_Max_It; //! Variable define max iterations
                int solver_verb = 0;

                bool residualHeat; //! Set whether to use residual heat stresses, default false

                std::vector<double> stackSeq;
                std::vector<std::vector<double>> sig;
                bool stress_Plot; //! Variable if stress is plotted out in vtk format
                bool stress_Plot_ascii;
                bool sigResized;
                std::string vtk_displacement_output;
                std::string vtk_stress_output;
                std::string vtk_properties_output;
                int intOrder; //! Gaussian order of integration (default value is 5)
                std::vector<int> nel; //! Number of elements in each direction
                bool isParallel; //! Boolean variable which indicates if run is parallel or not, this later will allow to run multiple solves in sequential mode even though helper.size() > 1
                bool setUp_Required;

                SolverInfo solverParameters; //! Instance of SolverInfo to define solver parameters and store solver results e.g. iterations, converged etc.

                double q; //pressure value

                bool failureCriterionDefined;

                int refineBaseGrid;

                double gravity = 9.80665;
                //! The main constructor
                /*!
                  Constructor of structuredGridModel, requiring only the MPI Helper defined by dune
                  */
                baseGridModel(Dune::MPIHelper& helper_, bool GenEO = false) : helper(helper_){

                    isParallel = false;
                    if (helper.size() > 1){
                        isParallel = true;
                    }

                    intOrder = config.get<int>("fe.intOrder",5);

                    if(helper.rank()==0)solver_verb = 1;
                    sigResized = false;
                    stress_Plot = true;
                    stress_Plot_ascii = false;

                    vtk_displacement_output = "baseModel_displacement";
                    vtk_stress_output = "baseModel_stress";
                    vtk_properties_output = "baseModel_properties";

                    residualHeat = false;
                    setUp_Required = true;
                    failureCriterionDefined = false;
                    refineBaseGrid = 0;
                };

                void inline setup(){
                    // setup -
                    LayerCake();
                    setUpMaterials();
                    computeElasticTensors();
                }


                void inline sizeStressVector(int nel){
                    sig.resize(nel);
                    for (unsigned int i = 0; i < sig.size(); i++){
                        sig[i].resize(6);
                    }
                }

                void inline setStress(int id, Dune::FieldVector<double,6>& element_stress){
                    for (int i = 0; i < 6; i++){
                        sig[id][i] = element_stress[i];
                    }
                }

                Dune::FieldVector<double,6> inline getStress(int id){
                    Dune::FieldVector<double,6> element_stress;
                    for (int i = 0; i < 6; i++){
                        element_stress[i] = sig[id][i];
                    }
                    return element_stress;
                }

                int inline getIntOrder() const {
                    return 5;
                }

                //! returns true when thermal strain is applied
                void inline setThermal(bool t) {
                    residualHeat = t;
                }

                int inline whichRegion(const FieldVec& x){
                    double Lx = 0.0;
                    int ans = 0;
                    bool flag = false;
                    for (unsigned int k = 0; k < R.size(); k++){
                        Lx += R[k].L[0];
                        if (x[0] < Lx && flag == false){ ans = k; flag = true; }
                    }
                    return ans;
                }

                int inline whichLayer(FieldVec& x, unsigned int r){
                    assert(r < R.size());
                    bool flag = false;
                    double z = 0.0;
                    int ans = 0;
                    for (unsigned int k = 0; k < R[r].layers.size(); k++){
                        z += R[r].layers[k].getThickness();
                        if (x[2] < z && flag == false){ ans = k; flag = true; }
                    }
                    return ans;
                }

                //! Function which evaluates the thermal stress tensor sig\_thermal on an element.
                void inline evaluateHeat(Tensor2& f,int id) const{
                    f = 0.0;
                    if(residualHeat != false){
                        double deltaT = -160.; //temperature difference
                        double alpha11 = -0.342 * 1e-6;
                        double alpha22 =  25.8  * 1e-6;
                        double alpha33 =  25.8  * 1e-6;
                        double alpharesin = 25.8 * 1e-6;
                        if(getMaterialTypeFromElement(id) == 0){
                            f[0] = alpharesin*deltaT;
                            f[1] = alpharesin*deltaT;
                            f[2] = alpharesin*deltaT;
                        }
                        else{
                            f[0] = alpha11*deltaT;
                            f[1] = alpha22*deltaT;
                            f[2] = alpha33*deltaT;
                        }
                    }
                }

                inline void evaluateWeight(Dune::FieldVector<double,3>& f, int id) const{
                    int pg = elemIndx2PG[id];
                    f[2] = - gravity * myMaterials[pg].density;
                }

                inline void setPressure(double q_new) {
                    q = q_new;
                }

                //! function which evaluates the Neumann boundary conditions at the point x
                inline void evaluateNeumann(const Dune::FieldVector<double,3> &x,Dune::FieldVector<double,3>& h,const  Dune::FieldVector<double,3>& normal) const{
                    h = 0.0;
                }

                //! \fn A member function which sets the solver type
                void inline setSolver(){
                    solver_type = 0; // UmfPack as default solver
                    if (helper.size() > 1){ // Parallel Solver
                        solver_type = 1;
                    }
                }

                int inline getElements(int i){
                    return nel[i];
                }

                /*! \brief A member taking which constructs the base layered laminate in untransformed coordinates
                 *
                 * LayerCake is a function which defines a flat three-dimensional tensor product grid and the grid's properties
                 * (e.g. periodicity); building  $\hat \Omega := [0,L] \times [0,W] \times [0,T]$,
                 * where $L$, $W$ and $T$ are the length, width and thickness of the laminate
                 * This function provides default values which can be overwritten by defining your own class
                 *
                 */
                void inline LayerCake(){

                    // Defines Default values for Base Class - can be overwritten by defining your own derived class

                    // Define Periodic Boundary Conditions
                    periodic[0] = false; periodic[1] = false; periodic[2] = false;

                    // Uniform Tensor Product Grid on [0,1]^3
                    std::vector<int> nel(3);
                    nel[0] = 50; nel[1] = 50; nel[2] = 50;
                    for (int i = 0; i < 3; i++){
                        double h = 1. / nel[i];
                        coords[i].resize(nel[i] + 1);
                        coords[i][0] = 0.0;
                        for (int j = 1; j <= nel[i]; j++){
                            coords[i][j] = coords[i][j-1] + h;
                        }
                    }
                }

                /*!
                 * \brief Function to read LayerCake from File
                 *
                 * The function provides an interface to read all the information needed to set up a LayerCake, i.e.
                 * size of tensor product grid, periodicity information, etc. from a .csv file.
                 * A user can also define the grid, layering and region properties manually.
                 *
                 */
                void inline LayerCakeFromFile(string pathToCsv){

                    ifstream csvFile(pathToCsv);
                    if(! csvFile) throw std::runtime_error("Could not open file " + pathToCsv);

                    string item;
                    int numRegions = 0;
                    int maxPly = 0;

                    string line;
                    getline(csvFile, line);
                    istringstream lineStream(line);

                    for (int j = 0; j < 2; j++){
                        // Obtain number of regions
                        getline(lineStream, item, ',');
                        if(helper.rank() == 0) std::cout << item << std::endl;
                    }

                    numRegions = atof(item.c_str());
                    if(helper.rank() == 0) std::cout << "Number of Regions = " << numRegions << std::endl;
                    R.resize(numRegions);

                    getline(csvFile, line);
                    istringstream lS2(line);

                    // Obtain max number of plies
                    getline(lS2, item, ',');
                    getline(lS2, item, ',');
                    maxPly = atof(item.c_str());
                    if(helper.rank() == 0) std::cout << "maxPly = " << maxPly << std::endl;

                    getline(csvFile, line);
                    istringstream lS3(line);
                    for (int j = 0; j < 4; j++){
                        // Obtain max number of plies
                        getline(lS3, item, ',');
                        if (j > 0){
                            periodic[j-1] = atof(item.c_str());
                        }
                    }
                    if(helper.rank() == 0) std::cout << "periodic = " << periodic << std::endl;

                    for(int j = 0; j < numRegions; j++){ // For each regions
                        //if(helper.rank() == 0) std::cout << "========= Region = " << j << std::endl;

                        R[j].setNumLayers(maxPly);
                        getline(csvFile, line); getline(csvFile, line);
                        int numParameters = 1;
                        istringstream lS4(line);
                        for (int ii = 0; ii < numParameters; ii++){
                            getline(lS4, item, ',');
                            if (ii == 0){
                                R[j].setType(atof(item.c_str()));
                                numParameters = R[j].numParameters;
                            }
                            else{
                                R[j].setParameter(ii-1,atof(item.c_str()));
                            }
                        } // For each parameter

                        for (int k = 0; k < maxPly; k++){
                            //if(helper.rank() == 0) std::cout << "Ply Number = " << k << std::endl;
                            getline(csvFile, line); // Obtain next row
                            istringstream lS5(line);
                            getline(lS5,item,',');
                            R[j].layers[k].setElements(atof(item.c_str()));

                            //if(helper.rank() == 0) std::cout << "Elements " << atof(item.c_str()) << std::endl;
                            getline(lS5,item,',');

                            //if(helper.rank() == 0) std::cout << "Material " << atof(item.c_str()) << std::endl;
                            R[j].layers[k].setMaterial(atof(item.c_str()) );
                            getline(lS5,item,',');
                            R[j].layers[k].setThickness(atof(item.c_str()));

                            //if(helper.rank() == 0) std::cout << "Thickness " << atof(item.c_str()) << std::endl;
                            getline(lS5,item,',');
                            R[j].layers[k].setOrientation(atof(item.c_str()));

                            //if(helper.rank() == 0) std::cout << "Orientation " << atof(item.c_str()) << std::endl;

                        } // for each ply
                    } // For each region
                }


                /*!
                 * \brief Function defining mesh grading
                 *
                 * This is a helper function that defines a graded mesh in the format used by ABAQUS
                 *
                 */
                std::vector<double> meshRatio(int N, double scaling, int direction){
                    if(direction == 0){ // ie no scaling
                        std::vector<double> h(N);
                        for(int i = 0; i < N; i++)
                            h[i] = 1./N;
                        return h;
                    }
                    else{
                        //meshing with refinement in one direction
                        if(direction==1){
                            double k = scaling;
                            if(scaling<1) k = 1.0/scaling;
                            std::vector<double> h(N);
                            double h0 = 1./N*1./(1.+(k-1)/2.);
                            for(int i = 0; i < N; i++){
                                h[i] = (i*(k-1)/(N-1)+1.0)*h0;
                            }
                            if(scaling<1){
                                std::vector<double> h2(N);
                                for(int i = 0; i< N; i++){
                                    h2[i]   = h[N-1-i];
                                }
                                return h2;
                            }
                            return h;


                        }
                        //meshing with refinement in both directions
                        else{
                            if((N/2)*2 != N)
                                N += 1; //ensure this number is even
                            if(N == 2){
                                return {0.5, 0.5};
                            }
                            int M = N/2;
                            std::vector<double> h(N), a(M);
                            int k2 = scaling;
                            if(scaling<1) k2 = 1.0/scaling;
                            for(int i = 0; i < M; i++){
                                a[i] = float(i)*(k2-1)/(M-1)+1.;
                            }
                            double h0 = 1./M*1./(1.+(k2-1)/2.);
                            for(int i = 0; i < M; i++){
                                h[i] = a[i]*h0/2.;
                                h[N-1 - i] = a[i]*h0/2.;
                            }
                            if(scaling<1){
                                auto h2 = h;
                                for(int i = 0; i< M; i++){
                                    h2[i]   = h[M-1-i];
                                    h2[M+i] = h[N-1-i];
                                }
                                return h2;
                            }
                            return h;
                        }
                    }
                }


                /* ! \brief Build geometry with mesh grading
                 *
                 *  The grading can be given using directions: 0,1,2 for no grading, grading in one direction or grading in both directions
                 *  and grading giving the strength of the grading, the grading needs to be given for each region of the mesh
                 *
                 */
                void inline GeometryBuilder(std::vector<Dune::FieldVector<int,3>> directions, std::vector<Dune::FieldVector<double,3>> grading){
                    // Builds Overall Geometry and Mesh from a set of regions
                    if(helper.rank() == 0) std::cout << "=== Building Geometry" << std::endl;

                    Dune::FieldVector<double,3> origin(0.0); // Define Origin, the default of this will be set to [0.0,0.0,0.0]

                    // Compute number of elements in x and y. Since structures are current prismatic and tensor product grid y elements must be the same so we check this
                    int numRegions = R.size();

                    int nelx = 0;
                    int nely = R[0].nel[1];

                    for (int i = 0; i < numRegions;i++){
                        nelx += R[i].nel[0];
                        assert(R[i].nel[1] == nely); //"Overall tensor product grid y must be equal across all regions, check geometry input file");
                    }

                    if(helper.rank() == 0) std::cout << "Elements in x direction " << std::endl;

                    coords[0].resize(nelx + 1);
                    coords[0][0] = origin[0];

                    int k = 1; // initialise counter
                    for (int i = 0; i < numRegions; i++){
                        auto hx = meshRatio(R[i].nel[0], grading[i][0], directions[i][0]); // This makes grid uniform across a region
                        R[i].nel[0] = hx.size(); //resize if mesh size not even
                        for (int j = 0; j < R[i].nel[0]; j++){
                            coords[0][k] = coords[0][k-1] + R[i].L[0] * hx[j];
                            k++; // Increment Counter in the x-direction
                        }
                    }

                    // === Prismatic mesh in y
                    coords[1].resize(nely + 1);
                    coords[1][0] = origin[0];
                    auto hy = meshRatio(R[0].nel[1], grading[0][1],directions[0][1]);
                    for (int j = 1; j < R[0].nel[1] + 1; j++){
                        coords[1][j] = coords[1][j-1] + R[0].L[1] * hy[j-1];
                    }

                    // Build Coordinates in z
                    // For now let's assume uniform thickness (not necessary though)
                    int nelz = 0;
                    for (unsigned int i = 0; i < R[0].layers.size(); i++){
                        // for each layer
                        nelz += R[0].layers[i].getElements();
                    }
                    coords[2].resize(nelz + 1);
                    coords[2][0] = origin[2];

                    int e = 0;
                    for (unsigned int i = 0; i < R[0].layers.size(); i++){
                        auto hz = meshRatio(R[0].layers[i].getElements(), grading[0][2], directions[0][2]);
                        for (int k = 0; k < R[0].layers[i].getElements(); k++){
                            e++;
                            coords[2][e] = coords[2][e-1] + R[0].layers[i].getThickness() * hz[k];
                        }
                    }
                    R[0].L[2] = coords[2][e];
                }


                void inline GeometryBuilder(){
                    // Builds Overall Geometry and Mesh from a set of regions
                    if(helper.rank() == 0) std::cout << "=== Building Geometry" << std::endl;

                    Dune::FieldVector<double,3> origin(0.0); // Define Origin, the default of this will be set to [0.0,0.0,0.0]

                    // Compute number of elements in x and y. Since structures are current prismatic and tensor product grid y elements must be the same so we check this
                    int numRegions = R.size();

                    int nelx = 0;
                    int nely = R[0].nel[1];

                    for (int i = 0; i < numRegions;i++){
                        nelx += R[i].nel[0];
                        assert(R[i].nel[1] == nely); //"Overall tensor product grid y must be equal across all regions, check geometry input file");
                    }

                    if(helper.rank() == 0) std::cout << "Elements in x direction " << std::endl;

                    coords[0].resize(nelx + 1);
                    coords[0][0] = origin[0];

                    int k = 1; // initialise counter

                    for (int i = 0; i < numRegions; i++){
                        double hx = R[i].L[0] / R[i].nel[0]; // This makes grid uniform across a region TO DO:: Update for non-uniform grids
                        for (int j = 0; j < R[i].nel[0]; j++){
                            coords[0][k] = coords[0][k-1] + hx;
                            k++; // Increment Counter in the x-direction
                        }
                    }

                    // === Prismatic mesh in y
                    coords[1].resize(nely + 1);
                    coords[1][0] = origin[0];
                    double hy = R[0].L[1] / R[0].nel[1];
                    for (int j = 1; j < R[0].nel[1] + 1; j++){
                        coords[1][j] = coords[1][j-1] + hy;
                    }

                    // Build Coordinates in z
                    // For now let's assume uniform thickness (not necessary though)
                    int nelz = 0;
                    for (unsigned int i = 0; i < R[0].layers.size(); i++){
                        // for each layer
                        nelz += R[0].layers[i].getElements();
                    }
                    coords[2].resize(nelz + 1);
                    coords[2][0] = origin[2];

                    int e = 0;
                    for (unsigned int i = 0; i < R[0].layers.size(); i++){
                        double hz = R[0].layers[i].getThickness() / R[0].layers[i].getElements();
                        for (int k = 0; k < R[0].layers[i].getElements(); k++){
                            e++;
                            coords[2][e] = coords[2][e-1] + hz;
                        }
                    }
                    R[0].L[2] = coords[2][e];
                }


                int inline getSizeCoords(int i){
                    return coords[i].size();
                }


                void setUpMaterials(){
                    numMaterials = 3;
                    myMaterials.resize(numMaterials);

                    // Material I: Isotropic Resin
                    myMaterials[0].setType(0); // Isotropic
                    myMaterials[0].setProp(0,10000.); // - E,  Young's Modulus
                    myMaterials[0].setProp(1,0.35); // - nu,   Poisson Ratio
                    myMaterials[0].setDensity(1.57e-5); // - rho,  Density  kg/mm^2

                    // Material II: Orthotropic, Composite Ply - AS4/8552
                    myMaterials[1].setType(1); // Orthotropic Composite
                    myMaterials[1].setProp(0,162000.); // E11
                    myMaterials[1].setProp(1,10000.); // E22
                    myMaterials[1].setProp(2,10000.); // E33
                    myMaterials[1].setProp(3,0.35); // nu12
                    myMaterials[1].setProp(4,0.35); // nu13
                    myMaterials[1].setProp(5,0.49); // nu23
                    myMaterials[1].setProp(6,5200.); // G12
                    myMaterials[1].setProp(7,5200.); // G13
                    myMaterials[1].setProp(8,3500.); // G23
                    myMaterials[1].setDensity(1.57e-5); // - rho, Density kg/mm^2

                    // Material I: Isotropic Resin
                    myMaterials[2].setType(0); // Isotropic
                    myMaterials[2].setProp(0,10000000.); // - E,  Young's Modulus
                    myMaterials[2].setProp(1,0.35); // - nu,   Poisson Ratio
                    myMaterials[2].setDensity(1.57e-5); // - rho,  Density  kg/mm^2
                }


                double inline FailureCriteria(Dune::FieldVector<double,6>& stress) const{
                    return stress[0];
                }

                //defaults to resing if not overwritten
                void assignMaterials(std::vector<int>& PGin_){
                    elemIndx2PG.resize(PGin_.size());  //material number following the numbering above
                    rot123.resize(PGin_.size());       //orientation of the material
                    for(int id = 0; id < PGin_.size(); id++){
                        rot123[id] = 0;
                        elemIndx2PG[id] = 1;
                    } // End for-loop of each element */
                }

                template <class GV>
                    void computeTensorsOnGrid(const GV& gv){
                        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
                        const typename GV::IndexSet& is = gv.indexSet();
                        C.resize(gv.size(0));

                        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it){ // loop through each element
                            int id = is.index(*it); // Get element id

                            int pg = elemIndx2PG[id];
                            C[id] = TensorRotation(baseTensor[pg],rot123[id],2);
                        } // End for-loop of each element */
                    }; // End Constructor


                template <class GV, class GV2>
                    void computeTensorsOnGrid(const GV& gv, const GV2& gv2){

                        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
                        const typename GV::IndexSet& is = gv.indexSet();
                        C.resize(gv.size(0));

                        //  Calculate Cijkl[id] in each element
                        auto it2 = gv2.template begin<0>();

                        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
                        { // loop through each element
                            int id = is.index(*it); // Get element id
                            int pg = elemIndx2PG[id];
                            C[id] = TensorRotation(baseTensor[pg],rot123[id],2);
                            it2++;
                        } // End for-loop of each element */
                    }; // End Constructor

                inline FieldVec gridTransformation(const FieldVec & x) const{
                    // Default grid transformation is the identity map
                    auto y = x;
                    return y;
                }

                std::array<int, 3> inline GridPartition(){
                    // Function returns default partition
                    std::array<int, 3> Partitioning;
                    std::fill(Partitioning.begin(),Partitioning.end(), 1);
                    return Partitioning;
                }

                void inline setPG(const YGV& ygv){
                    elemIndx2PG.resize(ygv.size(0));
                    rot123.resize(ygv.size(0));
                    // Loop Through Each element and Assign Physical Group and Rotation
                    for (const auto& eit : elements(ygv)){
                        int id = ygv.indexSet().index(eit);
                        auto point = eit.geometry().center();
                        // Which Layer?
                        int theRegion = whichRegion(point);
                        int theLayer = whichLayer(point,theRegion);
                        elemIndx2PG[id] = R[theRegion].layers[theLayer].getMaterial();
                        rot123[id] = R[theRegion].layers[theLayer].getOrientation();
                        for(int i = 0; i < eit.geometry().corners(); i++){  //Loop over nodes of the element
                            auto pc = eit.geometry().corner(i);
                            if(isMPC(pc)){ //Check for mpc constraints
                                elemIndx2PG[id] = 2;
                                rot123[id] = 0;
                            }
                        }
                    }
                }

                int inline getMaterialType(int whichMaterial){
                    assert(whichMaterial < numMaterials);
                    return myMaterials[whichMaterial].Type;
                }

                int inline setPhysicalGroup(FieldVec& x){
                    return 0;
                }

                void inline computeElasticTensors(){
                    // Calculate Elastic Tensors for each material
                    baseTensor.resize(numMaterials);
                    for (int i = 0; i < numMaterials; i++){
                        baseTensor[i] = Cijkl(myMaterials[i]);
                    }
                }


                Tensor4 TensorRotation(Tensor4& baseC, double theta, int i_dim) const
                {
                    Tensor4 locC = baseC;
                    Tensor4 R2(0.0),R2T(0.0);

                    //Set up 3D rotation matrix
                    double c = cos(theta*M_PI/180.0);
                    double s = sin(theta*M_PI/180.0);
                    Dune::FieldMatrix<double,3,3> A;
                    if(i_dim == 2) A = {{c,s,0},{-s,c,0},{0,0,1}};
                    if(i_dim == 1) A = {{c,0,-s},{0,1,0},{s,0,c}};
                    if(i_dim == 0) A = {{1,0,0},{0,c,s},{0,-s,c}};

                    //Transform to tensor rotation and transpose
                    R2 = getMatrix(A);
                    transpose<Tensor4,6>(R2,R2T);

                    //rotate the input matrix
                    locC.leftmultiply(R2);
                    locC.rightmultiply(R2T);
                    return locC;
                } // end eval

                Tensor4 TensorRotationLeft(Tensor4& C_local, double theta, int i_dim) const{
                    Tensor4 locC = C_local;
                    Tensor4 R2(0.0),R2T(0.0);

                    //Set up 3D rotation matrix
                    double c = cos(theta*M_PI/180.0);
                    double s = sin(theta*M_PI/180.0);
                    Dune::FieldMatrix<double,3,3> A;
                    if(i_dim == 2) A = {{c,s,0},{-s,c,0},{0,0,1}};
                    if(i_dim == 1) A = {{c,0,-s},{0,1,0},{s,0,c}};
                    if(i_dim == 0) A = {{1,0,0},{0,c,s},{0,-s,c}};

                    //Transform to tensor rotation
                    R2 = getMatrix(A);

                    //rotate the input matrix
                    locC.leftmultiply(R2);
                    return locC;
                }

                Tensor4 TensorRotationRight(Tensor4& C_local, double theta, int i_dim) const{
                    Tensor4 locC = C_local;
                    Tensor4 R2(0.0),R2T(0.0);

                    //set up 3D rotation matrix
                    double c = cos(theta*M_PI/180.0);
                    double s = sin(theta*M_PI/180.0);
                    Dune::FieldMatrix<double,3,3> A;
                    if(i_dim == 2) A = {{c,s,0},{-s,c,0},{0,0,1}};
                    if(i_dim == 1) A = {{c,0,-s},{0,1,0},{s,0,c}};
                    if(i_dim == 0) A = {{1,0,0},{0,c,s},{0,-s,c}};

                    //Transform to tensor rotation and transpose
                    R2 = getMatrix(A);
                    transpose<Tensor4,6>(R2,R2T);

                    //rotate the input matrix
                    locC.rightmultiply(R2T);
                    return locC;
                }

            private:
                //Turn 3D rotation matrix into a tensor rotation
                Tensor4 getMatrix(Dune::FieldMatrix<double,3,3> A) const{
                    Tensor4 R;
                    R[0][0] = A[0][0]*A[0][0]; R[0][1] = A[0][1]*A[0][1]; R[0][2] = A[0][2]*A[0][2];
                    R[1][0] = A[1][0]*A[1][0]; R[1][1] = A[1][1]*A[1][1]; R[1][2] = A[1][2]*A[1][2];
                    R[2][0] = A[2][0]*A[2][0]; R[2][1] = A[2][1]*A[2][1]; R[2][2] = A[2][2]*A[2][2];

                    R[0][3] = 2.0*A[0][1]*A[0][2]; R[0][4] = 2.0*A[0][0]*A[0][2]; R[0][5] = 2.0*A[0][0]*A[0][1];
                    R[1][3] = 2.0*A[1][1]*A[1][2]; R[1][4] = 2.0*A[1][0]*A[1][2]; R[1][5] = 2.0*A[1][0]*A[1][1];
                    R[2][3] = 2.0*A[2][1]*A[2][2]; R[2][4] = 2.0*A[2][0]*A[2][2]; R[2][5] = 2.0*A[2][0]*A[2][1];

                    R[3][0] = A[1][0]*A[2][0]; R[3][1] = A[1][1]*A[2][1]; R[3][2] = A[1][2]*A[2][2];
                    R[4][0] = A[0][0]*A[2][0]; R[4][1] = A[0][1]*A[2][1]; R[4][2] = A[0][2]*A[2][2];
                    R[5][0] = A[0][0]*A[1][0]; R[5][1] = A[0][1]*A[1][1]; R[5][2] = A[0][2]*A[1][2];

                    R[3][3] = A[1][1]*A[2][2]+A[1][2]*A[2][1]; R[3][4] = A[1][0]*A[2][2]+A[1][2]*A[2][0]; R[3][5] = A[1][0]*A[2][1]+A[1][1]*A[2][0];
                    R[4][3] = A[0][1]*A[2][2]+A[2][1]*A[0][2]; R[4][4] = A[0][0]*A[2][2]+A[2][0]*A[0][2]; R[4][5] = A[2][0]*A[0][1]+A[2][1]*A[0][0];
                    R[5][3] = A[0][1]*A[1][2]+A[0][2]*A[1][1]; R[5][4] = A[0][0]*A[1][2]+A[0][2]*A[1][0]; R[5][5] = A[0][0]*A[1][1]+A[0][1]*A[1][0];
                    return R;
                }

                //void inline evaluateTensor(int id, Tensor4& C_){    C_ = C[id];    }
            public:

                Tensor4 inline getElasticTensor(int id) const { return C[id]; }

                std::vector<double> returnTheta(const FieldVec& x) const{
                    FieldVec y, yph;
                    double h = 1e-10;
                    FieldVec xph = x; xph[0] += h/2.;
                    FieldVec xmh = x; xmh[0] -= h/2.;
                    y = gridTransformation(xph);
                    yph = gridTransformation(xmh);
                    auto angle = atan2((yph[2]-y[2])/h,(yph[0]-y[0])/h);
                    std::vector<double> theta(3);
                    theta[1] = angle;
                    return theta;
                }

                //! A function which returns whether a point x lies on a Dirichlet boundary for the degree of freedom dof.
                //! If the function returns false it assumes this node is on the Neumann boundary or an interior node.
                bool inline isDirichlet(FieldVec& x, int dof){
                    // Default isDirichlet function, returning homogeneous Dirichlet boundary conditioners
                    if (x[0] < 1e-6)
                        return true;
                    return false;
                }

                //! A function which returns whether a point is constrained by a multi-point-constraint
                bool inline isMPC(FieldVec& x){
                    return false;
                }

                //! returns the Dirichlet boundary value at a point for the degree of freedom dof.
                double inline evaluateDirichlet(const FieldVec& x, int dof) const{
                    // Return homogeneous boundary conditions
                    return 0.0;
                }

                double inline getMaterialTypeFromElement(int id) const{
                    int pg = elemIndx2PG[id];
                    return myMaterials[pg].Type;
                }

                double inline getOrientation(int id) const {return rot123[id];}

                template<class GO, class U, class GFS, class C,class MBE, class GV>
                    void inline postprocess(const GO& go, U& u, const GFS& gfs, const C& cg, const GV gv, MBE& mbe){
                        // Default does nothing
                    }

                template<class GO, class GO_EXT, class V, class GFS, class GFS_EXT, class C, class CON, class MBE, class LOP>
                    bool inline solve(GO& go, V& u, const GFS& gfs, const GFS_EXT& gfs_ext, C& cg, CON& constraints, MBE& mbe, LOP& lop){

                        bool flag = true;
                        if(solverParameters.UG == true){
                            typedef Dune::PDELab::ISTLBackend_NOVLP_CG_AMG_SSOR<GO> NOVLP_AMG;
                            NOVLP_AMG ls(go,500,2);
                            Dune::PDELab::StationaryLinearProblemSolver<GO,NOVLP_AMG,V> slp(go,ls,u,solverParameters.KrylovTol);
                            slp.apply();
                        }
                        else{
                            if(gfs.gridView().comm().size() > 1){ // == DEFAULT AVAILABLE PARALLEL SOLVERS
                                if(solverParameters.preconditioner == "GenEO"){
#if HAVE_SUITESPARSE and HAVE_ARPACKPP
                                    if(gfs.gridView().comm().rank() == 0) std::cout << "Starting solve with Geneo Preconditioner" << std::endl;
                                    typedef Dune::Geneo::CG_GenEO<V,GFS,GFS_EXT,GO,GO_EXT,LOP,CON,MBE> GenEO;
                                    int cells = std::pow(gfs.gridView().size(0),1./3.);
                                    GenEO mysolver(u,gfs,gfs_ext,lop,constraints,mbe,solverParameters,helper, cells);
                                    mysolver.apply();
                                    return flag;
#else
                                    if(gfs.gridView().comm().rank() == 0) std::cout << "GenEO solver requires UMFPack and Arpack." << std::endl;
#endif
                                }
                                if(solverParameters.preconditioner == "OneLevel_AdditiveSchwarz"){
#if HAVE_SUITESPARSE
                                    typedef Dune::PDELab::ISTLBackend_OVLP_CG_UMFPack<GFS, C> PAR_UMF;  //CG with UMFPack
                                    PAR_UMF ls(gfs,cg,solverParameters.MaxIt,solverParameters.verb);
                                    Dune::PDELab::StationaryLinearProblemSolver<GO,PAR_UMF,V> slp(go,ls,u,solverParameters.KrylovTol);
                                    slp.apply();
                                    return flag;
#else
                                    if(gfs.gridView().comm().rank() == 0) std::cout << "PAR_UMFPack solver not available" << std::endl;
#endif
                                }
                            }
                            if(solverParameters.solver == "Hypre"){
#if HAVE_HYPRE
                                typedef typename GO::Jacobian Matrix;
                                Matrix a(go);
                                auto b = u;
                                go.residual(u,b);
                                b *= -1.0;
                                HypreSolver<GFS, Matrix, V> solver(gfs,u,a,b, solverParameters.KrylovTol);
                                go.jacobian(u,a);
                                solver.solve(u);
                                return flag;
#else
                                if(gfs.gridView().comm().rank() == 0) std::cout << "Hypre solver not available" << std::endl;
#endif
                            }
                            if(gfs.gridView().comm().size() == 1){ // == SEQUENTIAL SOLVER
                                if(solverParameters.solver == "UMFPack"){
#if HAVE_SUITESPARSE
                                    typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack SEQ_UMF;
                                    SEQ_UMF ls(solver_verb); //Direct solver
                                    Dune::PDELab::StationaryLinearProblemSolver<GO, SEQ_UMF,V> slp(go,ls,u,1e-6); // note tolerance is ignored as Direct solver
                                    slp.apply();
#else
                                    if(gfs.gridView().comm().rank() == 0) std::cout << "SEQ_UMFPack solver not available" << std::endl;
#endif
                                    return flag;
                                }
                                if(solverParameters.preconditioner == "AMG"){
                                    typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> SEQ_CG_AMG_SSOR;
                                    SEQ_CG_AMG_SSOR ls(500,verbosity);
                                    Dune::PDELab::StationaryLinearProblemSolver<GO,SEQ_CG_AMG_SSOR,V> slp(go,ls,u,tolerance);
                                    slp.apply();
                                }
                            }
                        }

                        return flag;
                    }


            private:
                bool GenEO;

        };

    }
}

#endif
