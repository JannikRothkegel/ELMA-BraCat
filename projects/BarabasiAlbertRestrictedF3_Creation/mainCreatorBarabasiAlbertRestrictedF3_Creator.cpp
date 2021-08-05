

#include <cstring>
#include <stdlib.h> //for atoi
#include <unistd.h> //for getopt

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/utility/DepthIteratorPredicates.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include "UpdaterCreatorBarabasiAlbertRestrictedF3.h"

int main(int argc, char* argv[]) {
    try {
        // FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
        // FeatureFeatureAttributes<> is equivalent to FeatureAttributes<int32_t>
        //typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>) Features;
        typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>, FeatureExcludedVolumeSc<>) Features;

        typedef ConfigureSystem<VectorInt3, Features, 3> Config;
        typedef Ingredients<Config> Ing;
        Ing myIngredients;

        std::string filename = "BarabasiAlbert.bfm";
        uint32_t numbermonomers = 128;

        int option_char(0);
        int box(128);

        //read in options by getopt
        while ((option_char = getopt(argc, argv, "f:n:b:h")) != EOF) {
            switch (option_char) {
                case 'f':
                    filename = optarg;
                    break;
                case 'n':
                    numbermonomers = atof(optarg);
                    break;
                case 'b':
                    box = atoi(optarg);
                    break;
                case 'h': std::cout << "This program creates a BarabasiAlbert branched structure with restricted func f <= 3 -> Dendrimer." << std::endl;
                    std::cout << "It only adds one node with exactly one edge to the graph. The creation of functionality of nodes is limited to 3." << std::endl;
                    std::cout << "The creation is done under full excluded volume conditions." << std::endl;

                default:
                    std::cerr << "Usage: " << argv[0] << " [-f filename] [-n numbermonomers(=128)] [-b box size] \n";
                    return 0;
            }
        }
        std::cout << "write out options: filename=" << filename << " maxmonomers=" << numbermonomers << " in box " << box << std::endl;
        
        //seed the globally available random number generators
        RandomNumberGenerators randomNumbers;
        randomNumbers.seedAll();

        myIngredients.setName(filename);

        UpdaterCreatorBarabasiAlbertRestrictedF3<Ing> U_BA(myIngredients, numbermonomers, box, box, box);

        U_BA.initialize();
        U_BA.execute();
        U_BA.cleanup();

        //output
        AnalyzerWriteBfmFile<Ing> Write(filename, myIngredients);
        Write.initialize();
        Write.execute();
        Write.cleanup();


        //TaskManager taskmanager;
        //taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
        //here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
        //(other than for latticeOccupation, valid bonds, frozen monomers...)
        //taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

        //taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));

        //taskmanager.initialize();
        //taskmanager.run(max_mcs/save_interval);
        //taskmanager.cleanup();

    } catch (std::exception& err) {
        std::cerr << err.what();
    }
    return 0;

}

