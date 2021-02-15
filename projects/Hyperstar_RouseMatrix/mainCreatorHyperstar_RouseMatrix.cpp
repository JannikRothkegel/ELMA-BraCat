

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
#include <LeMonADE/analyzer/AnalyzerWriteBfmFileSubGroup.h>

#include "UpdaterCreatorHyperstar.h"
#include "Analyzer_EigenvaluesRouseMatrix.h"


int main(int argc, char* argv[])
{
  try{
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
        // FeatureFeatureAttributes<> is equivalent to FeatureAttributes<int32_t>
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 500> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	std::string filename="Hyperstar.bfm";
        uint32_t maxfunctionality = 7;
        
	int option_char(0);
	int box(128);
        
        int numSamples=1;

	//read in options by getopt
	while ((option_char = getopt (argc, argv, "f:g:b:s:h"))  != EOF){
		switch (option_char)
		{
		case 'f':
			filename=optarg;
			break;
		case 'g':
			maxfunctionality = atof(optarg);
			break;
		case 'b':
			box = atoi(optarg);
			break;
                case 's':
			numSamples = atoi(optarg);
			break;        
		case 'h': std::cout << "This program creates a hyperstar branched structure." << std::endl;
                          //std::cout << "It uses no excluded volume condition for creation anf functionality of nodes is limited to 500." << std::endl;
                          std::cout << "It evaluates the ideal Gaussian radius of gyration" << std::endl;
                          //std::cout << "Option 'p' represents the walking probability proportional to the walking rate~walking steps/time between reaction event" << std::endl;
                          //std::cout << "Option 't' represents the walking probability over tertiary bond e.g. t=1.0 no restriction and t=0.0 full barrier." << std::endl;
                
		default:
			std::cerr << "Usage: " << argv[0] << " [-f filename] [-g generation(=7)] [-s number_of_samples(=1)] [-b box size] \n";
			return 0;
		}
	}
	std::cout << "write out options: filename=" << filename <<  " generation=" << maxfunctionality << " in box " << box <<std::endl;
        std::cout << "Number of samples for ensemble average: " << numSamples <<std::endl;

	//seed the globally available random number generators
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    myIngredients.setName(filename);


    UpdaterCreatorHyperstar<Ing> U_CW(myIngredients, maxfunctionality, box, box, box);

    Analyzer_EigenvaluesRouseMatrix<Ing> A_EVRM(myIngredients, ".");

        for(int i = 0; i < numSamples; i++)
        {
        U_CW.initialize();
        U_CW.execute();

        A_EVRM.initialize();
        A_EVRM.execute();

        }
        U_CW.cleanup();
        A_EVRM.cleanup();



	//TaskManager taskmanager;
	//taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	//taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	//taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	//taskmanager.initialize();
	//taskmanager.run(max_mcs/save_interval);
	//taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

