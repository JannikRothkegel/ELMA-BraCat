

#include <cstring>
#include <stdlib.h> //for atoi
#include <unistd.h> //for getopt

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
//#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

#include "UpdaterCreateDendrimerGeneral_FGS.h"


int main(int argc, char* argv[])
{
  try{
	//std::string infile;
	std::string outfile;
	uint32_t max_mcs=0;
	uint32_t save_interval=0;
	
	outfile="outfile.bfm";
	
	int option_char(0);

	int generation = 1;
	int functionality = 3;
	int spacer_length = 1;
	int box_size_x = 128;
	int box_size_y = 128;
	int box_size_z = 128;

		//read in options by getopt
		while ((option_char = getopt (argc, argv, "o:f:g:s:x:y:z:h"))  != EOF){
			switch (option_char)
			{

			case 'o':
					outfile=optarg;
					break;

			case 's':
					spacer_length = atoi(optarg);
					break;
			case 'f':
					functionality = atoi(optarg);

					if(functionality > 6)
						{
						throw std::runtime_error("ERROR: Functionality larger than 6. This is maybe not suitable for the BFM. Exiting...");
						}

					break;

			case 'g':
					generation = atoi(optarg);

				break;
			case 'x':
					box_size_x = atoi(optarg);
					break;
			case 'y':
					box_size_y = atoi(optarg);
					break;
			case 'z':
					box_size_z = atoi(optarg);
					break;

			case 'h':
			default:
				std::cerr << "General Creation of Dendrimer with Generation G, Spacer length S, Functionality F<7" << std::endl;
				std::cerr << "Usage: " << argv[0] << "  [-o filenameOutput] [-f functionality] [-g generation] [-s spacer length] [-x box size X] [-y box size Y] [-z box size Z] \n";
				return 0;
			}
		}
	std::cout << "write out options:  filenameOutput=" << outfile <<  " spacer length:" << spacer_length << " generation:" << generation << " functionality:" << functionality << " box size X" << box_size_x  << " box size Y" << box_size_y << " box size Z" << box_size_z << std::endl;


	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	

	// FeatureExcludedVolumeSc<> is equivalent to FeatureExcludedVolumeSc<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes<>,FeatureExcludedVolumeSc<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;
	//myIngredients.setName(outfile);

	UpdaterCreateDendrimerGeneral_FGS<Ing> UC_FGS(myIngredients, generation, spacer_length, functionality, box_size_x, box_size_y, box_size_z);
	UC_FGS.initialize();
	UC_FGS.execute();
	UC_FGS.cleanup();

	//taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>( outfile,myIngredients));
	AnalyzerWriteBfmFile<Ing> AWBFM( outfile,myIngredients);
	AWBFM.initialize();
	AWBFM.execute();
	AWBFM.cleanup();

	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

