

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>


#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
#include "MixedEventSkimmer/MixedEventAnalyzer.h"
#include "CorrectedSpectra/CorrectedSpectraMaker.h"

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);


	TaskFactory::registerTaskRunner<FemtoDstSkimmer>( "FemtoDstSkimmer" );
	TaskFactory::registerTaskRunner<MixedEventAnalyzer>( "MixedEventAnalyzer" );
	TaskFactory::registerTaskRunner<CorrectedSpectraMaker>( "CorrectedSpectraMaker" );

	TaskEngine engine( argc, argv );

	return 0;
}
