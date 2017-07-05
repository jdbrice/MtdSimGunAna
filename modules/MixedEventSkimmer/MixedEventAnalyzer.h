#ifndef MIXED_EVENT_ANALYZER_H
#define MIXED_EVENT_ANALYZER_H

#include "MixedEventSkimmer.h"

// Analyzers
#include "FemtoDstSkimmer/PairHistogramMaker.h"
#include "FemtoDstSkimmer/TrackHistogramMaker.h"
#include "FemtoDstSkimmer/MtdHistogramMaker.h"

class MixedEventAnalyzer : public MixedEventSkimmer
{
protected:

	PairHistogramMaker phm;
	TrackHistogramMaker thm;
	MtdHistogramMaker mhm;

	bool trackQA, pairQA, mtdQA;

public:
	virtual const char* classname() const { return "MixedEventAnalyzer"; }
	MixedEventAnalyzer() {}
	~MixedEventAnalyzer() {}

	virtual void initialize(){
		MixedEventSkimmer::initialize();

		trackQA = config.getBool( "TrackHistogramMaker:enable", true );
		pairQA  = config.getBool( "PairHistogramMaker:enable", true );
		mtdQA   = config.getBool( "MtdHistogramMaker:enable", true );

		if ( pairQA )  phm.setup( config, "PairHistogramMaker", book );
		if ( trackQA ) thm.setup( config, "TrackHistogramMaker", book );
		if ( mtdQA )   mhm.setup( config, "MtdHistogramMaker", book );

		LOG_IF_F( INFO, pairQA,  "PairHistogramMaker ENABLED" );
		LOG_IF_F( INFO, trackQA, "TrackHistogramMaker ENABLED" );
		LOG_IF_F( INFO, mtdQA,   "MtdHistogramMaker ENABLED" );
	}

protected:

	virtual void analyzePairs( FemtoTrackContainer &tc, vector<FemtoTrackContainer> &buff ){

		if ( trackQA ) thm.analyze( tc._proxy );
		if ( mtdQA )   mhm.analyze( tc._proxy );

		if ( pairQA ){
			
			phm.setEvent( _rEvent.get() );
			for ( auto btc : buff ){
				if ( tc._runId == btc._runId && tc._eventId == btc._eventId ) continue;
				phm.analyze( tc._proxy, btc._proxy );
			}

			// phm.fillAggregates();
		}
		
	}
	
};


#endif