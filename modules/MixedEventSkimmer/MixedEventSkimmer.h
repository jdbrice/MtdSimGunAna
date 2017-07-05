#ifndef MIXED_EVENT_SKIMMER_H
#define MIXED_EVENT_SKIMMER_H

#include "TreeAnalyzer.h"

#include "HistoBins.h"

#include "FemtoDstFormat/FemtoTrackContainer.h"
#include "Filters/TrackFilter.h"
#include "Filters/MtdFilter.h"

#include "vendor/loguru.h"

#include <map>
#include <tuple>

// ROOT
#include "TRandom3.h"
#include "TNamed.h"

class MixedEventSkimmer : public TreeAnalyzer
{
protected:
	
	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	size_t bufferSize;
	map< tuple<char, size_t, size_t>, vector<FemtoTrackContainer> > buffer;
	map< tuple<size_t, size_t, size_t>, bool> trackPresent;
	HistoBins h_vtxZ, h_gRefMult;

	vector<FemtoTrackContainer> nullBuffer;

	TRandom3 rander;

	bool rejectSingleTrack = false;
	TrackFilter _trackFilter;
	MtdFilter _mtdFilter;


public:
	virtual const char* classname() const { return "MixedEventSkimmer"; }
	MixedEventSkimmer() {}
	~MixedEventSkimmer() {}

	virtual void overrideConfig(){
		// set the number of Event Loops to 2!
		config.set( nodePath + ":nEventLoops", "2" );
	}
	

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		bufferSize = config.getInt( "bufferSize", 10 );

		h_vtxZ.load( config, nodePath + ".Mix.vtxZ" );
		h_gRefMult.load( config, nodePath + ".Mix.gRefMult" );

		LOG_F(1, "h_vtxZ: %s", h_vtxZ.toString().c_str() );
		LOG_F(1, "h_gRefMult: %s", h_gRefMult.toString().c_str() );

		rander.SetSeed(0);
		_trackFilter.load( config, nodePath + ".TrackFilter" );
		_mtdFilter.load( config, nodePath + ".MtdFilter" );

		rejectSingleTrack = config.getBool( nodePath + ".Mix:rejectSingleTrack", true );

	}

protected:


	vector<FemtoTrackContainer> &currentBuffer( int charge ){
		FemtoEvent *event = _rEvent.get();

		if ( nullptr == event ) return nullBuffer;
		if ( charge != -1 && charge != 1 ) return nullBuffer;

		size_t iVtx = h_vtxZ.findBin( event->mPrimaryVertex_mX3 );
		size_t iGRM = h_gRefMult.findBin( event->mGRefMult );

		// creates if DNE
		return buffer[ std::make_tuple( charge, iVtx, iGRM ) ];

	}

	void addToBuffer( vector<FemtoTrackContainer> &vtc, FemtoTrackContainer tc ){
		
		auto t = std::make_tuple( tc._runId, tc._eventId, tc._track->mId );
		if ( trackPresent.count( t ) > 0 && trackPresent[ t ] == true ) return;

		trackPresent[ t ] = true;
		if ( vtc.size() < bufferSize )
			vtc.push_back( tc );
		else {
			// replace in proper size of buffer
			size_t index = rander.Integer( bufferSize );
			trackPresent[ std::make_tuple( vtc[ index ]._runId, vtc[ index ]._eventId, vtc[ index ]._track->mId ) ] = false;
			vtc[ index ] = tc;
		}
	}

	bool allBuffersFull(  ){
		for (auto bf : buffer ){
			if ( bf.second.size() < bufferSize )
				return false;
		}
		return true;
	}


	virtual void analyzePairs( FemtoTrackContainer &tc, vector<FemtoTrackContainer> &buff ){

	}


	virtual void analyzeEvent(){
		FemtoEvent *event = _rEvent.get();
		if ( nullptr == event ) return;

		vector<FemtoTrackContainer> &pos = currentBuffer( 1 );
		vector<FemtoTrackContainer> &neg = currentBuffer( -1 );

		FemtoTrackContainer tc;
		tc._runId = event->mRunId;
		tc._eventId = event->mEventId;

		size_t nTracks = _rTracks.N();
		if ( rejectSingleTrack && nTracks < 2 ) return;

		for (size_t i = 0; i < nTracks; i++ ){
			tc.assemble( i, _rTracks, _rMtdPid );

			if ( _trackFilter.fail( tc._proxy ) ) continue;
			if ( _mtdFilter.fail( tc._proxy ) ) continue;

			if ( tc._track->charge() == 1 ){
				if ( iEventLoop == 1 ){
					analyzePairs( tc, neg ); // ULS pairs
					analyzePairs( tc, pos ); // LS pairs
				}
				addToBuffer( pos, tc );
				
			}
			else if ( tc._track->charge() == -1 ){
				
				if ( iEventLoop == 1 ){
					analyzePairs( tc, pos ); // ULS pairs
					analyzePairs( tc, neg ); // LS pairs
				}
				addToBuffer( neg, tc );
			}
		}

		if ( iEventLoop == 0 && allBuffersFull() ){
			LOG_F( INFO, "Buffers full, skipping to next loop, iEvent=%llu", iEvent );
			finishCurrentEventLoopEarly = true;
		}


		// analyzeBuffers( pos, neg );
	}

	virtual void postEventLoop(){
		LOG_F( 0, "# of Buffers: %lu", buffer.size() );

		LOG_F( 0, "Entries in buffer[vtx, grm ]" );
		for ( auto kv : buffer ){
			LOG_F( 0, "Entries in buffer[ %d, %lu, %lu ] = %lu", std::get<0>(kv.first), std::get<1>(kv.first), std::get<2>(kv.first), kv.second.size() );
			// LOG_F( 0, "RunId=%lu, EvnetId=%lu", kv.second[0]._runId, kv.second[0]._eventId );
		}

		
		if ( iEventLoop == 1 ){
			book->clone( "lsp_pT_mass", "ls_pT_mass" );
			book->get( "ls_pT_mass" )->Add( book->get( "lsn_pT_mass" ) );
		}
		if ( config.getInt("jobIndex") == 0 && iEventLoop == 1 ){

			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
		
	}
};



#endif