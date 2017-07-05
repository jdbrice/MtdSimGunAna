#ifndef MTD_FILTER_H
#define MTD_FILTER_H
#include "XmlConfig.h"
#include "XmlRange.h"

// FemtoDst
#include "FemtoDstFormat/FemtoTrackProxy.h"

// STL
#include <map>
#include <set>

#include "vendor/loguru.h"

class MtdFilter
{
public:
	MtdFilter() {}
	MtdFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~MtdFilter() {}

	map<string, XmlRange> range;
	vector<string> failCutsList;
	set<string> failCuts;
	bool inverted = false;

	void load( XmlConfig &_cfg, string _nodePath ){
		vector<string> rngPaths=_cfg.childrenOf( _nodePath, "Range" );
		for ( auto p : rngPaths ){
			XmlRange rng( _cfg, p );

			range[ _cfg[ p + ":name" ] ] = rng;
			LOG_F( INFO, "[%s] = %s", _cfg[(p+":name")].c_str(), rng.toString().c_str() );
		}

		failCutsList = _cfg.getStringVector( _nodePath + ":fail" );
		for ( auto fc : failCutsList ){
			failCuts.insert( fc );
		}
	}



	vector<string> passPID( FemtoTrackProxy &_proxy ){
		vector<string> failedCuts;

		if ( range.count( "deltaY" ) ){
			bool in_range = range[ "deltaY" ].inInclusiveRange( _proxy._mtdPid->mDeltaY * _proxy._track->charge() );
			if ( !in_range ) {
				failedCuts.push_back( "deltaY" );
			}
		}
			
		if ( range.count( "deltaZ" ) ){
			bool in_range = range[ "deltaZ" ].inInclusiveRange( _proxy._mtdPid->mDeltaZ );
			if ( !in_range ) {
				failedCuts.push_back( "deltaZ" );
			}
		}
		
		if ( range.count( "deltaTimeOfFlight" ) ){
			bool in_range = range[ "deltaTimeOfFlight" ].inInclusiveRange( _proxy._mtdPid->mDeltaTimeOfFlight );
			if ( !in_range ) {
				failedCuts.push_back( "deltaTimeOfFlight" );
			}
		}

		if ( range.count( "cell" ) ){
			bool in_range = range[ "cell" ].inInclusiveRange( _proxy._mtdPid->cell() );
			if ( !in_range ) {
				failedCuts.push_back( "cell" );
			}
		}


		return failedCuts;
	}

	bool pass( FemtoTrackProxy &_proxy ){
		
		
		if ( range.count( "matchFlag" ) && !range[ "matchFlag" ].inInclusiveRange( _proxy._mtdPid->mMatchFlag ) )
			return false;
		if ( range.count( "triggerFlag" ) && !range["triggerFlag"].inInclusiveRange( _proxy._mtdPid->mTriggerFlag ) )
			return false;


		// check PID
		vector<string> fCuts = passPID( _proxy );
		if ( fCuts.size() != failCuts.size() )
			return false;

		for ( auto fCut : fCuts ){
			if ( failCuts.count( fCut ) == 0 )
				return false;
		}


		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}
	
};


#endif