/*
* @Author: Daniel
* @Date:   2017-05-05 14:22:49
* @Last Modified by:   Daniel
* @Last Modified time: 2017-05-26 10:16:12
*/
#ifndef LOW_PT_MUON_FILTER_H
#define LOW_PT_MUON_FILTER_H

#include "XmlConfig.h"
#include "XmlRange.h"

// FemtoDst
#include "FemtoDstFormat/FemtoTrackProxy.h"

// STL
#include <map>

#include "vendor/loguru.h"

class LowPtMuonFilter {
public:
	LowPtMuonFilter() {}
	LowPtMuonFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~LowPtMuonFilter() {}

	void load( XmlConfig &_cfg, string _nodePath ){
		vector<string> rngPaths=_cfg.childrenOf( _nodePath, "Range" );
		for ( auto p : rngPaths ){
			XmlRange rng( _cfg, p );

			range[ _cfg[ p + ":name" ] ] = rng;
			LOG_F( INFO, "[%s] = %s", _cfg[(p+":name")].c_str(), rng.toString().c_str() );
		}

	}


	map<string, XmlRange> range;

	bool pass( FemtoTrackProxy &_proxy ){

		double nHitsRatio = fabs(_proxy._track->mNHitsFit) / _proxy._track->mNHitsMax;

		if ( range.count( "nSigmaPion" ) && !range["nSigmaPion"].inInclusiveRange( _proxy._track->nSigmaPion() ) )
			return false;

		


		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}

};

#endif