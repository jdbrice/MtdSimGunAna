/*
* @Author: Daniel
* @Date:   2017-05-05 14:22:49
* @Last Modified by:   Daniel
* @Last Modified time: 2017-05-26 10:16:12
*/
#ifndef TRACK_FILTER_H
#define TRACK_FILTER_H

#include "XmlConfig.h"
#include "XmlRange.h"

// FemtoDst
#include "FemtoDstFormat/FemtoTrackProxy.h"

// STL
#include <map>

#include "vendor/loguru.h"

class TrackFilter {
public:
	TrackFilter() {}
	TrackFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~TrackFilter() {}

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

		if ( range.count( "pT" ) && !range[ "pT" ].inInclusiveRange( _proxy._track->mPt ) )
			return false;
		if ( range.count( "eta" ) && !range[ "eta" ].inInclusiveRange( _proxy._track->mEta ) )
			return false;
		if ( range.count( "nHitsFit" ) && !range[ "nHitsFit" ].inInclusiveRange( abs(_proxy._track->mNHitsFit) ) )
			return false;
		if ( range.count( "nHitsDedx" ) && !range[ "nHitsDedx" ].inInclusiveRange( _proxy._track->mNHitsDedx) )
			return false;
		if ( range.count( "nHitsRatio" ) && !range[ "nHitsRatio" ].inInclusiveRange( nHitsRatio ) )
			return false;
		if ( range.count( "nSigmaPion" ) && !range["nSigmaPion"].inInclusiveRange( _proxy._track->nSigmaPion() ) )
			return false;

		if ( nullptr != _proxy._helix && range.count( "dca" ) && !range["dca"].inInclusiveRange( _proxy._helix->mDCA ) )
			return false;


		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}

};

#endif