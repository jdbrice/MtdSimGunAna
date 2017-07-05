#ifndef PID_LR_H
#define PID_LR_H

#include "IObject.h"
#include "XmlConfig.h"
#include "XmlFunction.h"

using namespace jdb;

#include <map>
#include <string>
#include <limits>
using namespace jdb;


// Project
#include "CandidateTrack.h"
#include "CandidateTrackMtdPidTraits.h"

class PidLR : public IObject
{
public:
	virtual const char* classname() const { return "PidLR"; }
	PidLR() {}
	~PidLR() {}

	void addPDF( XmlConfig &_cfg, string nodePath ) {

		XmlFunction xf1;
		xf1.set( _cfg, nodePath );
		string name = _cfg.getXString( nodePath + ":name" );

		pdfRatio[ name ] = xf1;

		INFOC( "Loaded: [" << name <<"] == " << xf1.toString() );
		xf1.getTF1()->Write();
	}


	double eval( CandidateTrack *_track, CandidateTrackMtdPidTraits * _mtd ){
		if ( nullptr == _track || nullptr == _mtd ) return std::numeric_limits<double>::quiet_NaN();

		double Y = 1.0;
		for ( auto nf : pdfRatio ){
			DEBUGC( "name " << nf.first );
			DEBUGC( "(before) Y = " << Y );
			if ("nSigmaPion" == nf.first ){
				DEBUGC( "x = " << _track->mNSigmaPion );
				double y = nf.second.eval( _track->nSigmaPion() );
				if ( y < 0 ) y = 100;	// something wrong, must be background
				Y *= y;
			}
			if ("dca" == nf.first ){
				DEBUGC( "x = " << _track->gDCA() );
				double y = nf.second.eval( _track->gDCA() );
				if ( y < 0 ) y = 100;	// something wrong, must be background
				Y *= y;
			}
			if ("mtdCell" == nf.first ){
				DEBUGC( "x = " << _mtd->cell( ) );
				double y = nf.second.eval( _mtd->cell( ) );
				if ( y < 0 ) y = 100;	// something wrong, must be background
				Y *= y;
			}
			if ("mtdDeltaY" == nf.first ){
				DEBUGC( "x = " << _mtd->mDeltaY * _track->charge() );
				double y = nf.second.eval( _mtd->mDeltaY * _track->charge() );
				if ( y < 0 ) y = 100;	// something wrong, must be background
				Y *= y;
			}
			if ("mtdDeltaZ" == nf.first ){
				DEBUGC( "x = " << _mtd->mDeltaZ );
				double y = nf.second.eval( _mtd->mDeltaZ );
				if ( y < 0 ) y = 100;	// something wrong, must be background
				Y *= y;
			}
			if ("mtdDeltaTOF" == nf.first ){
				DEBUGC( "x = " << _mtd->mDeltaTimeOfFlight );
				double y = nf.second.eval( _mtd->mDeltaTimeOfFlight );
				if ( y < 0 ) y = 100;	// something wrong, must be background
				Y *= y;
			}
			DEBUGC( "(now) Y = " << Y );
		}

		DEBUGC( "Y = " << Y );
		double LR = ( 1.0 - Y ) / ( 1.0 + Y );
		DEBUGC( "LR = " << LR );
		return LR;
	}

	int nPdfs(  ){
		return pdfRatio.size();
	}




protected:

	map<string, XmlFunction> pdfRatio;



	
};


#endif