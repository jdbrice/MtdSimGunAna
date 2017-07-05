#ifndef TRACK_HISTOGRAM_MAKER_H
#define TRACK_HISTOGRAM_MAKER_H

#include "IObject.h"
#include "HistoBook.h"
#include "XmlConfig.h"

#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include <memory>

#include "vendor/loguru.h"

class TrackHistogramMaker : public IObject
{
protected:
	XmlConfig _xfg;
	shared_ptr<HistoBook> _book = nullptr;
	string _nodePath;

	size_t nPos, nNeg;

public:
	virtual const char* classname() const { return "TrackHistogramMaker"; }
	TrackHistogramMaker() {}
	TrackHistogramMaker( XmlConfig &_xfg, string _np, shared_ptr<HistoBook> _book ) {
		setup( _xfg, _np, _book );
	}
	~TrackHistogramMaker() {}

	void setup( XmlConfig &_xfg, string _np, shared_ptr<HistoBook> _book ) {
		this->_xfg      = _xfg;
		this->_nodePath = _np;
		this->_book     = _book;

		if ( nullptr != this->_book ){
			LOG_F( INFO, "Making TrackHistogramMaker histos @ %s", (this->_nodePath + ".histograms").c_str() );
			this->_book->cd();
			this->_book->makeAll( this->_xfg, this->_nodePath + ".histograms" );
		}

		nPos = nNeg = 0;
	}


	void analyze( FemtoTrackProxy &_tp ){
		fill( _tp, "" );
		if ( _tp._track->charge() > 0 ){
			nPos++;
			fill( _tp, "pos_" );
		}
		else {
			nNeg++;
			fill( _tp, "neg_" );
		}
	}

	void fill( FemtoTrackProxy &_tp, string prefix ){
		string dir = "track/";
		_book->fill ( dir + prefix + "mPt", _tp._track->mPt );
		_book->fill ( dir + prefix + "mEta", _tp._track->mEta );
		_book->fill ( dir + prefix + "mPhi", _tp._track->mPhi );
		_book->fill ( dir + prefix + "mId", _tp._track->mId );
		_book->fill ( dir + prefix + "mDedx", _tp._track->dEdx() );
		_book->fill ( dir + prefix + "mNHitsFit", _tp._track->mNHitsFit );
		_book->fill ( dir + prefix + "mNHitsMax", _tp._track->mNHitsMax );
		_book->fill ( dir + prefix + "mNHitsDedx", _tp._track->mNHitsDedx );
		_book->fill ( dir + prefix + "mNSigmaPion", _tp._track->nSigmaPion() );
		_book->fill ( dir + prefix + "mNSigmaKaon", _tp._track->nSigmaKaon() );
		_book->fill ( dir + prefix + "mNSigmaProton", _tp._track->nSigmaProton() );
		_book->fill ( dir + prefix + "mNSigmaElectron", _tp._track->nSigmaElectron() );
		_book->fill ( dir + prefix + "mDCA", _tp._track->gDCA() );
		if ( _tp._helix )
			_book->fill ( dir + prefix + "mHelixDCA", _tp._helix->mDCA );
	}


	void fillAggregates(){
		_book->fill( "track/nNeg_nPos", nPos, nNeg );
		nPos = nNeg = 0;
	}
	
};




#endif