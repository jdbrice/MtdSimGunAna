#ifndef MTD_HISTOGRAM_MAKER_H
#define MTD_HISTOGRAM_MAKER_H

#include "IObject.h"
#include "HistoBook.h"
#include "XmlConfig.h"

#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include <memory>

#include "vendor/loguru.h"

class MtdHistogramMaker : public IObject
{
protected:
	XmlConfig _xfg;
	shared_ptr<HistoBook> _book = nullptr;
	string _nodePath;

public:
	virtual const char* classname() const { return "MtdHistogramMaker"; }
	MtdHistogramMaker() {}
	MtdHistogramMaker( XmlConfig &_xfg, string _np, shared_ptr<HistoBook> _book ) {
		setup( _xfg, _np, _book );
	}
	~MtdHistogramMaker() {}

	void setup( XmlConfig &_xfg, string _np, shared_ptr<HistoBook> _book ) {
		this->_xfg      = _xfg;
		this->_nodePath = _np;
		this->_book     = _book;

		if ( nullptr != this->_book ){
			LOG_F( INFO, "Making MtdHistogramMaker histos @ %s", (this->_nodePath + ".histograms").c_str() );
			this->_book->cd();
			this->_book->makeAll( this->_xfg, this->_nodePath + ".histograms" );
		}
	}


	void analyze( FemtoTrackProxy &_tp ){
		fill( _tp, "" );
		if ( _tp._track->charge() > 0 ){
			fill( _tp, "pos_" );
		}
		else {
			fill( _tp, "neg_" );
		}
	}

	void fill( FemtoTrackProxy &_tp, string prefix ){
		string dir = "mtd/";
		_book->fill( dir + prefix + "mDeltaY", _tp._mtdPid->mDeltaY );
		_book->fill( dir + prefix + "mDeltaZ", _tp._mtdPid->mDeltaZ );
		_book->fill( dir + prefix + "mDeltaTimeOfFlight", _tp._mtdPid->mDeltaTimeOfFlight );
		_book->fill( dir + prefix + "mMatchFlag", _tp._mtdPid->mMatchFlag );
		_book->fill( dir + prefix + "mMtdHitChan", _tp._mtdPid->mMtdHitChan );
		_book->fill( dir + prefix + "mTriggerFlag", _tp._mtdPid->mTriggerFlag );

		_book->fill( dir + prefix + "mBackleg", _tp._mtdPid->backleg() );
		_book->fill( dir + prefix + "mModule", _tp._mtdPid->module() );
		_book->fill( dir + prefix + "mCell", _tp._mtdPid->cell() );
		
	}


	void fillAggregates(){
	}
	
};




#endif