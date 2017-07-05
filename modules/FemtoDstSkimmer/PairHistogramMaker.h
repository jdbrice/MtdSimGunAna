#ifndef PAIR_HISTOGRAM_MAKER_H
#define PAIR_HISTOGRAM_MAKER_H

#include "IObject.h"
#include "HistoBook.h"
#include "XmlConfig.h"

#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include <memory>

#include "vendor/loguru.h"

class PairHistogramMaker : public IObject
{
protected:
	XmlConfig _xfg;
	shared_ptr<HistoBook> _book = nullptr;
	string _nodePath;

	FemtoEvent *_event = nullptr;

	// aggregate quantities
	size_t nLsp, nLsn, nUls;
public:
	virtual const char* classname() const { return "PairHistogramMaker"; }
	PairHistogramMaker() {
		this->_nodePath = "";
	}
	PairHistogramMaker( XmlConfig &_xfg, string _np, shared_ptr<HistoBook> _book ) {
		setup( _xfg, _np, _book );
	}
	~PairHistogramMaker() {}

	void setup( XmlConfig &_xfg, string _np, shared_ptr<HistoBook> _book ) {
		this->_xfg      = _xfg;
		this->_nodePath = _np;
		this->_book     = _book;

		if ( nullptr != this->_book ){
			LOG_F( INFO, "Making PairHistogramMaker histos @ %s", (this->_nodePath + ".histograms").c_str() );
			this->_book->cd();
			this->_book->makeAll( this->_xfg, this->_nodePath + ".histograms" );
		}

		nLsn = nLsp = nUls = 0;
	}

	void setEvent( FemtoEvent *_event ){
		this->_event = _event;
	}

	void analyze( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2 ){
		if ( nullptr == this->_book) return;

		// if ( _tp1._track->mPt < 1.55 || _tp2._track->mPt < 1.55 ) return;

		int chargeSum = _tp1._track->charge() + _tp2._track->charge();
		FemtoTrackProxy *ltp = &_tp1;
		FemtoTrackProxy *sltp = &_tp2;

		if ( fabs(_tp2._track->mPt) > fabs( _tp1._track->mPt ) ){
			ltp  = &_tp2;
			sltp = &_tp1;
		}
		
		fill( *ltp, *sltp, "" );
		if ( chargeSum == 2 ){
			fill( *ltp, *sltp, "lsp_" );
			nLsp++;
		}
		else if ( chargeSum == -2 ){
			fill( *ltp, *sltp, "lsn_" );
			nLsn++;
		}
		else if ( chargeSum == 0 ){
			fill( *ltp, *sltp, "uls_" );
			nUls++;
		}


	}

	void fill( 	FemtoTrackProxy &_ltp, FemtoTrackProxy &_sltp, string prefix = "" ){

		TLorentzVector lv1 = _ltp._track->lv( 0.105 );
		TLorentzVector lv2 = _sltp._track->lv( 0.105 );
		TLorentzVector lv = lv1 + lv2;


		_book->fill( prefix + "mass", lv.M() );
		if ( "lsp_" == prefix || "lsn_" == prefix )
			_book->fill( "ls_mass", lv.M() );
		
		// if ( 1 == _ltp._mtdPid->mTriggerFlag && 1 == _sltp._mtdPid->mTriggerFlag )
		// 	_book->fill( prefix + "mass_tf1_1", lv.M() );
		// else if ( _ltp._mtdPid->mTriggerFlag > 1 && _sltp._mtdPid->mTriggerFlag > 1 )
		// 	_book->fill( prefix + "mass_tf1_1", lv.M() );
		// else if ( (1 == _ltp._mtdPid->mTriggerFlag && _sltp._mtdPid->mTriggerFlag > 1) || (_ltp._mtdPid->mTriggerFlag > 1 && 1 == _sltp._mtdPid->mTriggerFlag) )
		// 	_book->fill( prefix + "mass_tf1_N", lv.M() );

		_book->fill( prefix + "pT", lv.Pt() );
		_book->fill( prefix + "pT_lpT", lv1.Pt(), lv.Pt() );
		_book->fill( prefix + "pT_slpT", lv2.Pt(), lv.Pt() );

		_book->fill( prefix + "pT_mass", lv.M(), lv.Pt() );
		if ( "lsp_" == prefix || "lsn_" == prefix )
			_book->fill( "ls_pT_mass", lv.M(), lv.Pt() );
		if ( nullptr != this->_event ){
			_book->fill( prefix + "vtxZ_mass", lv.M(), this->_event->mPrimaryVertex_mX3 );
			_book->fill( prefix + "gRefMult_mass", lv.M(), this->_event->mGRefMult );
			_book->fill( prefix + "mass_runIndex", this->_event->mRunIndex, lv.M() );
			if ( "lsp_" == prefix || "lsn_" == prefix )
				_book->fill( "ls_mass_runIndex", this->_event->mRunIndex, lv.M() );
		}

		_book->fill( prefix + "mass_deltaPhi", fabs( lv1.DeltaPhi( lv2) ), lv.M() );
		_book->fill( prefix + "pT_deltaPhi", fabs( lv1.DeltaPhi( lv2) ), lv.Pt() );
		if ( "lsp_" == prefix || "lsn_" == prefix )
			_book->fill( "ls_deltaPhi", fabs( lv1.DeltaPhi( lv2) ) );

		_book->fill( prefix + "deltaPhi", fabs( lv1.DeltaPhi( lv2) ) );
		_book->fill( prefix + "deltaEta", fabs( lv1.PseudoRapidity() - lv2.PseudoRapidity() ) );

		_book->fill( prefix + "eta_phi", lv.Phi(), lv.PseudoRapidity() );
		_book->fill( prefix + "eta_leta", lv1.PseudoRapidity(), lv.PseudoRapidity() );
		_book->fill( prefix + "eta_sleta", lv2.PseudoRapidity(), lv.PseudoRapidity() );

	}

	void fillAggregates(  ){

		_book->fill( "uls_N", nUls );
		_book->fill( "lsp_N", nLsp );
		_book->fill( "lsn_N", nLsn );

		nLsn = nLsp = nUls = 0;
	}
	
};

#endif