#ifndef FEMTO_DST_SKIMMER_H
#define FEMTO_DST_SKIMMER_H

#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMcTrack.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

// Analyzers
#include "FemtoDstSkimmer/PairHistogramMaker.h"
#include "FemtoDstSkimmer/TrackHistogramMaker.h"
#include "FemtoDstSkimmer/MtdHistogramMaker.h"

#include "Filters/TrackFilter.h"
#include "Filters/MtdFilter.h"
#include "Filters/EventFilter.h"

#define LOGURU_WITH_STREAMS 1
#include "vendor/loguru.h"



class PlcNode {
public:
	static vector<string> GEANTNAMES;
	PlcNode *genisis(){
		if ( nullptr == this->_parent )
			return this;
		return this->_parent->genisis();
	}
	PlcNode *_parent = nullptr;
	vector<PlcNode*> _children;
	vector<PlcNode*> siblings(){
		if ( this->_parent )
			return this->_parent->_children;
		return vector<PlcNode*>();
	}

	bool _matched = false;
	bool _mtdTruth = false;
	size_t _index = 0;


	// FemtoTrackProxy *_self;
	FemtoMcTrack * _self;

	void print_tree( int n = 0){

		string prefix = "";
		for ( int  i = 0; i < n; i++ ){
			prefix += "|";
			if ( i == n-1 )
				prefix += "--";
			else 
				prefix += "  ";
		}

		string matchStr = "";
		if ( this->_matched )
			matchStr = " (RC) ";
		string mtdStr = "";
		if ( this->_mtdTruth )
			mtdStr = " (MTD) ";
		// string k = "(" + dts( this->_self->mPt ) + ", " + dts( this->_self->mEta ) + dts( this->_self->mPhi ) + ")";
		string k="";
		// LOG_S( 2 )
		cerr << prefix << PlcNode::GEANTNAMES[ this->_self->mGeantPID ] << matchStr << mtdStr << std::setprecision(2) << "(" << this->_self->mPt << "," << this->_self->mEta << "," << this->_self->mPhi << ")" << endl;
		size_t nKids = this->_children.size();
		for ( size_t i = 0; i < nKids; i++ ){
			auto plc = this->_children[i];
			plc->print_tree( n+1 );
		}
	}
};

vector<string> PlcNode::GEANTNAMES;

class FemtoDstSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMcTrack> _rMcTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	vector<PlcNode*> plcTree;

public:
	virtual const char* classname() const {return "FemtoDstSkimmer";}
	FemtoDstSkimmer() {}
	~FemtoDstSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMcTracks.setup( chain, "McTracks" );
		_rHelices.setup( chain, "Helices" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		PlcNode::GEANTNAMES = config.getStringVector( "bins.GeantBins" );

	}

protected:

	void makeMcTree(){
		_event = _rEvent.get();

		size_t nMcTracks = _rMcTracks.N();
		plcTree.clear();

		// McTracks are already sorted such that a child never has a smaller index than a parent
		for (size_t i = 0; i < nMcTracks; i++){
			FemtoMcTrack *mcTrack = _rMcTracks.get( i );
			PlcNode* plc= new PlcNode();
			plc->_self = mcTrack;
			plc->_index = i;

			if ( plc->_self->mParentIndex >= 0 ){
				plc->_parent = plcTree[ plc->_self->mParentIndex ];
				plc->_parent->_children.push_back( plc );
			}

			plcTree.push_back( plc );
		}

		// if the McTracks were not ordered we would need this second loop
		
		// for (auto plc : plcTree ){
		// 	if ( plc->_self->mParentIndex >= 0 ){
		// 		plc->_parent = plcTree[ plc->_self->mParentIndex ];
		// 		plc->_parent->_children.push_back( plc );
		// 	}
		// }

		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rMtdPid, _rMcTracks );

			if ( _proxy._track->mMcIndex >= 0 )
				plcTree[ _proxy._track->mMcIndex ]->_matched = true;

			if ( _proxy._track->mMtdPidTraitsIndex >= 0 && _proxy._mtdPid){
				plcTree[_proxy._mtdPid->mIdTruth]->_mtdTruth = true;
			}
		}


		// for (auto plc : plcTree ){
		// 	if ( plc->_parent == nullptr ){
		// 		plc->print_tree();
		// 	}
		// }
	}


	int mc_charge( FemtoMcTrack *mcTrack ){
		if ( 8 == mcTrack->mGeantPID || 5 == mcTrack->mGeantPID )
			return 1;
		if ( 9 == mcTrack->mGeantPID || 6 == mcTrack->mGeantPID )
			return -1;
		return 0;
	}
	bool is_pion( FemtoMcTrack *mcTrack ){
		if( 8 == mcTrack->mGeantPID || 9 == mcTrack->mGeantPID )
			return true;
		return false;
	}
	bool is_muon( FemtoMcTrack *mcTrack ){
		if( 5 == mcTrack->mGeantPID || 6 == mcTrack->mGeantPID )
			return true;
		return false;
	}
	bool is_primary( PlcNode *plc ){
		if ( nullptr == plc->_parent )
			return true;
		return false;
	}

	// pass in the idTruth PlcNode
	bool punch_through( PlcNode * plc ){

		if ( is_muon( plc->_self ) )
			return false;

		return true;
	}

	string cs( int c ){
		if ( c > 0 )
			return "p";
		if ( c < 0 )
			return "m";

		return "0";
	}

	virtual void analyzeEvent(){
	
		_event = _rEvent.get();
		makeMcTree();
		
		size_t nTracks = _rTracks.N();
		size_t nMtdPid = _rMtdPid.N();
		FemtoTrackProxy _proxy;
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rMtdPid, _rMcTracks );

			if ( _proxy._mcTrack ){
				book->fill( "McIndex", _proxy._track->mMcIndex );
				book->fill( "PtResolution", (_proxy._track->mPt - _proxy._mcTrack->mPt ) / _proxy._mcTrack->mPt   );
				book->fill( "PtResolutionVsPt", _proxy._mcTrack->mPt, (_proxy._track->mPt - _proxy._mcTrack->mPt ) / _proxy._mcTrack->mPt   );
			}

			if ( _proxy._mcTrack && _proxy._mtdPid ){

				int bl = _proxy._mtdPid->backleg();
				int cell = _proxy._mtdPid->cell();
				int module = _proxy._mtdPid->module();

				if ( bl == 7 || bl == 23 ) {
					// LOG_F( INFO, "Skipping BL == 7, 23" );
					continue;
				}

				book->fill( "mtdMatchFlag", _proxy._mtdPid->mMatchFlag );
				// if ( _proxy._mtdPid->mMatchFlag != 7 ) {
				// 	LOG_F( INFO, "Skipping Match Flag 7" );
				// 	continue;
				// }

				// LOG_F( INFO, "mGeantPID = %d", _proxy._mcTrack->mGeantPID  );

				if ( 8 == _proxy._mcTrack->mGeantPID ){
					book->fill( "pip_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
					book->fill( "pip_mDeltaY_mBL", _proxy._track->mPt, _proxy._mtdPid->mDeltaY, bl );
					book->fill( "pip_mDeltaY_mMod", _proxy._track->mPt, _proxy._mtdPid->mDeltaY, module );
					book->fill( "pip_mDeltaY_mCell", _proxy._track->mPt, _proxy._mtdPid->mDeltaY, cell );
				}

				if ( 9 == _proxy._mcTrack->mGeantPID ){
					LOG_F( INFO, "charge = %d", _proxy._track->charge() );
					book->fill( "pim_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
					book->fill( "pim_mDeltaY_mBL", _proxy._track->mPt, _proxy._mtdPid->mDeltaY, bl );
					book->fill( "pim_mDeltaY_mMod", _proxy._track->mPt, _proxy._mtdPid->mDeltaY, module );
					book->fill( "pim_mDeltaY_mCell", _proxy._track->mPt, _proxy._mtdPid->mDeltaY, cell );
				}


				book->fill("idTruthVsIndex", _proxy._track->mMcIndex - (_proxy._mtdPid->mIdTruth) );

				PlcNode *plcMatch = plcTree[ _proxy._track->mMcIndex ];
				PlcNode *plcHit   = plcTree[ _proxy._mtdPid->mIdTruth ];

				if ( punch_through( plcHit ) ){
					book->fill( "punch_pi" + cs(_proxy._track->charge()) + "_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
				} else {
					book->fill( "sec_mu" + cs(_proxy._track->charge()) + "_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
				}

				// if ( _proxy._track->mMcIndex - (_proxy._mtdPid->mIdTruth-1) == 0 ){
				// 	if ( 8 == _proxy._mcTrack->mGeantPID )
				// 		book->fill( "tpip_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );

				// 	if ( 9 == _proxy._mcTrack->mGeantPID )
				// 		book->fill( "tpim_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
				// } else {
				// 	if ( 8 == _proxy._mcTrack->mGeantPID )
				// 		book->fill( "dpip_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );

				// 	if ( 9 == _proxy._mcTrack->mGeantPID )
				// 		book->fill( "dpim_mDeltaY", _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
				// }
			}

		}
	}


	virtual void postEventLoop(){
		TreeAnalyzer::postEventLoop();

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
		
	}
	
};

#endif
