#ifndef FEMTO_DST_SKIMMER_H
#define FEMTO_DST_SKIMMER_H

#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMcTrack.h"
#include "FemtoDstFormat/FemtoMcVertex.h"
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


class FemtoTracklet {
	public:
		FemtoTracklet( ) {
			mPt = 0;
			mEta = 0;
			mPhi = 0;
		}

		FemtoTracklet( FemtoTrack * t ){
			mPt  = t->mPt;
			mEta = t->mEta;
			mPhi = t->mPhi;
		}
		
		FemtoTracklet( FemtoMcTrack * t ){
			mPt  = t->mPt;
			mEta = t->mEta;
			mPhi = t->mPhi;
		}

		void set( FemtoMcTrack * t ){
			mPt  = t->mPt;
			mEta = t->mEta;
			mPhi = t->mPhi;
		}
		void set( FemtoTrack * t ){
			mPt  = t->mPt;
			mEta = t->mEta;
			mPhi = t->mPhi;
		}



	Float_t   mPt;        // primary track px
	Float_t   mEta;       // primary track py
	Float_t   mPhi;       // primary track pz
};

vector<string> PlcNode::GEANTNAMES;

class FemtoDstSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMcTrack> _rMcTracks;
	TClonesArrayReader<FemtoMcVertex> _rMcVertices;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	vector<PlcNode*> plcTree;

	bool isRealData = false;

public:
	virtual const char* classname() const {return "FemtoDstSkimmer";}
	FemtoDstSkimmer() {}
	~FemtoDstSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMcTracks.setup( chain, "McTracks" );
		_rMcVertices.setup( chain, "McVertices" );
		_rHelices.setup( chain, "Helices" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		PlcNode::GEANTNAMES = config.getStringVector( "bins.GeantBins" );

	}

protected:

	void makeMcTree(){
		_event = _rEvent.get();

		size_t nMcTracks = _rMcTracks.N();
		plcTree.clear();

		// Generate the plcNodes linked list
		// NB. McTracks are already sorted such that a child never has a smaller index than a parent
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

		// Now search for the match states
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

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();

		// for ( int iBL = 0; iBL < 30; iBL++ ){
		// 	for ( int iCell = 0; iCell < 12; iCell++ ){
		// 		book->clone( "pos_mDeltaY_mPt", "pos_mDeltaY_mPt_mBL"+ts(iBL) +"_mCell" + ts(iCell) );
		// 	}
		// }

	}

	virtual void fillHistograms( FemtoTrackProxy &_proxy ){
		if ( _proxy._track == nullptr ) return;
		string prefix = "pos_";
		if ( _proxy._track->charge() < 0 )
			prefix = "neg_";

		FemtoTracklet tl;
		if ( nullptr != _proxy._mcTrack )
			tl.set( _proxy._mcTrack );
		else 
			tl.set( _proxy._track );

		if ( tl.mPt <= 0.1 ){
			// LOG_F( INFO, "what mPt=%f", _proxy._track->mPt );
			return;
		}

		if ( _proxy._mcTrack ){
			
			FemtoMcVertex *stopVertex = nullptr;
			if ( _proxy._mcTrack->mStopVertexIndex >= 0 )
				stopVertex = _rMcVertices.get( _proxy._mcTrack->mStopVertexIndex );
			FemtoMcVertex *startVertex = nullptr;
			if ( _proxy._mcTrack->mStartVertexIndex >= 0 )
				startVertex = _rMcVertices.get( _proxy._mcTrack->mStartVertexIndex );
			
			PlcNode *pNode = plcTree[ _proxy._track->mMcIndex ];
			if ( pNode->_parent == nullptr && pNode->_children.size() == 2 && ( pNode->_children[0]->_self->mGeantPID==5 || pNode->_children[0]->_self->mGeantPID==6 || pNode->_children[1]->_self->mGeantPID==5 || pNode->_children[1]->_self->mGeantPID==6 ) ){

			}
		}

		// all tracks
		book->fill( prefix + "TPC_mEta_mPt", tl.mPt, tl.mEta );
		book->fill( prefix + "TPC_mPhi_mPt", tl.mPt, tl.mPhi );
		book->fill( prefix + "TPC_mEta_mPhi", tl.mPhi, tl.mEta );

		book->fill( prefix + "mNHitsFit_mPt", tl.mPt, _proxy._track->mNHitsFit );

		if ( fabs(_proxy._track->mNHitsFit) < 15 ) return;

		if ( _proxy._mtdPid ){

			int bl     = _proxy._mtdPid->backleg() + 1;
			int cell   = _proxy._mtdPid->cell() + 1;
			int strip  = _proxy._mtdPid->mMtdHitChan - bl * 60 + 1;
			int module = _proxy._mtdPid->module() + 1;

			if ( bl == 7 || bl == 23 ) {
				// LOG_F( INFO, "Skipping BL == 7, 23" );
				return;
			}

			book->fill( prefix + "mMtdMatchFlag_mPt", tl.mPt, _proxy._mtdPid->mMatchFlag );

			if ( _proxy._mtdPid->mMatchFlag > 1 ) return;

			book->fill( prefix + "mEta_mPt", tl.mPt, tl.mEta );
			book->fill( prefix + "mPhi_mPt", tl.mPt, tl.mPhi );
			book->fill( prefix + "mEta_mPhi", tl.mPhi, tl.mEta );
			book->fill( prefix + "mBL_mPhi", tl.mPhi, bl );

			
			book->fill( prefix + "MtdHitMap", bl, strip );

			book->fill( prefix + "mDeltaY_mPt", tl.mPt, _proxy._mtdPid->mDeltaY );
			book->fill( prefix + "mDeltaZ_mPt", tl.mPt, _proxy._mtdPid->mDeltaZ );
			book->fill( prefix + "mCell_mPt", tl.mPt, cell );

			// book->fill( prefix + "mDeltaY_mPt_mStrip" + ts(strip), tl.mPt, _proxy._mtdPid->mDeltaY );
			book->fill( prefix + "mDeltaY_mPt_mCell" + ts(cell), tl.mPt, _proxy._mtdPid->mDeltaY );

			// book->fill( prefix + "mDeltaY_mPt_mBL" +ts(bl) + "_mCell" + ts(cell), tl.mPt, _proxy._mtdPid->mDeltaY );

			book->fill( prefix + "mDeltaY_mBL_mCell" + ts(cell), bl, _proxy._mtdPid->mDeltaY );
			book->fill( prefix + "mDeltaY_mBL_mStrip" + ts(strip), bl, _proxy._mtdPid->mDeltaY );

		
		}

	}

	virtual void analyzeEvent(){
		_event = _rEvent.get();
		
		size_t nMcTracks = _rMcTracks.N();
		if ( nMcTracks > 0 )
			makeMcTree();
		else 
			isRealData = 0;
		
		size_t nTracks = _rTracks.N();
		size_t nMtdPid = _rMtdPid.N();
		FemtoTrackProxy _proxy;
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rMtdPid, _rMcTracks );
			fillHistograms( _proxy );
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
