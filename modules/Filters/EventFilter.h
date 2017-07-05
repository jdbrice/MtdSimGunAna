#ifndef EVENT_FILTER_H
#define EVENT_FILTER_H

class EventFilter
{
public:
	EventFilter() {}
	EventFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~EventFilter() {}

	map<string, XmlRange> range;

	void load( XmlConfig &_cfg, string _nodePath ){
		
		vector<string> rngPaths=_cfg.childrenOf( _nodePath, "Range" );
		for ( auto p : rngPaths ){
			XmlRange rng( _cfg, p );

			range[ _cfg[ p + ":name" ] ] = rng;
			LOG_F( INFO, "[%s] = %s", _cfg[(p+":name")].c_str(), rng.toString().c_str() );
		}
	}

	bool pass( FemtoEvent *_event ){
		if ( range.count( "bin16" ) && !range[ "bin16" ].inInclusiveRange( _event->mBin16 ) )
			return false;
		return true;
	}

	bool fail( FemtoEvent *_event ){
		return !pass( _event );
	}



	
};


#endif