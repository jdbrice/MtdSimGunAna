#ifndef CORRECTED_SPECTRA_MAKER_H
#define CORRECTED_SPECTRA_MAKER_H

#include "HistoAnalyzer.h"
#include "XmlRange.h"

#define LOGURU_WITH_STREAMS 1
#include "vendor/loguru.h"

class CorrectedSpectraMaker : public HistoAnalyzer
{
protected:
	XmlRange pt_range;
	HistoBins mass_signal_bins;
public:
	CorrectedSpectraMaker() {}
	~CorrectedSpectraMaker() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();

		LOG_F( INFO, "initialize" );

		pt_range.loadConfig( config, "p.Range[0]" );
		LOG_F( INFO, "pt_range: %s", pt_range.toString().c_str()  );

		mass_signal_bins.load( config, "bins.mass" );
	}

	TH1 * geometric_mean( TH1* lsp, TH1* lsn, string name ){
		LOG_SCOPE_FUNCTION( 1 );
		TH1 * gm = (TH1*) lsp->Clone( name.c_str() );
		if ( nullptr == lsp || nullptr == lsn ) {
			LOG_F(ERROR, "NULL like sign histos");
			return nullptr;
		}

		int nBins = lsn->GetNbinsX();
		for ( int iBin = 1; iBin <= nBins; iBin++ ){
			
			double Npp = lsp->GetBinContent( iBin );
			double Nmm = lsn->GetBinContent( iBin );

			if ( 0 == Npp || 0 == Nmm ){
				gm->SetBinContent( iBin, 0 );
				gm->SetBinError( iBin, 0 );
				// LOG_F( ERROR, "Zero Bin in GM!" );
				continue;
			}

			double Epp = lsp->GetBinError( iBin ) / Npp;
			double Emm = lsn->GetBinError( iBin ) / Nmm;

			double N = 2 * sqrt( Npp * Nmm ); // DOES include the factor of 2!!!
			// Q = x^n  
			// ->
			// dQ/Q = n * dx/x
			double E = (Epp+Emm) * 0.5 * N;
			LOG_S( 1 ) << "B++ = " << Npp << " +/- " << Epp << endl;
			LOG_S( 1 ) << "B-- = " << Nmm << " +/- " << Emm << endl;
			LOG_S( 1 ) << "gm = " << N << " +/- " << E << endl;


			gm->SetBinContent( iBin, N );
			gm->SetBinError( iBin, E );
		}
		return gm;
	}


	TH1* project( string name, string ds, string nname){
		LOG_SCOPE_FUNCTION( INFO );
		LOG_F( INFO, "projecting %s (%s) from %f -> %f", name.c_str(), ds.c_str(), pt_range.min, pt_range.max );
		TH2 * h2 = get<TH2D>( name, ds );
		int b1 = h2->GetYaxis()->FindBin( pt_range.min );
		int b2 = h2->GetYaxis()->FindBin( pt_range.max );

		LOG_F( INFO, "(%f -> %f) = bins(%d -> %d)", pt_range.min, pt_range.max, b1, b2 );

		TH1* h1 = h2->ProjectionX( (nname + "_raw").c_str(), b1, b2 );
		h1 = (TH1*)h1->Rebin( mass_signal_bins.nBins(), nname.c_str(), mass_signal_bins.bins.data() );
		return h1;
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION(INFO);
		book->cd();
		TH1 * se_uls = project( "uls_pT_mass", "same_event", "se_uls" );
		TH1 * se_lsp = project( "lsp_pT_mass", "same_event", "se_lsp" );
		TH1 * se_lsn = project( "lsn_pT_mass", "same_event", "se_lsn" );

		TH1 * me_uls = project( "uls_pT_mass", "mixed_event", "me_uls" );
		TH1 * me_lsp = project( "lsp_pT_mass", "mixed_event", "me_lsp" );
		TH1 * me_lsn = project( "lsn_pT_mass", "mixed_event", "me_lsn" );

		TH1 * se_gm = geometric_mean( se_lsp, se_lsn, "se_gm" );
		se_gm->SetLineColor(kRed);

		TH1 * me_gm = geometric_mean( me_lsp, me_lsn, "me_gm" );
		me_gm->SetLineColor(kBlack);
		TH1 * acc_corr = (TH1*)me_uls->Clone("acc_corr");
		acc_corr->Divide( me_gm );

		TH1 * bg = (TH1*)se_gm->Clone("bg");
		bg->Multiply( acc_corr );


		TH1 * sig = (TH1*)se_uls->Clone( "sig" );
		sig->Add( bg, -1 );
	}


protected:
};


#endif