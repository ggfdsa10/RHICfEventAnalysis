#ifndef RHICfAsymmetry_hh
#define RHICfAsymmetry_hh

#include <vector> 

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"

#include "RHICfOptContainer.hh"
#include "RHICfEventDstReader.hh"
#include "RHICfTableMaker.hh"
#include "RHICfParticleMaker.hh"

#include "RHICfMassFitting.hh"
#include "RHICfBinning.hh"
#include "RHICfDilutionFactor.hh"
#include "RHICfPolarization.hh"

using namespace std;

static const int kBeamMetNum = 2;
static const int kBeamTOPRefNum = 2;

class RHICfAsymmetry : public RHICfOptContainer, RHICfTableMaker
{
    public:
        RHICfAsymmetry();
        ~RHICfAsymmetry();

        int Init();
        void Calculate();
    

    private:
        void AsymmetryPi0();
        void SystematicErrorCalculatePi0();
        
        void DrawAsymmetryGraph();
        TString GetSystemString();

        TGraphErrors* GetANGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx);
        TGraphErrors* GetANBkgGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx);
        TGraphErrors* GetANSubtGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx);

        void InitGraph();
        double RelativeLuminosity(int fillIdx);
        double BeamPolarization(int fillIdx);
        double BeamPolarizationError(int fillIdx);

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;

        RHICfBinning* mBinning[kRunNum][kBeamMetNum][kBeamTOPRefNum];
        RHICfMassFitting* mMassFitting[kRunNum][kBeamMetNum][kBeamTOPRefNum];
        RHICfDilutionFactor* mDilution[kRunNum][kBeamMetNum][kBeamTOPRefNum];
        RHICfPolarization* mPolarization[kRunNum][kBeamMetNum][kBeamTOPRefNum];
        
        TString mCondName;

        vector<TGraphErrors*> mGraphAN_pT[kRunNum][kTypeNum][kDLENum][5]; // [beamHit, beamScan, beam method summed]
        vector<TGraphErrors*> mGraphAN_xF[kRunNum][kTypeNum][kDLENum][5]; // [beamHit, beamScan, beam method summed]
        vector<TGraphErrors*> mGraphAN_Bkg_pT[kRunNum][kTypeNum][kDLENum][5]; // [beamHit, beamScan, beam method summed]
        vector<TGraphErrors*> mGraphAN_Bkg_xF[kRunNum][kTypeNum][kDLENum][5]; // [beamHit, beamScan, beam method summed]
        vector<TGraphErrors*> mGraphAN_Subt_pT[kRunNum][kTypeNum][kDLENum][5];
        vector<TGraphErrors*> mGraphAN_Subt_xF[kRunNum][kTypeNum][kDLENum][5];

        vector<TGraphErrors*> mGraphANSysErr_Calcul_subt_pT[kRunNum][kTypeNum][kDLENum][4][4][2]; // [Beam Method][dilution error, polarization error, B/S ratio error, DLE error][two point]
        vector<TGraphErrors*> mGraphANSysErr_Calcul_subt_xF[kRunNum][kTypeNum][kDLENum][4][4][2]; // [Beam Method][dilution error, polarization error, B/S ratio error, DLE error][two point]
        vector<TGraphErrors*> mGraphANSysErr_Subt_pT[kRunNum][kTypeNum][kDLENum][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]
        vector<TGraphErrors*> mGraphANSysErr_Subt_xF[kRunNum][kTypeNum][kDLENum][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]

        vector<TGraphErrors*> mGraphAN_pTSummary[kDLENum];
        vector<TGraphErrors*> mGraphAN_xFSummary[kDLENum];
        vector<TGraphErrors*> mGraphAN_pTSummarySysErr[kDLENum][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]
        vector<TGraphErrors*> mGraphAN_xFSummarySysErr[kDLENum][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]

        vector<TGraphErrors*> mGraphAN_pTSummary_merge[kDLENum];
        vector<TGraphErrors*> mGraphAN_xFSummary_merge[kDLENum];
        vector<TGraphErrors*> mGraphAN_pTSummarySysErr_merge[kDLENum][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]
        vector<TGraphErrors*> mGraphAN_xFSummarySysErr_merge[kDLENum][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]

        vector<TGraphErrors*> mGraphAN_pTSummary_DLEmerge[3];
        vector<TGraphErrors*> mGraphAN_xFSummary_DLEmerge[3];
        vector<TGraphErrors*> mGraphAN_pTSummarySysErr_DLEmerge[3][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]
        vector<TGraphErrors*> mGraphAN_xFSummarySysErr_DLEmerge[3][5]; // [dilution Err, polaarization err, B/S err, beam center err, total err]

};

#endif
