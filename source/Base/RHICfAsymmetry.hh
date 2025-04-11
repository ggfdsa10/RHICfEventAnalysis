#ifndef RHICfAsymmetry_hh
#define RHICfAsymmetry_hh

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TCanvas.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"
#include "RHICfMassFitting.hh"
#include "RHICfDilutionFactor.hh"

using namespace std;

class RHICfAsymmetry
{
    public:
        RHICfAsymmetry();
        ~RHICfAsymmetry();

        void Init();
        void InitAsymmetryData();

        void SetBinning(RHICfBinning* binning);
        void SetMassFitting(RHICfMassFitting* massFitting);
        void SetDilution(RHICfDilutionFactor* dilution);

        void FillPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double angle, bool isSpinUp);
        void FillBkgPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double angle, bool isSpinUp);

        void CalculateAN();
        void SaveAsymmetryData();

        int GetAsymmetryTableFlag();
        double GetPolNum(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand);
        double GetBkgPolNum(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand);

        TGraphErrors* GetANGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx);
        TGraphErrors* GetANGraph(bool isPtGraph, int anType, int dleIdx, int binIdx);
        TGraphErrors* GetANSummaryGraph(bool isPtGraph, int anType, int dleIdx);

    private:
        void AsymmetryPi0();
        void AsymmetryNeutron();

        double GetInclusiveAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetInclusiveANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetInclusiveANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

        double GetBkgAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetBkgANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetBkgANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

        double GetExclusiveAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetExclusiveANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetExclusiveANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

        bool RejectPoint(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

        void InitGraph();
        double RelativeLuminosity(int fillIdx);
        double BeamPolarization(int fillIdx);

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;
        RHICfMassFitting* mMassFitting;
        RHICfDilutionFactor* mDilution;

        vector<vector<double> > mRightAsymmetry[kFillNum][kTypeNum][kDLENum];
        vector<vector<double> > mLeftAsymmetry[kFillNum][kTypeNum][kDLENum];
        vector<vector<double> > mBkgRightAsymmetry[kFillNum][kTypeNum][kDLENum];
        vector<vector<double> > mBkgLeftAsymmetry[kFillNum][kTypeNum][kDLENum];

        vector<TGraphErrors*> mGraphAN_Global_pT[kRunNum][kTypeNum][kDLENum][3]; // [inclusive, bkg, exclusive]
        vector<TGraphErrors*> mGraphAN_Global_xF[kRunNum][kTypeNum][kDLENum][3]; // [inclusive, bkg, exclusive]
        vector<TGraphErrors*> mGraphAN_pT[kDLENum][3]; // [inclusive, bkg, exclusive]
        vector<TGraphErrors*> mGraphAN_xF[kDLENum][3]; // [inclusive, bkg, exclusive]
        TGraphErrors* mGraphAN_pTSummary[kDLENum][3];
        TGraphErrors* mGraphAN_xFSummary[kDLENum][3];
};

#endif
