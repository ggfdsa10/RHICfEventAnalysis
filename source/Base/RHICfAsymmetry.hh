#ifndef RHICfAsymmetry_hh
#define RHICfAsymmetry_hh

#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"
#include "RHICfDilutionFactor.hh"

using namespace std;

class RHICfAsymmetry
{
    public:
        RHICfAsymmetry();
        ~RHICfAsymmetry();

        void Init();
        void InitData();

        void SetBinning(RHICfBinning* binning);
        void SetDilution(RHICfDilutionFactor* dilution);

        void FillPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand);
        void FillBkgPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand);

        void Calculate();
        void SaveAsymmetryData();

    private:
        double RelativeLuminosity(int fillIdx);
        double BeamPolarization(int fillIdx);

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;
        RHICfDilutionFactor* mDilution;

        vector<vector<double> > mRightAsymmetry[kFillNum][kTypeNum][kDLENum];
        vector<vector<double> > mLeftAsymmetry[kFillNum][kTypeNum][kDLENum];
        vector<vector<double> > mBkgRightAsymmetry[kFillNum][kTypeNum][kDLENum];
        vector<vector<double> > mBkgLeftAsymmetry[kFillNum][kTypeNum][kDLENum];
};

#endif
