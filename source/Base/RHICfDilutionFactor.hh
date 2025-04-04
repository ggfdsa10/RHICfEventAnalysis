#ifndef RHICfDilutionFactor_hh
#define RHICfDilutionFactor_hh

#include "TH1D.h"
#include "TMath.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"

using namespace std;

class RHICfDilutionFactor
{

    public:
        RHICfDilutionFactor();
        ~RHICfDilutionFactor();

        void Init();
        void InitHist();

        void SetBinning(RHICfBinning* binning);

        void FillAngle(int run, int type, int dle, int ptIdx, int xfIdx, double angle);

        void Calculate();
        void SaveDilutionData();

        int GetDilutionTableFlag();
        double GetDilutionFactor(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

        TH1D* GetDilutionHist(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

    private:
        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;

        vector<vector<TH1D*> > mDilutionHist[kRunNum][kTypeNum][kDLENum];
        vector<vector<double> > mDilutionFactor[kRunNum][kTypeNum][kDLENum];
};

#endif
