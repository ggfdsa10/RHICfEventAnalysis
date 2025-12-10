#ifndef RHICfSimProcessRatio_hh
#define RHICfSimProcessRatio_hh

#include "TH1D.h"
#include "TMath.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"

using namespace std;

class RHICfSimProcessRatio
{
    public:
        RHICfSimProcessRatio(TString tableName="");
        ~RHICfSimProcessRatio();

        void Init();
        void InitSimProcessRatioData();
        void InitHist();

        void SetBinning(RHICfBinning* binning);

        void Calculate();
        void SaveSimProcessRatioData();

        int GetSimProcessRatioTableFlag();

    private:
        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;
        TString mTableName;

};

#endif
