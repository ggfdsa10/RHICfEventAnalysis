#ifndef RHICfSimAnalysis_hh
#define RHICfSimAnalysis_hh

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#include "RHICfOptContainer.hh"
#include "RHICfSimDstReader.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"
#include "RHICfParticleMaker.hh"
#include "StRHICfMiniSimDst.h"

using namespace std;

class RHICfSimAnalysis : public RHICfOptContainer
{
    public:
        RHICfSimAnalysis();
        ~RHICfSimAnalysis();

        int Init();
        int Calculate();
        int Finish();

    private:
        int GetTruthParFlag(int type);

        RHICfBinning* mBinning;

        RHICfSimDstReader* mSimDstReader[rGeneratorNum][kRunNum];
        StRHICfMiniSimDst* mMiniSimDst;

        int mCalRunType;
        const int mModelNum = 4;

};

#endif
