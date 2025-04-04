#ifndef RHICfFigureDrawing_hh
#define RHICfFigureDrawing_hh

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"

#include "RHICfBinning.hh"
#include "RHICfMassFitting.hh"
#include "RHICfDilutionFactor.hh"

using namespace std;

class RHICfFigureDrawing
{
    public:
        RHICfFigureDrawing();
        ~RHICfFigureDrawing();

        void Init();
        void SetBinning(RHICfBinning* binning);
        void SetMassFitting(RHICfMassFitting* fitting);
        void SetDilution(RHICfDilutionFactor* dilution);

        void DrawMassHist();
        void DrawBinningHist();
        void DrawDilutionHist();

    private:
        TString GetSystemString();

        RHICfOptContainer* mOptContainer;
        RHICfBinning* mBinning;
        RHICfMassFitting* mFitting;
        RHICfDilutionFactor* mDilution;

};

#endif
