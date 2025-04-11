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
#include "RHICfAsymmetry.hh"

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
        void SetAsymmetry(RHICfAsymmetry* asymmetry);

        void DrawMassHist();
        void DrawBinningHist();
        void DrawDilutionHist();
        void DrawAsymmetryGraph();

    private:
        TString GetSystemString();
        // int GetMarkerStyleIdx(int idx);
        // int GetMarkerStyle(int idx);

        RHICfOptContainer* mOptContainer;
        RHICfBinning* mBinning;
        RHICfMassFitting* mFitting;
        RHICfDilutionFactor* mDilution;
        RHICfAsymmetry* mAsymmetry;

};

#endif
