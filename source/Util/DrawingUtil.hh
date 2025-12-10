#ifndef DrawingUtil_hh
#define DrawingUtil_hh

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "TSystem.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"

using namespace std;

class DrawingUtil
{
    struct InfoBase
    {
        TString title;
        int cIdx;
        int objIdx;
        int colorIdx;
        int markerIdx;

        void Clear(){
            title = "";
            cIdx = 0;
            objIdx = 0;
            colorIdx = -1;
            markerIdx = -1;
        };
    };

    struct InfoHist
    {
        TH1D* hist;
        InfoBase base;
    };

    struct InfoGraph
    {
        TGraphErrors* graph;
        InfoBase base;
    };

    struct InfoText
    {
        TString text;
        bool persistency;
        int cIdx;
        int objIdx;
        double x;
        double y;
        double size;
        double font;
        TObject* obj;
        TString opt;

        void Clear(){
            text = "";
            persistency = false;
            cIdx = 0;
            objIdx = 0;
            x = -1;
            y = -1;
            size = -1;
            font = -1;
            obj = 0;
            opt = "";
        };
    };
    
    public:
        DrawingUtil(TString figurePath="./");
        ~DrawingUtil();

        void Init();
        void Clear(TString opt="");
        void DrawHist(TString opt="");
        void DrawHistWithRatio(int baseIdx, TString opt="", TString yTitle="");
        void DrawGraph(TString opt="");

        void SaveFigure(TString name);

        void SetCanvas(int divideNum, TString name="", double size=700.);
        void SetXLabelName(int i, TString name);

        void SetPalette(int idx){mColorPaletteIdx = idx;}
        void SetMarkerFill(int idx){mMarkerStyleIdx = idx;}

        void SetTH1DRatio(int num, double* y, double base, TString title="", int cIdx=-1, int colorIdx=-1, int markerIdx=-1);
        void SetTH1D(TH1D* hist, TString title="", int cIdx=-1, int colorIdx=-1, int markerIdx=-1);
        TH1D* GetTH1D(int cIdx, int objIdx);

        void SetGraph(int num, double* x, double* y, TString title="", int cIdx=-1, int colorIdx=-1, int markerIdx=-1);
        void SetGraphRatio(int num, double* x, double* y, double* yBase, TString title="", int cIdx=-1, int colorIdx=-1, int markerIdx=-1);
        void SetGraphError(int num, double* x, double* y, double* ex, double* ey, TString title="", int cIdx=-1, int colorIdx=-1, int markerIdx=-1);
        TGraphErrors* SetGraph(TString title="", int cIdx=-1, int colorIdx=-1, int markerIdx=-1);
        TGraphErrors* GetGraph(int cIdx, int objIdx);

        void SetText(bool isPersistant, TString text, int cIdx=-1, double x=-1, double y=-1, double size=-1, int font=-1);
        void OffPersistantText(int c, int idx);
        void SetLegend(bool isPersistant, TObject* obj, TString text, int cIdx=-1, TString opt="");

        TCanvas* GetCanvas();
        TLatex* GetLatex();
        TLegend* GetLegend();
        TLine* GetLine();

    private:
        int GetColor(int idx);
        int GetMarker(int idx);

        void InitHist(int& cNum);
        void InitGraph(int& cNum);

        void DrawText(int cNum);
        void DrawLegend(int cNum);

        // ************ Drawing Options ******************
        // n, norm : hiostogram normalization
        // xlog, ylog, zlog : log drawing 
        // ratio : ratio 0 to 1 
        // **********************************************
        bool InitDrawOption(TString opt);
        bool GetDrawFlag(TString opt);
        TString GetDrawValue(TString opt);
        
        bool GetTitle(TString inputTitle, TString& name, TString& xTitle, TString& yTitle);
        void NormalizedTH1D(TH1D* hist, double scale=-1.);

        TCanvas* mCanvas = nullptr;
        vector<TString> mDrawOptArr;

        vector<InfoHist> mHist;
        vector<InfoGraph> mGraph;

        vector<InfoText> mLatex;
        vector<InfoText> mLegend;

        int mColorPaletteIdx;
        int mMarkerStyleIdx;
        vector<pair<int, TString> > mLabelName;

        TString mFigurePath;
};

#endif
