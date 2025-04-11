#ifndef RHICfSystematicError_hh
#define RHICfSystematicError_hh

#include "TGraphErrors.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"

class RHICfSystematicError
{
    public:
        RHICfSystematicError();
        ~RHICfSystematicError();

        void Init();
        void InitDLESystematicData();

        void SetBinning(RHICfBinning* binning);

        void CalculateDLESystematic();
        void SaveDLESystematicData();

        int GetDLESystematicTableFlag();

        void FillPolSDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, bool blueSpinUp);
        void FillPolDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int bbcSE, int bbcLE, bool blueSpinUp);
        void FillPolNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int btof, int bbcSE, int bbcLE, int bbcSW, int bbcLW, bool blueSpinUp);
        
        double GetPolSDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int polIdx);
        double GetPolDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ddleIdx, int polIdx);
        double GetPolNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ndleIdx1, int ndleIdx2, int ndleIdx3, int polIdx);

    private:
        void InitGraph();
        double GetRawAsymmetrySDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx);
        double GetRawAsymmetryDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ddleIdx);
        double GetRawAsymmetryNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ndleIdx1, int ndleIdx2, int ndleIdx3);

        double RelativeLuminosity(int fillIdx);

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;

        static const int mDDLEBins = 3;
        static const int mNDLEBins = 3;
        int mBTofMult_thrBoundary[4];
        int mBBCSmallEast_thrBoundary[4];
        int mBBCSmallWest_thrBoundary[4];
        int mBBCLargeEast_thrBoundary[4];
        int mBBCLargeWest_thrBoundary[4];

        vector<vector<double> > mBeamSpinSDLE[2][kFillNum][kTypeNum];
        vector<vector<double> > mBeamSpinDDLE[2][mDDLEBins][kFillNum][kTypeNum];
        vector<vector<double> > mBeamSpinNDLE[2][mNDLEBins][mNDLEBins][mNDLEBins][kFillNum][kTypeNum];

        vector<TGraphErrors*> mGraphAsymSDLEPt[kFillNum][kTypeNum];
        vector<TGraphErrors*> mGraphAsymSDLEXf[kFillNum][kTypeNum];
        vector<TGraphErrors*> mGraphAsymDDLEPt[mDDLEBins][kFillNum][kTypeNum];
        vector<TGraphErrors*> mGraphAsymDDLEXf[mDDLEBins][kFillNum][kTypeNum];
        vector<TGraphErrors*> mGraphAsymNDLEPt[mNDLEBins][mNDLEBins][mNDLEBins][kFillNum][kTypeNum];
        vector<TGraphErrors*> mGraphAsymNDLEXf[mNDLEBins][mNDLEBins][mNDLEBins][kFillNum][kTypeNum];

};

#endif
