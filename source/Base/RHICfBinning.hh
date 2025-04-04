#ifndef RHICfBinning_hh
#define RHICfBinning_hh

#include "TH1D.h"
#include "TH2D.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"

class RHICfBinning
{
    public:
        RHICfBinning();
        ~RHICfBinning();

        void Init();

        void FillKinematics(int runIdx, int typeIdx, int dleIdx, double pt, double xf);

        void Binning();
        void SaveBinningData();

        int GetBinningTableFlag();
        int GetPtBinNum(int runIdx, int typeIdx, int dleIdx);
        int GetXfBinNum(int runIdx, int typeIdx, int dleIdx);
        double GetPtBinBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx);
        double GetXfBinBoundary(int runIdx, int typeIdx, int dleIdx, int xfIdx);

        double GetPtBinMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetXfBinMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

        int GetGlobalPtBinNum();
        int GetGlobalXfBinNum();
        double GetGlobalPtBinBoundary(int ptIdx);
        double GetGlobalXfBinBoundary(int xfIdx);

        TH1D* GetPtHist(int runIdx, int typeIdx, int dleIdx);
        TH1D* GetXfHist(int runIdx, int typeIdx, int dleIdx);
        TH2D* GetKinematicsHist(int runIdx, int typeIdx, int dleIdx);

    private:
        void Pi0GlobalBinning();

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;

        TH1D* mPtHist[kRunNum][kTypeNum][kDLENum]; 
        TH1D* mXfHist[kRunNum][kTypeNum][kDLENum]; 
        TH2D* mPtXfBinHist[kRunNum][kTypeNum][kDLENum]; 
        TH2D* mKinematicsHist[kRunNum][kTypeNum][kDLENum];

        vector<double> mGlobalPtBoundary;
        vector<double> mGlobalXfBoundary;
        vector<double> mPtBoundary[kRunNum][kTypeNum][kDLENum];
        vector<double> mXfBoundary[kRunNum][kTypeNum][kDLENum];
        vector<vector<double> > mPtMean[kRunNum][kTypeNum][kDLENum];
        vector<vector<double> > mXfMean[kRunNum][kTypeNum][kDLENum];

        int mLocalBinNum;
        int mPtBins;
        int mXfBins;

};

#endif
