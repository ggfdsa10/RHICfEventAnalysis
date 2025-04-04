#ifndef RHICfMassFitting_hh
#define RHICfMassFitting_hh

#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

#include "RHICfOptContainer.hh"
#include "RHICfTableMaker.hh"
#include "RHICfBinning.hh"

using namespace std;

static const int kMassFitParNum = 10;

class RHICfMassFitting
{
    struct MassFitResult
    {
        double allMean;
        double allSigma;
        double allPar[kMassFitParNum];
        vector<vector<double> > mean; // mass signal fit
        vector<vector<double> > sigma; // mass signal fit
        vector<vector<double> > massAllCounts; //intergral counts for mass All within 3 sigma
        vector<vector<double> > massSignalCounts; //intergral counts for mass Signal within 3 sigma
        vector<vector<double> > massBkgCounts;  //intergral counts for mass Bkg within 3 sigma
        vector<vector<double> > par[kMassFitParNum];

        void Resize(int ptNum, int xfNum)
        {
            allMean = 0.;
            allSigma = 0.;
            memset(allPar, 0., sizeof(allPar));
            mean.resize(ptNum, vector<double>(xfNum));
            sigma.resize(ptNum, vector<double>(xfNum));
            massAllCounts.resize(ptNum, vector<double>(xfNum));
            massSignalCounts.resize(ptNum, vector<double>(xfNum));
            massBkgCounts.resize(ptNum, vector<double>(xfNum));
            for(int i=0; i<kMassFitParNum; i++){par[i].resize(ptNum, vector<double>(xfNum));}
        };
    };

    public:
        RHICfMassFitting();
        ~RHICfMassFitting();

        void Init();
        void InitHist();

        void SetBinning(RHICfBinning* binning);

        void FillMass(int run, int type, int dle, double mass);
        void FillMass(int run, int type, int dle, int ptIdx, int xfIdx, double mass);

        void Fitting();
        void SaveMassData();

        int GetMassTableFlag();
        double GetMassMean(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);
        double GetMassSigma(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);
        double GetMassLowerBoundary(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);
        double GetMassUpperBoundary(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);

        double GetMassBkgLowerBoundary(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);
        double GetMassBkgUpperBoundary(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);

        double GetMassAllCounts(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);
        double GetMassSignalCounts(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);
        double GetMassBkgCounts(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);

        TH1D* GetMassHist(int run, int type, int dle, int ptIdx=-1, int xfIdx=-1);

    private:
        double* FitMassHist(TH1D* hist);
        double Pi0SignalFitter(double* x, double* par);
        double Pi0BkgFitter(double* x, double* par);
        double Pi0MassFitter(double* x, double* par);

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;

        TF1* mMassFitter;
        TF1* mSignalFitter;
        TF1* mBkgFitter;

        TH1D* mMassHistAll[kRunNum][kTypeNum][kDLENum]; // [run][type][dle]
        vector<vector<TH1D*> > mMassHistKinematic[kRunNum][kTypeNum][kDLENum]; // [run][type][dle][pt][xf]

        MassFitResult mMassFitResults[kRunNum][kTypeNum][kDLENum];
};

#endif
