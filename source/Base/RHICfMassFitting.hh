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
        double allMassAllCounts;
        double allMassSignalCounts;
        double allMassBkgCounts;
        double allMassPeakDifference;
        double allMassFitChi2;
        double allPar[kMassFitParNum];
        vector<vector<double> > mean; // mass signal fit
        vector<vector<double> > sigma; // mass signal fit
        vector<vector<double> > massAllCounts; //intergral counts for mass All within 3 sigma
        vector<vector<double> > massSignalCounts; //intergral counts for mass Signal within 3 sigma
        vector<vector<double> > massBkgCounts;  //intergral counts for mass Bkg within 3 sigma
        vector<vector<double> > massPeakDifference;
        vector<vector<double> > massFitChi2;
        vector<vector<double> > par[kMassFitParNum];

        void Resize(int ptNum, int xfNum)
        {
            allMean = 0.;
            allSigma = 0.;
            allMassAllCounts = 0.;
            allMassSignalCounts = 0.;
            allMassBkgCounts = 0.;
            allMassPeakDifference = 0.;
            allMassFitChi2 = 0.;
            memset(allPar, 0., sizeof(allPar));
            mean.resize(ptNum, vector<double>(xfNum));
            sigma.resize(ptNum, vector<double>(xfNum));
            massAllCounts.resize(ptNum, vector<double>(xfNum));
            massSignalCounts.resize(ptNum, vector<double>(xfNum));
            massBkgCounts.resize(ptNum, vector<double>(xfNum));
            massPeakDifference.resize(ptNum, vector<double>(xfNum));
            massFitChi2.resize(ptNum, vector<double>(xfNum));
            for(int i=0; i<kMassFitParNum; i++){par[i].resize(ptNum, vector<double>(xfNum));}
        };
    };

    public:
        RHICfMassFitting();
        ~RHICfMassFitting();

        void Init();
        void InitMassData();
        void InitHist();

        void SetBinning(RHICfBinning* binning);

        void FillMass(int runIdx, int typeIdx, int dleIdx, double mass);
        void FillMass(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double mass);

        void Fitting();
        void SaveMassData();

        int GetMassTableFlag();
        double GetMassMean(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassSigma(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassLowerBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassUpperBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassBkgLowerBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassBkgUpperBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassAllCounts(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassSignalCounts(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassBkgCounts(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassPeakDifference(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassFitPeak(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassFitChi2(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        TH1D* GetMassHist(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

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
