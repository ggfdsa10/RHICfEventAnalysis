#ifndef RHICfMassFitting_hh
#define RHICfMassFitting_hh

#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TFile.h"
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
        double allMassAllCounts[2];    // integral counst, integral error
        double allMassSignalCounts[2]; // integral counst, integral error
        double allMassBkgCounts[2];    // integral counst, integral error
        double allMassFitChi2;
        double allPar[kMassFitParNum];
        double allParErr[kMassFitParNum];
        vector<vector<double> > mean; // mass signal fit
        vector<vector<double> > sigma; // mass signal fit
        vector<vector<double> > massAllCounts[2]; //intergral counts for mass All within 3 sigma and its error
        vector<vector<double> > massSignalCounts[2]; //intergral counts for mass Signal within 3 sigma and its error
        vector<vector<double> > massBkgCounts[2];  //intergral counts for mass Bkg within 3 sigma and its error
        vector<vector<double> > massFitChi2;
        vector<vector<double> > par[kMassFitParNum];
        vector<vector<double> > parErr[kMassFitParNum];

        void Resize(int ptNum, int xfNum)
        {
            allMean = 0.;
            allSigma = 0.;
            for(int i=0; i<2; i++){
                allMassAllCounts[i] = 0.;
                allMassSignalCounts[i] = 0.;
                allMassBkgCounts[i] = 0.;
            }
            allMassFitChi2 = 0.;
            memset(allPar, 0., sizeof(allPar));
            memset(allParErr, 0., sizeof(allParErr));
            mean.resize(ptNum, vector<double>(xfNum));
            sigma.resize(ptNum, vector<double>(xfNum));
            for(int i=0; i<2; i++){
                massAllCounts[i].resize(ptNum, vector<double>(xfNum));
                massSignalCounts[i].resize(ptNum, vector<double>(xfNum));
                massBkgCounts[i].resize(ptNum, vector<double>(xfNum));
            }
            massFitChi2.resize(ptNum, vector<double>(xfNum));
            for(int i=0; i<kMassFitParNum; i++){par[i].resize(ptNum, vector<double>(xfNum));}
            for(int i=0; i<kMassFitParNum; i++){parErr[i].resize(ptNum, vector<double>(xfNum));}
        };
    };

    public:
        RHICfMassFitting(TString tableName="");
        ~RHICfMassFitting();

        void Init();
        void InitMassData();
        void InitHist();

        void SetBinning(RHICfBinning* binning);

        void FillMass(int runIdx, int typeIdx, int dleIdx, double mass);
        void FillMass(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double mass);

        void Fitting();
        void SaveMassData();

        bool IsExistMassHistFile();
        int GetMassTableFlag();
        double GetMassMean(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassSigma(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassLowerBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassUpperBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassBkgLowerBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassBkgUpperBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassAllCounts(bool isErr, int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassSignalCounts(bool isErr, int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetMassBkgCounts(bool isErr, int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetBSRatio(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double GetBSRatioErr(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        double GetMassFitChi2(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);
        double* GetMassFitPar(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        TH1D* GetMassHist(int runIdx, int typeIdx, int dleIdx, int ptIdx=-1, int xfIdx=-1);

        // new parameter getting  functions for Gaussian Process Regression method 
        double GetBSRatio_GPR(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetBSRatioUpper_GPR(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);
        double GetBSRatioLower_GPR(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx);

    private:
        bool FitMassHist(TH1D* hist);
        void FindFitRange(TH1D* hist, double& lower, double& upper);

        void InitMassHistFile();
        
        double Pi0SignalFitter(double* x, double* par);
        double Pi0BkgFitter(double* x, double* par);
        double Pi0MassFitter(double* x, double* par);

        RHICfOptContainer* mOptContainer;
        RHICfTableMaker* mTableMaker;
        RHICfBinning* mBinning;
        TString mTableName;

        int mSigFitParNum;
        int mBkgFitParNum;
        TF1* mMassFitter;
        TF1* mSignalFitter_drawing;
        TF1* mBkgFitter_drawing;

        TF1* mSignalFitter;
        TF1* mBkgFitter;     
        TFitResultPtr mFitResult;   

        TH1D* mMassHistAll[kRunNum][kTypeNum][kDLENum]; // [run][type][dle]
        vector<vector<TH1D*> > mMassHistKinematic[kRunNum][kTypeNum][kDLENum]; // [run][type][dle][pt][xf]

        MassFitResult mMassFitResults[kRunNum][kTypeNum][kDLENum];

        bool mIsFindMassFile;
        TFile* mMassFile;
        
};

#endif
