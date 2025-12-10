#include "RHICfMassFitting.hh"

#include <TVirtualFitter.h>
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TObjArray.h"
#include "TObjString.h"

RHICfMassFitting::RHICfMassFitting(TString tableName) : mTableName(tableName) 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
    mMassFitter = 0;
    mSignalFitter_drawing = 0;
    mBkgFitter_drawing = 0;
    mSignalFitter = 0;
    mBkgFitter = 0;

    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                mMassHistAll[run][type][dle] = 0;
                mMassHistKinematic[run][type][dle].clear();
                mMassFitResults[run][type][dle].Resize(0, 0);
            }
        }
    }
}

RHICfMassFitting::~RHICfMassFitting()
{
}

void RHICfMassFitting::Init()
{
    // mass fitter
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        double fitLowerBoundary = 20.;
        double fitUpperBoundary = 200.;
        mSigFitParNum = 3;
        mBkgFitParNum = 7;

        if(!mMassFitter){
            mMassFitter = new TF1("MassFitter", this, &RHICfMassFitting::Pi0MassFitter, fitLowerBoundary, fitUpperBoundary, kMassFitParNum);
            mMassFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
            mMassFitter -> SetLineColor(kBlack);
            mMassFitter -> SetLineWidth(1.5);
        }
        if(!mSignalFitter_drawing){
            mSignalFitter_drawing = new TF1("SignalFitter_drawing", this, &RHICfMassFitting::Pi0SignalFitter, fitLowerBoundary, fitUpperBoundary, mSigFitParNum);
            mSignalFitter_drawing -> SetParameters(0, 0, 0);
            mSignalFitter_drawing -> SetLineColor(kRed);
            mSignalFitter_drawing -> SetLineWidth(1.5);
        }
        if(!mSignalFitter){
            mSignalFitter = new TF1("SignalFitter", this, &RHICfMassFitting::Pi0SignalFitter, fitLowerBoundary, fitUpperBoundary, mSigFitParNum);
            mSignalFitter -> SetParameters(0, 0, 0);
        }
        if(!mBkgFitter_drawing){
            mBkgFitter_drawing = new TF1("BkgFitter_drawing", this, &RHICfMassFitting::Pi0BkgFitter, fitLowerBoundary, fitUpperBoundary, mBkgFitParNum);
            mBkgFitter_drawing -> SetParameters(0, 0, 0, 0, 0, 0, 0);
            mBkgFitter_drawing -> SetLineColor(kBlue);
            mBkgFitter_drawing -> SetLineWidth(1.5);
        }        
        if(!mBkgFitter){
            mBkgFitter = new TF1("BkgFitter", this, &RHICfMassFitting::Pi0BkgFitter, fitLowerBoundary, fitUpperBoundary, mBkgFitParNum);
            mBkgFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0);
        }        
    }

    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        // to be updated
    }   

    if(mTableName.Sizeof() == 1){mTableName = mOptContainer->GetTableSubName();}
    mTableMaker -> InitTable("MassFit"+mTableName);

    InitMassHistFile();

    cout << "RHICfMassFitting::Init() -- Done." << endl;
}

void RHICfMassFitting::InitMassData()
{
    int flag = GetMassTableFlag();
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptBinNum = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, 1, -1, -1)-1;
                int xfBinNum = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, -1, 1, -1)-1;
                if(ptBinNum > 0 && xfBinNum > 0){
                    mMassFitResults[run][type][dle].Resize(ptBinNum, xfBinNum);
                }

                if(flag == kExistTable){
                    mMassFitResults[run][type][dle].allMean = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 0);
                    mMassFitResults[run][type][dle].allSigma = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 1);
                    mMassFitResults[run][type][dle].allMassAllCounts[0] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 2);
                    mMassFitResults[run][type][dle].allMassAllCounts[1] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 3);
                    mMassFitResults[run][type][dle].allMassSignalCounts[0] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 4);
                    mMassFitResults[run][type][dle].allMassSignalCounts[1] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 5);
                    mMassFitResults[run][type][dle].allMassBkgCounts[0] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 6);
                    mMassFitResults[run][type][dle].allMassBkgCounts[1] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 7);
                    mMassFitResults[run][type][dle].allMassFitChi2 = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 8);

                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 9+i);
                    }

                    for(int pt=0; pt<ptBinNum; pt++){
                        for(int xf=0; xf<xfBinNum; xf++){
                            mMassFitResults[run][type][dle].mean[pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 0);
                            mMassFitResults[run][type][dle].sigma[pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 1);
                            mMassFitResults[run][type][dle].massAllCounts[0][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 2);
                            mMassFitResults[run][type][dle].massAllCounts[1][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 3);
                            mMassFitResults[run][type][dle].massSignalCounts[0][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 4);
                            mMassFitResults[run][type][dle].massSignalCounts[1][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 5);
                            mMassFitResults[run][type][dle].massBkgCounts[0][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 6);
                            mMassFitResults[run][type][dle].massBkgCounts[1][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 7);
                            mMassFitResults[run][type][dle].massFitChi2[pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 8);

                            for(int i=0; i<kMassFitParNum; i++){
                                mMassFitResults[run][type][dle].par[i][pt][xf] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, pt, xf, 9+i);
                            }
                        }
                    }
                }
                else if(flag == kExistPartOfTable){
                    mMassFitResults[run][type][dle].allMean = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 0);
                    mMassFitResults[run][type][dle].allSigma = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 1);
                    mMassFitResults[run][type][dle].allMassAllCounts[0] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 2);
                    mMassFitResults[run][type][dle].allMassAllCounts[1] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 3);
                    mMassFitResults[run][type][dle].allMassSignalCounts[0] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 4);
                    mMassFitResults[run][type][dle].allMassSignalCounts[1] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 5);
                    mMassFitResults[run][type][dle].allMassBkgCounts[0] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 6);
                    mMassFitResults[run][type][dle].allMassBkgCounts[1] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 7);
                    mMassFitResults[run][type][dle].allMassFitChi2 = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 8);
                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = mTableMaker->GetTableData("MassFit"+mTableName, run, type, dle, -1, -1, 9+i);
                    }

                }
            }
        }
    }
}

void RHICfMassFitting::InitHist()
{
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        // mass histogram
        int massBins;
        double massLowerBoundary;
        double massUpperBoundary;

        for(int run=0; run<kRunNum; run++){
            TString runName = mOptContainer -> GetRunTypeName(run);
            for(int type=0; type<kTypeNum; type++){
                for(int dle=0; dle<kDLENum; dle++){
                    TString dleName = mOptContainer -> GetDLEName(dle);
                    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
                        massBins = 75;
                        massLowerBoundary = 0.;
                        massUpperBoundary = 300.;

                        if(!mMassHistAll[run][type][dle] && !mIsFindMassFile){
                            mMassHistAll[run][type][dle] = new TH1D(Form("pi0MassHistAll_%s_type%i_%s", runName.Data(), type, dleName.Data()), "", massBins, massLowerBoundary, massUpperBoundary);
                            mMassHistAll[run][type][dle] -> SetStats(0);
                            mMassHistAll[run][type][dle] -> SetTitle(Form("%s Invariant Mass; M_{#gamma#gamma} [MeV/c^{2}]; Counts", (mOptContainer->GetParticleRunName()).Data()));
                        }
                        mMassHistAll[run][type][dle] -> Clear("ICESM");
                        mMassHistAll[run][type][dle] -> Sumw2();
                        
                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);

                        if(ptNum > 0 && xfNum > 0){
                            if(mMassHistKinematic[run][type][dle].size() == 0){
                                mMassHistKinematic[run][type][dle].resize(ptNum, vector<TH1D*>(xfNum));
                            }
                            for(int pt=0; pt<ptNum; pt++){
                                for(int xf=0; xf<xfNum; xf++){
                                    if(!mMassHistKinematic[run][type][dle][pt][xf] && !mIsFindMassFile){
                                        mMassHistKinematic[run][type][dle][pt][xf] = new TH1D(Form("pi0MassHist_%s_type%i_%s_pt%i_xf%i", runName.Data(), type, dleName.Data(), pt, xf), "", massBins, massLowerBoundary, massUpperBoundary);
                                        mMassHistKinematic[run][type][dle][pt][xf] -> SetStats(0);
                                        mMassHistKinematic[run][type][dle][pt][xf] -> SetTitle(Form("%s Invariant Mass; M_{#gamma#gamma} [MeV/c^{2}]; Counts", (mOptContainer->GetParticleRunName()).Data()));
                                        mMassHistKinematic[run][type][dle][pt][xf] -> Sumw2();
                                    }
                                }
                            }
                            if(mMassFitResults[run][type][dle].mean.size() == 0){
                                mMassFitResults[run][type][dle].Resize(ptNum, xfNum);  
                            }
                        }
                    }
                }
            }
        }
    }
}

void RHICfMassFitting::SetBinning(RHICfBinning* binning)
{
    mBinning = binning;
}

void RHICfMassFitting::FillMass(int runIdx, int typeIdx, int dleIdx, double mass)
{
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        if(mMassHistAll[runIdx][typeIdx][dleIdx]){
            mMassHistAll[runIdx][typeIdx][dleIdx] -> Fill(mass);
        }
    }
    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        // to be updated 
    }
}

void RHICfMassFitting::FillMass(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double mass)
{
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        if(mMassHistKinematic[runIdx][typeIdx][dleIdx][ptIdx][xfIdx]){
            mMassHistKinematic[runIdx][typeIdx][dleIdx][ptIdx][xfIdx] -> Fill(mass);
        }
    }
    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        // to be updated 
    }
}

void RHICfMassFitting::Fitting()
{
    TVirtualFitter::SetMaxIterations(100000);
    double bkg_xZero[2];
    double bkg_xSign[2];

    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                if(mOptContainer->GetParticleRunIdx() == kPi0Run){

                    FitMassHist(mMassHistAll[run][type][dle]);

                    double allLowerBoundary = mSignalFitter->GetParameter(0) - 3.*mSignalFitter->GetParameter(1);
                    double allUpperBoundary = mSignalFitter->GetParameter(0) + 3.*mSignalFitter->GetParameter(1);
                    mMassFitResults[run][type][dle].allMean = mSignalFitter->GetParameter(0);
                    mMassFitResults[run][type][dle].allSigma = mSignalFitter->GetParameter(1);
                    mMassFitResults[run][type][dle].allMassFitChi2 = mMassFitter -> GetChisquare()/mMassFitter->GetNDF();
                    mMassFitResults[run][type][dle].allMassAllCounts[0] = mMassFitter -> Integral(allLowerBoundary, allUpperBoundary);
                    mMassFitResults[run][type][dle].allMassSignalCounts[0] = mSignalFitter -> Integral(allLowerBoundary, allUpperBoundary);
                    mMassFitResults[run][type][dle].allMassBkgCounts[0] = mBkgFitter -> Integral(allLowerBoundary, allUpperBoundary);

                    TMatrixDSym covMatrix = mFitResult->GetCovarianceMatrix();
                    TMatrixDSym sigMatrix;
                    TMatrixDSym bkgMatrix;
                    covMatrix.GetSub(0, 2, sigMatrix);
                    covMatrix.GetSub(3, kMassFitParNum-1, bkgMatrix);

                    double allIntergralErr = mMassFitter->IntegralError(allLowerBoundary, allUpperBoundary, mFitResult->GetParams(), covMatrix.GetMatrixArray());
                    double sigIntegralErr = mSignalFitter->IntegralError(allLowerBoundary, allUpperBoundary, mSignalFitter->GetParameters() , sigMatrix.GetMatrixArray());
                    double bkgIntegralErr = mBkgFitter->IntegralError(allLowerBoundary, allUpperBoundary, mBkgFitter->GetParameters() , bkgMatrix.GetMatrixArray());
                    mMassFitResults[run][type][dle].allMassAllCounts[1] = allIntergralErr;
                    mMassFitResults[run][type][dle].allMassSignalCounts[1] = sigIntegralErr;
                    mMassFitResults[run][type][dle].allMassBkgCounts[1] = bkgIntegralErr;
            
                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = mMassFitter->GetParameter(i);
                    }

                    int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                    int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                    for(int pt=0; pt<ptNum; pt++){
                        for(int xf=0; xf<xfNum; xf++){
                            FitMassHist(mMassHistKinematic[run][type][dle][pt][xf]);

                            double lowerBoundary = mSignalFitter->GetParameter(0) - 3.*mSignalFitter->GetParameter(1);
                            double upperBoundary = mSignalFitter->GetParameter(0) + 3.*mSignalFitter->GetParameter(1);
                            mMassFitResults[run][type][dle].mean[pt][xf] = mSignalFitter->GetParameter(0);
                            mMassFitResults[run][type][dle].sigma[pt][xf] = mSignalFitter->GetParameter(1);
                            mMassFitResults[run][type][dle].massFitChi2[pt][xf] = mMassFitter -> GetChisquare()/mMassFitter->GetNDF();

                            mMassFitResults[run][type][dle].massAllCounts[0][pt][xf] = mMassFitter -> Integral(lowerBoundary, upperBoundary);
                            mMassFitResults[run][type][dle].massSignalCounts[0][pt][xf] = mSignalFitter -> Integral(lowerBoundary, upperBoundary);


                            TMatrixDSym covMatrix_ptxf = mFitResult->GetCovarianceMatrix();
                            TMatrixDSym sigMatrix_ptxf;
                            TMatrixDSym bkgMatrix_ptxf;
                            covMatrix_ptxf.GetSub(0, 2, sigMatrix_ptxf);
                            covMatrix_ptxf.GetSub(3, kMassFitParNum-1, bkgMatrix_ptxf);

                            double allIntergralErr_ptxf = mMassFitter->IntegralError(allLowerBoundary, allUpperBoundary, mFitResult->GetParams(), covMatrix_ptxf.GetMatrixArray());
                            double sigIntegralErr_ptxf = mSignalFitter->IntegralError(allLowerBoundary, allUpperBoundary, mSignalFitter->GetParameters() , sigMatrix_ptxf.GetMatrixArray());
                            mMassFitResults[run][type][dle].massAllCounts[1][pt][xf] = allIntergralErr_ptxf;
                            mMassFitResults[run][type][dle].massSignalCounts[1][pt][xf] = sigIntegralErr_ptxf;

                            memset(bkg_xZero, 0., sizeof(bkg_xZero));
                            memset(bkg_xSign, 0., sizeof(bkg_xSign));
                            bkg_xZero[0] = mBkgFitter -> GetX(0., lowerBoundary, mSignalFitter->GetParameter(0));
                            bkg_xSign[0] = mBkgFitter -> Eval(bkg_xZero[0]+0.5) - mBkgFitter -> Eval(bkg_xZero[0]-0.5);
                            bkg_xZero[1] = mBkgFitter -> GetX(0., mSignalFitter->GetParameter(0), upperBoundary);
                            bkg_xSign[1] = mBkgFitter -> Eval(bkg_xZero[1]+0.5) - mBkgFitter -> Eval(bkg_xZero[1]-0.5);
                            
                            double bkgCount = 0.;
                            double bkgCountErr = 0.;
                            bool isXzeroIntegral = false;
                            for(int step=0; step<2; step++){
                                double bkg_integralBound = -1;
                                if(lowerBoundary < bkg_xZero[step] && bkg_xZero[step] < upperBoundary){
                                    bkg_integralBound = bkg_xZero[step];
                                    isXzeroIntegral = true;
                                }
                                if(bkg_integralBound < 0){continue;}

                                double bkgIntegral = 0.;
                                double bkgIntegralErr = 0.;
                                if(bkg_xSign[step] > 0.){
                                    bkgIntegral = mBkgFitter -> Integral(bkg_xZero[step], upperBoundary);
                                    bkgIntegralErr = mBkgFitter->IntegralError(bkg_xZero[step], allUpperBoundary, mBkgFitter->GetParameters() , bkgMatrix_ptxf.GetMatrixArray());
                                }
                                if(bkg_xSign[step] < 0.){
                                    bkgIntegral = mBkgFitter -> Integral(lowerBoundary, bkg_xZero[step]);
                                    bkgIntegralErr = mBkgFitter->IntegralError(lowerBoundary, bkg_xZero[step], mBkgFitter->GetParameters() , bkgMatrix_ptxf.GetMatrixArray());
                                }
                                bkgCount += bkgIntegral;
                                bkgCountErr += (bkgIntegralErr*bkgIntegralErr);
                            }
                            bkgCountErr = sqrt(bkgCountErr);

                            if(!isXzeroIntegral){
                                bkgCount = mBkgFitter -> Integral(lowerBoundary, upperBoundary);
                                bkgCountErr = mBkgFitter->IntegralError(lowerBoundary, upperBoundary, mBkgFitter->GetParameters() , bkgMatrix_ptxf.GetMatrixArray());
                            }
                            mMassFitResults[run][type][dle].massBkgCounts[0][pt][xf] = bkgCount;
                            mMassFitResults[run][type][dle].massBkgCounts[1][pt][xf] = bkgCountErr;

                            if(mMassFitResults[run][type][dle].massAllCounts[0][pt][xf] < 10){
                                mMassFitResults[run][type][dle].massAllCounts[1][pt][xf] = 0.;
                                mMassFitResults[run][type][dle].massSignalCounts[1][pt][xf] = 0.;
                                mMassFitResults[run][type][dle].massBkgCounts[1][pt][xf] = 0.;
                            }

                            for(int i=0; i<kMassFitParNum; i++){
                                mMassFitResults[run][type][dle].par[i][pt][xf] = mMassFitter->GetParameter(i);
                            }
                        }
                    }
                }
                if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
                    // to be updated 
                }
            }
        }   
    }
    cout << "RHICfMassFitting::Fitting() -- mass fitting has done. " << endl;
}

void RHICfMassFitting::SaveMassData()
{
    vector<RHICfTableMaker::TableData> table;
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                RHICfTableMaker::TableData data;
                data.runIdx = run;
                data.typeIdx = type;
                data.dleIdx = dle;
                data.ptIdx = -1;
                data.xfIdx = -1;

                data.values.clear();
                data.values.push_back(GetMassMean(run, type, dle));
                data.values.push_back(GetMassSigma(run, type, dle));
                data.values.push_back(GetMassAllCounts(false, run, type, dle));
                data.values.push_back(GetMassAllCounts(true, run, type, dle));
                data.values.push_back(GetMassSignalCounts(false, run, type, dle));
                data.values.push_back(GetMassSignalCounts(true, run, type, dle));
                data.values.push_back(GetMassBkgCounts(false, run, type, dle));
                data.values.push_back(GetMassBkgCounts(true, run, type, dle));
                data.values.push_back(GetMassFitChi2(run, type, dle));

                for(int i=0; i<kMassFitParNum; i++){
                    data.values.push_back(mMassFitResults[run][type][dle].allPar[i]);
                }

                table.push_back(data);

                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        RHICfTableMaker::TableData data_kinematics;
                        data_kinematics.runIdx = run;
                        data_kinematics.typeIdx = type;
                        data_kinematics.dleIdx = dle;
                        data_kinematics.ptIdx = pt;
                        data_kinematics.xfIdx = xf;

                        data_kinematics.values.clear();
                        data_kinematics.values.push_back(GetMassMean(run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassSigma(run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassAllCounts(false, run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassAllCounts(true, run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassSignalCounts(false, run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassSignalCounts(true, run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassBkgCounts(false, run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassBkgCounts(true, run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassFitChi2(run, type, dle, pt, xf));

                        for(int i=0; i<kMassFitParNum; i++){
                            data_kinematics.values.push_back(mMassFitResults[run][type][dle].par[i][pt][xf]);
                        }
                        table.push_back(data_kinematics);
                    }
                }       
            }
        }
    }
    mTableMaker -> SaveTable("MassFit"+mTableName, table);

    if(!mIsFindMassFile){
        mMassFile ->cd();

        for(int run=0; run<kRunNum; run++){
            if(run != mOptContainer->GetRunType()){continue;}
            
            for(int type=0; type<kTypeNum; type++){
                for(int dle=0; dle<kDLENum; dle++){
                    int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                    int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                    if(ptNum > 0 && xfNum > 0){
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                mMassHistKinematic[run][type][dle][pt][xf] -> Write();
                            }
                        }
                    }
                }
            }
        }
        mMassFile -> Close();

    }
}

bool RHICfMassFitting::IsExistMassHistFile(){return mIsFindMassFile;}

int RHICfMassFitting::GetMassTableFlag()
{    
    int particleRunIdx = mOptContainer -> GetParticleRunIdx();
    if(particleRunIdx == kNeutronRun || particleRunIdx == kGammaRun){return kExistTable;}
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("MassFit"+mTableName, runIdx, 0, 3, 1, 1, 0) >= 0){
        return kExistTable;
    }
    if(mTableMaker -> GetTableData("MassFit"+mTableName, runIdx, 0, 3, -1, -1, 0) >= 0){
        return kExistPartOfTable;
    }
    return kNotExist;
}

double RHICfMassFitting::GetMassMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMean;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].mean[ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassSigma(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allSigma;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].sigma[ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassLowerBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return GetMassMean(runIdx, typeIdx, dleIdx, ptIdx, xfIdx) - 3. * GetMassSigma(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassUpperBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return GetMassMean(runIdx, typeIdx, dleIdx, ptIdx, xfIdx) + 3. * GetMassSigma(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassBkgLowerBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return GetMassMean(runIdx, typeIdx, dleIdx, ptIdx, xfIdx) - 5. * GetMassSigma(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassBkgUpperBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return GetMassMean(runIdx, typeIdx, dleIdx, ptIdx, xfIdx) + 5. * GetMassSigma(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassAllCounts(bool isErr, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassAllCounts[isErr];}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massAllCounts[isErr][ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassSignalCounts(bool isErr, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassSignalCounts[isErr];}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massSignalCounts[isErr][ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassBkgCounts(bool isErr, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassBkgCounts[isErr];}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massBkgCounts[isErr][ptIdx][xfIdx];
}

double RHICfMassFitting::GetBSRatio(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double bsRatio = GetMassBkgCounts(false, runIdx, typeIdx, dleIdx, ptIdx, xfIdx)/GetMassSignalCounts(false, runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(bsRatio < 0){bsRatio = 0.;}
    return bsRatio;
}

double RHICfMassFitting::GetBSRatioErr(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double signalCounts = GetMassSignalCounts(false, runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double signalCountErr = GetMassSignalCounts(true, runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgCounts = GetMassBkgCounts(false, runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgCountErr = GetMassBkgCounts(true, runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bsRatio = GetBSRatio(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double err = fabs(bsRatio)*sqrt(pow(double(signalCountErr/signalCounts), 2.) + pow(double(bkgCountErr/bkgCounts), 2.));
    return err;
}

double RHICfMassFitting::GetMassFitChi2(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassFitChi2;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massFitChi2[ptIdx][xfIdx];
}

double* RHICfMassFitting::GetMassFitPar(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double par[10];
    if(ptIdx == -1 && xfIdx == -1){
        for(int i=0; i<kMassFitParNum; i++){
            par[i] = mMassFitResults[runIdx][typeIdx][dleIdx].allPar[i];
        }
        return par;
    }
    for(int i=0; i<kMassFitParNum; i++){
        par[i] = mMassFitResults[runIdx][typeIdx][dleIdx].par[i][ptIdx][xfIdx];
    }
    return par;
}

TH1D* RHICfMassFitting::GetMassHist(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        return mMassHistAll[runIdx][typeIdx][dleIdx];
    }
    return mMassHistKinematic[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

double RHICfMassFitting::GetBSRatio_GPR(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mMassFitResults[runIdx][typeIdx][dleIdx].massAllCounts[0][ptIdx][xfIdx];
}

double RHICfMassFitting::GetBSRatioLower_GPR(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mMassFitResults[runIdx][typeIdx][dleIdx].massAllCounts[1][ptIdx][xfIdx];
}

double RHICfMassFitting::GetBSRatioUpper_GPR(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mMassFitResults[runIdx][typeIdx][dleIdx].massSignalCounts[0][ptIdx][xfIdx];
}

bool RHICfMassFitting::FitMassHist(TH1D* hist)
{
    for(int i=0; i<kMassFitParNum; i++){
        mMassFitter -> SetParameter(i, 0.);
        mMassFitter -> SetParError(i, 0.);
        if(i < mSigFitParNum){
            mSignalFitter -> SetParameter(i, 0.);
            mSignalFitter -> SetParError(i, 0.);
        }
        if(i < mBkgFitParNum){
            mBkgFitter -> SetParameter(i, 0.);
            mBkgFitter -> SetParError(i, 0.);
        }
    }
    mMassFitter -> SetChisquare(0.);
    
    double entryPeak = hist -> GetBinContent(int(135./hist->GetBinWidth(1)));
    double entryBkg = hist -> GetBinContent(int(50./hist->GetBinWidth(1)));
    if(entryPeak < 8 || entryBkg > entryPeak){return 0;}

    double fitLower, fitUpper;
    FindFitRange(hist, fitLower, fitUpper);

    mMassFitter -> SetRange(fitLower, fitUpper);
    mSignalFitter_drawing -> SetRange(fitLower, fitUpper);
    mBkgFitter_drawing -> SetRange(fitLower, fitUpper);

    mMassFitter -> SetParameter(0, 135);
    mMassFitter -> SetParameter(1, 10);
    mMassFitter -> SetParameter(2, hist->GetMaximum()/2.);
    mMassFitter -> SetParameter(3, 5);
    mMassFitter -> SetParameter(4, 290);
    mMassFitter -> SetParameter(9, 1.);
    mMassFitter -> SetParLimits(0, 125, 145);
    mMassFitter -> SetParLimits(1, 5, 15);

    mFitResult = hist -> Fit(mMassFitter, "S, QR");
    mFitResult = hist -> Fit(mMassFitter, "S, EQR");

    const double* par = mFitResult -> GetParams();
    const double* parErr = mFitResult -> GetErrors();

    for(int i=0; i<kMassFitParNum; i++){
        if(i < mSigFitParNum){mSignalFitter_drawing -> FixParameter(i, par[i]);}
        if(i < mBkgFitParNum){mBkgFitter_drawing -> FixParameter(i, par[i+mSigFitParNum]);}
    }

    double bkg_xZero1 = mBkgFitter_drawing -> GetX(0., par[0]-3.*par[1], par[0]);
    double bkg_xZero2 = mBkgFitter_drawing -> GetX(0., par[0], par[0]+3.*par[1]);
    double bkgCount = mBkgFitter_drawing->Integral(par[0]-3.*par[1], par[0]+3.*par[1]);

    bool reFitFlag = false;
    if(bkgCount < 0){reFitFlag = true;}
    if(par[0]-3.*par[1] < bkg_xZero1 && bkg_xZero1 < par[0]+3.*par[1]){reFitFlag = true;}
    if(par[0]-3.*par[1] < bkg_xZero2 && bkg_xZero2 < par[0]+3.*par[1]){reFitFlag = true;}
    if(reFitFlag){
        mMassFitter -> SetRange(fitLower-3., fitUpper+3.);
        mFitResult = hist -> Fit(mMassFitter, "S, EQR");
    }
    const double* par2 = mFitResult -> GetParams();
    const double* parErr2 = mFitResult -> GetErrors();

    for(int i=0; i<kMassFitParNum; i++){
        if(i < mSigFitParNum){mSignalFitter_drawing -> FixParameter(i, par2[i]);}
        if(i < mBkgFitParNum){mBkgFitter_drawing -> FixParameter(i, par2[i+mSigFitParNum]);}
    }
    hist -> Fit(mSignalFitter_drawing, "BQR+");
    hist -> Fit(mBkgFitter_drawing, "BQR+");

    for(int i=0; i<kMassFitParNum; i++){
        if(i < mSigFitParNum){
            mSignalFitter -> SetParameter(i, par2[i]);
            mSignalFitter -> SetParError(i, parErr2[i]);
        }
        if(i < mBkgFitParNum){
            mBkgFitter -> SetParameter(i, par2[i+mSigFitParNum]);
            mBkgFitter -> SetParError(i, parErr2[i+mSigFitParNum]);
        }
    }
    return 1;
}

void RHICfMassFitting::FindFitRange(TH1D* hist, double& lower, double& upper)
{
    int startBin1 = int(50/hist->GetBinWidth(1));
    int lowerBin = int(100/hist->GetBinWidth(1));
    lower = 0.;
    for(int bin=startBin1; bin<lowerBin; bin++){
        double entrySum = hist -> Integral(startBin1, bin+1);
        if(entrySum > 5.){
            lower = bin*hist->GetBinWidth(1);
            break;
        }
    }
    if(lower < 10){
        lowerBin = 103.;
    }

    int startBin2 = int(280/hist->GetBinWidth(1));
    int UpperBin = int(170/hist->GetBinWidth(1));
    upper = 0.;
    for(int bin=startBin2; bin>UpperBin; bin--){
        double entrySum = hist -> Integral(bin-1, startBin2);
        if(entrySum > 5.){
            upper = bin*hist->GetBinWidth(1);
            break;
        }
    }
    if(upper < 10){
        upper = 170.;
    }
    // lower = 20.;
    // upper = 250.;
}

void RHICfMassFitting::InitMassHistFile()
{

    TString dataPath = mOptContainer->GetDataPath();
    TString runName = mOptContainer->GetRunTypeName(mOptContainer->GetRunType());
    TString conditionName = mOptContainer->GetConditionName();
    TString particleType = mOptContainer->GetParticleRunName();
    TString fileName = dataPath+"/MassHist_"+runName+"_"+conditionName+"_"+particleType+".root";


    // Find a directory
    TList *listOfDirs;
    TObject *objDir;

    TSystemDirectory dir("dir", dataPath);
    listOfDirs = dir.GetListOfFiles();
    TIter next(listOfDirs);
    
    mIsFindMassFile = false;
    while((objDir = next())){
        TSystemFile* dirPtr = dynamic_cast<TSystemFile*>(objDir);
        if(dirPtr && !dirPtr->IsDirectory()){
            TString dirName = dirPtr->GetName();
            if(dirName.Index(fileName) != -1){
                cout << "RHICfMassFitting::InitMassHistFile()-- Find a Mass histogram file " << dirName << endl;
                mIsFindMassFile = true;
            }
        }
    }

    if(mIsFindMassFile){
        mMassFile = new TFile(fileName, "read");
        for(int run=0; run<kRunNum; run++){
            TString runName = mOptContainer -> GetRunTypeName(run);
            for(int type=0; type<kTypeNum; type++){
                for(int dle=0; dle<kDLENum; dle++){
                    TString dleName = mOptContainer -> GetDLEName(dle);
                    mMassHistAll[run][type][dle] = (TH1D*)mMassFile -> Get(Form("pi0MassHistAll_%s_type%i_%s", runName.Data(), type, dleName.Data()));

                    int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                    int xfNum = mBinning -> GetXfBinNum(run, type, dle);

                    if(ptNum > 0 && xfNum > 0){
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                mMassHistKinematic[run][type][dle][pt][xf] = (TH1D*)mMassFile -> Get(Form("pi0MassHist_%s_type%i_%s_pt%i_xf%i", runName.Data(), type, dleName.Data(), pt, xf));
                            }
                        }
                    }
                }
            }
        }
    }
    if(!mIsFindMassFile){
        mMassFile = new TFile(fileName, "recreate");
    }
}

double RHICfMassFitting::Pi0SignalFitter(double* x, double* par)
{
    double signal = par[2]* TMath::Gaus(x[0], par[0], par[1]);
    return signal;
}

double RHICfMassFitting::Pi0BkgFitter(double* x, double* par)
{
    double bkg = pow(x[0]-par[0],2)*pow(x[0]-par[1],2)*(par[2]*pow(x[0],3) + par[3]*pow(x[0],2) + par[4]*x[0] + par[5]) + par[6];
    // double bkg = pow(x[0]-par[0],2)*pow(x[0]-par[1],2)*par[2]*pow(x[0],3) + par[6];
    return bkg;
}

// 137.043 9.42142 65.6383         107.974 248.737 7.86411e-12 -4.72089e-09 9.40123e-07 -6.18809e-05 0.559052 

double RHICfMassFitting::Pi0MassFitter(double* x, double* par)
{
    double massfit = Pi0SignalFitter(x, par) + Pi0BkgFitter(x, &par[3]);
    return massfit;
}