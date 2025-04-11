#include "RHICfMassFitting.hh"

#include <TVirtualFitter.h>

RHICfMassFitting::RHICfMassFitting() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
    mMassFitter = 0;
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
        double fitLowerBoundary = 40.;
        double fitUpperBoundary = 220.;

        if(!mMassFitter){
            mMassFitter = new TF1("MassFitter", this, &RHICfMassFitting::Pi0MassFitter, fitLowerBoundary, fitUpperBoundary, kMassFitParNum);
            mMassFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
            mMassFitter -> SetLineColor(kBlack);
            mMassFitter -> SetLineWidth(1.5);
        }
        if(!mSignalFitter){
            mSignalFitter = new TF1("SignalFitter", this, &RHICfMassFitting::Pi0SignalFitter, fitLowerBoundary, fitUpperBoundary, 3);
            mSignalFitter -> SetParameters(0, 0, 0);
            mSignalFitter -> SetLineColor(kRed);
            mSignalFitter -> SetLineWidth(1.5);
        }
        if(!mBkgFitter){
            mBkgFitter = new TF1("BkgFitter", this, &RHICfMassFitting::Pi0BkgFitter, fitLowerBoundary, fitUpperBoundary, 7);
            mBkgFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0);
            mBkgFitter -> SetLineColor(kBlue);
            mBkgFitter -> SetLineWidth(1.5);
        }        
    }

    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        // to be updated
    }   

    mTableMaker -> InitTable("Mass");
    cout << "RHICfMassFitting::Init() -- Done." << endl;
}

void RHICfMassFitting::InitMassData()
{
    int flag = GetMassTableFlag();
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptBinNum = mTableMaker->GetTableData("Binning", run, type, dle, 1, -1, -1)-1;
                int xfBinNum = mTableMaker->GetTableData("Binning", run, type, dle, -1, 1, -1)-1;
                if(ptBinNum > 0 && xfBinNum > 0){
                    mMassFitResults[run][type][dle].Resize(ptBinNum, xfBinNum);
                }

                if(flag == kExistTable){
                    mMassFitResults[run][type][dle].allMean = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 0);
                    mMassFitResults[run][type][dle].allSigma = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 1);
                    mMassFitResults[run][type][dle].allMassAllCounts = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 2);
                    mMassFitResults[run][type][dle].allMassSignalCounts = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 3);
                    mMassFitResults[run][type][dle].allMassBkgCounts = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 4);
                    mMassFitResults[run][type][dle].allMassPeakDifference = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 5);
                    mMassFitResults[run][type][dle].allMassFitChi2 = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 6);
                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 7+i);
                    }

                    for(int pt=0; pt<ptBinNum; pt++){
                        for(int xf=0; xf<xfBinNum; xf++){
                            mMassFitResults[run][type][dle].mean[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 0);
                            mMassFitResults[run][type][dle].sigma[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 1);
                            mMassFitResults[run][type][dle].massAllCounts[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 2);
                            mMassFitResults[run][type][dle].massSignalCounts[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 3);
                            mMassFitResults[run][type][dle].massBkgCounts[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 4);
                            mMassFitResults[run][type][dle].massPeakDifference[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 5);
                            mMassFitResults[run][type][dle].massFitChi2[pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 6);

                            for(int i=0; i<kMassFitParNum; i++){
                                mMassFitResults[run][type][dle].par[i][pt][xf] = mTableMaker->GetTableData("Mass", run, type, dle, pt, xf, 7+i);
                            }
                        }
                    }
                }
                else if(flag == kExistPartOfTable){
                    mMassFitResults[run][type][dle].allMean = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 0);
                    mMassFitResults[run][type][dle].allSigma = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 1);
                    mMassFitResults[run][type][dle].allMassAllCounts = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 2);
                    mMassFitResults[run][type][dle].allMassSignalCounts = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 3);
                    mMassFitResults[run][type][dle].allMassBkgCounts = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 4);
                    mMassFitResults[run][type][dle].allMassPeakDifference = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 5);
                    mMassFitResults[run][type][dle].allMassFitChi2 = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 6);
                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 7+i);
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
                        massBins = 100;
                        massLowerBoundary = 0.;
                        massUpperBoundary = 300.;

                        if(!mMassHistAll[run][type][dle]){
                            mMassHistAll[run][type][dle] = new TH1D(Form("pi0MassHistAll_%s_type%i_%s", runName.Data(), type, dleName.Data()), "", massBins, massLowerBoundary, massUpperBoundary);
                            mMassHistAll[run][type][dle] -> SetStats(0);
                            mMassHistAll[run][type][dle] -> SetTitle(Form("%s Invariant Mass; M_{#gamma#gamma} [MeV/c^{2}]; Counts", (mOptContainer->GetParticleRunName()).Data()));
                        }
                        mMassHistAll[run][type][dle] -> Clear("ICESM");
                        
                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);

                        if(ptNum > 0 && xfNum > 0){
                            if(mMassHistKinematic[run][type][dle].size() == 0){
                                mMassHistKinematic[run][type][dle].resize(ptNum, vector<TH1D*>(xfNum));
                            }
                            for(int pt=0; pt<ptNum; pt++){
                                for(int xf=0; xf<xfNum; xf++){
                                    if(!mMassHistKinematic[run][type][dle][pt][xf]){
                                        mMassHistKinematic[run][type][dle][pt][xf] = new TH1D(Form("pi0MassHist_%s_type%i_%s_pt%i_xf%i", runName.Data(), type, dleName.Data(), pt, xf), "", massBins, massLowerBoundary, massUpperBoundary);
                                        mMassHistKinematic[run][type][dle][pt][xf] -> SetStats(0);
                                        mMassHistKinematic[run][type][dle][pt][xf] -> SetTitle(Form("%s Invariant Mass; M_{#gamma#gamma} [MeV/c^{2}]; Counts", (mOptContainer->GetParticleRunName()).Data()));
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
    TVirtualFitter::SetMaxIterations(50000);

    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                if(mOptContainer->GetParticleRunIdx() == kPi0Run){
                    if(mMassHistAll[run][type][dle]->GetEntries() < 200){continue;}

                    double* par = FitMassHist(mMassHistAll[run][type][dle]);

                    double allLowerBoundary = par[0] - 3.*par[1];
                    double allUpperBoundary = par[0] + 3.*par[1];
                    int allLowerBin = mMassHistAll[run][type][dle] -> FindBin(allLowerBoundary);
                    int allUpperBin = mMassHistAll[run][type][dle] -> FindBin(allUpperBoundary);
                    mMassHistAll[run][type][dle] -> GetXaxis() -> SetRangeUser(allLowerBin, allUpperBin); 
                    int allPeakBin = mMassHistAll[run][type][dle]->GetMaximumBin();
                    double allPeakEntry = mMassHistAll[run][type][dle]->GetBinContent(allPeakBin);
                    double allFitPeakEntry = mMassFitter -> Eval(mMassHistAll[run][type][dle]->GetBinCenter(allPeakBin));
                    double allDiffentPeakEntry = fabs(allPeakEntry - allFitPeakEntry);

                    mMassHistAll[run][type][dle] -> GetXaxis() -> UnZoom();

                    mMassFitResults[run][type][dle].allMean = par[0];
                    mMassFitResults[run][type][dle].allSigma = par[1];
                    mMassFitResults[run][type][dle].allMassAllCounts = mMassFitter -> Integral(allLowerBoundary, allUpperBoundary);
                    mMassFitResults[run][type][dle].allMassSignalCounts = mSignalFitter -> Integral(allLowerBoundary, allUpperBoundary);
                    mMassFitResults[run][type][dle].allMassBkgCounts = mBkgFitter -> Integral(allLowerBoundary, allUpperBoundary);
                    mMassFitResults[run][type][dle].allMassPeakDifference = allDiffentPeakEntry;
                    mMassFitResults[run][type][dle].allMassFitChi2 = mMassFitter -> GetChisquare();

                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = par[i];
                    }

                    int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                    int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                    for(int pt=0; pt<ptNum; pt++){
                        for(int xf=0; xf<xfNum; xf++){
                            if(mMassHistKinematic[run][type][dle][pt][xf]->GetEntries() < 200){continue;}

                            double* par = FitMassHist(mMassHistKinematic[run][type][dle][pt][xf]);

                            double lowerBoundary = par[0] - 3.*par[1];
                            double upperBoundary = par[0] + 3.*par[1];
                            int lowerBin = mMassHistKinematic[run][type][dle][pt][xf] -> FindBin(lowerBoundary);
                            int upperBin = mMassHistKinematic[run][type][dle][pt][xf] -> FindBin(upperBoundary);
                            mMassHistKinematic[run][type][dle][pt][xf] -> GetXaxis() -> SetRangeUser(lowerBin, upperBin); 
                            int peakBin = mMassHistKinematic[run][type][dle][pt][xf]->GetMaximumBin();
                            double peakEntry = mMassHistKinematic[run][type][dle][pt][xf]->GetBinContent(peakBin);
                            double fitPeakEntry = mMassFitter -> Eval(mMassHistKinematic[run][type][dle][pt][xf]->GetBinCenter(peakBin));
                            double diffentPeakEntry = fabs(peakEntry - fitPeakEntry);

                            mMassHistKinematic[run][type][dle][pt][xf] -> GetXaxis() -> UnZoom();

                            mMassFitResults[run][type][dle].mean[pt][xf] = par[0];
                            mMassFitResults[run][type][dle].sigma[pt][xf] = par[1];
                            mMassFitResults[run][type][dle].massAllCounts[pt][xf] = mMassFitter -> Integral(par[0]-3.*par[1], par[0]+3.*par[1]);
                            mMassFitResults[run][type][dle].massSignalCounts[pt][xf] = mSignalFitter -> Integral(par[0]-3.*par[1], par[0]+3.*par[1]);
                            mMassFitResults[run][type][dle].massBkgCounts[pt][xf] = mBkgFitter -> Integral(par[0]-3.*par[1], par[0]+3.*par[1]);
                            mMassFitResults[run][type][dle].massPeakDifference[pt][xf] = diffentPeakEntry;
                            mMassFitResults[run][type][dle].massFitChi2[pt][xf] = mMassFitter -> GetChisquare();

                            for(int i=0; i<kMassFitParNum; i++){
                                mMassFitResults[run][type][dle].par[i][pt][xf] = par[i];
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
                data.values.push_back(GetMassAllCounts(run, type, dle));
                data.values.push_back(GetMassSignalCounts(run, type, dle));
                data.values.push_back(GetMassBkgCounts(run, type, dle));
                data.values.push_back(GetMassPeakDifference(run, type, dle));
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
                        data_kinematics.values.push_back(GetMassAllCounts(run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassSignalCounts(run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassBkgCounts(run, type, dle, pt, xf));
                        data_kinematics.values.push_back(GetMassPeakDifference(run, type, dle, pt, xf));
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
    mTableMaker -> SaveTable("Mass", table);
}

int RHICfMassFitting::GetMassTableFlag()
{    
    int particleRunIdx = mOptContainer -> GetParticleRunIdx();
    if(particleRunIdx == kNeutronRun || particleRunIdx == kGammaRun){return kExistTable;}
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Mass", runIdx, 0, 3, 1, 1, 0) >= 0){
        return kExistTable;
    }
    if(mTableMaker -> GetTableData("Mass", runIdx, 0, 3, -1, -1, 0) >= 0){
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

double RHICfMassFitting::GetMassAllCounts(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassAllCounts;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massAllCounts[ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassSignalCounts(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassSignalCounts;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massSignalCounts[ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassBkgCounts(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassBkgCounts;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massBkgCounts[ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassPeakDifference(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassPeakDifference;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massPeakDifference[ptIdx][xfIdx];
}

double RHICfMassFitting::GetMassFitPeak(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        for(int i=0; i<kMassFitParNum; i++){
            mMassFitter -> SetParameter(i, mMassFitResults[runIdx][typeIdx][dleIdx].allPar[i]);
        }
        return mMassFitter -> Eval(134.);
    }
    for(int i=0; i<kMassFitParNum; i++){
        mMassFitter -> SetParameter(i, mMassFitResults[runIdx][typeIdx][dleIdx].par[i][ptIdx][xfIdx]);
    }
    return mMassFitter -> Eval(134.);              
}

double RHICfMassFitting::GetMassFitChi2(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){return mMassFitResults[runIdx][typeIdx][dleIdx].allMassFitChi2;}
    return mMassFitResults[runIdx][typeIdx][dleIdx].massFitChi2[ptIdx][xfIdx];
}


TH1D* RHICfMassFitting::GetMassHist(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        return mMassHistAll[runIdx][typeIdx][dleIdx];
    }
    return mMassHistKinematic[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

double* RHICfMassFitting::FitMassHist(TH1D* hist)
{
    mMassFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    mSignalFitter -> SetParameters(0, 0, 0);
    mBkgFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0);

    mMassFitter -> SetParameter(0, 135);
    mMassFitter -> SetParameter(1, 10);
    mMassFitter -> SetParameter(2, hist->GetMaximum()/2.);
    mMassFitter -> SetParameter(3, 5);
    mMassFitter -> SetParameter(4, 290);
    mMassFitter -> SetParLimits(0, 125, 145);
    mMassFitter -> SetParLimits(1, 5, 15);

    hist -> Fit(mMassFitter, "QR");

    static double par[10];
    mMassFitter -> GetParameters(par);

    mSignalFitter -> FixParameter(0, par[0]);
    mSignalFitter -> FixParameter(1, par[1]);
    mSignalFitter -> FixParameter(2, par[2]);
    hist -> Fit(mSignalFitter, "QR+");

    mBkgFitter -> FixParameter(0, par[3]);
    mBkgFitter -> FixParameter(1, par[4]);
    mBkgFitter -> FixParameter(2, par[5]);
    mBkgFitter -> FixParameter(3, par[6]);
    mBkgFitter -> FixParameter(4, par[7]);
    mBkgFitter -> FixParameter(5, par[8]);
    mBkgFitter -> FixParameter(6, par[9]);
    hist -> Fit(mBkgFitter, "QR+");

    return par;
}

double RHICfMassFitting::Pi0SignalFitter(double* x, double* par)
{
    double signal = par[2]* TMath::Gaus(x[0], par[0], par[1], kTRUE);
    return signal;
}

double RHICfMassFitting::Pi0BkgFitter(double* x, double* par)
{
    double bkg = pow(x[0]-par[0],2)*pow(x[0]-par[1],2)*(par[2]*pow(x[0],3) + par[3]*pow(x[0],2) + par[4]*x[0] + par[5]) + par[6];
    return bkg;
}

double RHICfMassFitting::Pi0MassFitter(double* x, double* par)
{
    double massfit = Pi0SignalFitter(x, par) + Pi0BkgFitter(x, &par[3]);
    return massfit;
}
