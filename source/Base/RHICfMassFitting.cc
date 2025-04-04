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
        double massFitLowerBoundary = 25.;
        double massFitUpperBoundary = 200.;
        double bkgFitLowerBoundary = 25.;
        double bkgFitUpperBondary = 200.;

        if(!mMassFitter){
            mMassFitter = new TF1("MassFitter", this, &RHICfMassFitting::Pi0MassFitter, massFitLowerBoundary, massFitUpperBoundary, kMassFitParNum);
            mMassFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
            mMassFitter -> SetLineColor(kBlack);
            mMassFitter -> SetLineWidth(1.5);
        }
        if(!mSignalFitter){
            mSignalFitter = new TF1("SignalFitter", this, &RHICfMassFitting::Pi0SignalFitter, massFitLowerBoundary, massFitUpperBoundary, 3);
            mSignalFitter -> SetParameters(0, 0, 0);
            mSignalFitter -> SetLineColor(kRed);
            mSignalFitter -> SetLineWidth(1.5);
        }
        if(!mBkgFitter){
            mBkgFitter = new TF1("BkgFitter", this, &RHICfMassFitting::Pi0BkgFitter, bkgFitLowerBoundary, bkgFitUpperBondary, 7);
            mBkgFitter -> SetParameters(0, 0, 0, 0, 0, 0, 0);
            mBkgFitter -> SetLineColor(kBlue);
            mBkgFitter -> SetLineWidth(1.5);
        }        
    }

    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        // to be updated
    }   

    mTableMaker -> InitTable("Mass");
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
                        massBins = 250;
                        massLowerBoundary = 0.;
                        massUpperBoundary = 250.;

                        if(!mMassHistAll[run][type][dle]){
                            mMassHistAll[run][type][dle] = new TH1D(Form("pi0MassHistAll_%s_type%i_%s", runName.Data(), type, dleName.Data()), "", massBins, massLowerBoundary, massUpperBoundary);
                            mMassHistAll[run][type][dle] -> SetStats(0);
                            mMassHistAll[run][type][dle] -> SetTitle(Form("%s Invariant Mass; M_{#gamma#gamma} [MeV/c^{2}]; Counts", (mOptContainer->GetParticleRunName()).Data()));
                            mMassFitResults[run][type][dle].Resize(0, 0);  
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
                            mMassFitResults[run][type][dle].Resize(ptNum, xfNum);  
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

void RHICfMassFitting::FillMass(int run, int type, int dle, double mass)
{
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        if(mMassHistAll[run][type][dle]){
            mMassHistAll[run][type][dle] -> Fill(mass);
        }
    }
    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        // to be updated 
    }
}

void RHICfMassFitting::FillMass(int run, int type, int dle, int ptIdx, int xfIdx, double mass)
{
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        if(mMassHistKinematic[run][type][dle][ptIdx][xfIdx]){
            mMassHistKinematic[run][type][dle][ptIdx][xfIdx] -> Fill(mass);
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

                    mMassFitResults[run][type][dle].allMean = par[0];
                    mMassFitResults[run][type][dle].allSigma = par[1];
                    for(int i=0; i<kMassFitParNum; i++){
                        mMassFitResults[run][type][dle].allPar[i] = par[i];
                    }

                    int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                    int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                    for(int pt=0; pt<ptNum; pt++){
                        for(int xf=0; xf<xfNum; xf++){
                            if(mMassHistKinematic[run][type][dle][pt][xf]->GetEntries() < 200){continue;}

                            double* par = FitMassHist(mMassHistKinematic[run][type][dle][pt][xf]);

                            mMassFitResults[run][type][dle].mean[pt][xf] = par[0];
                            mMassFitResults[run][type][dle].sigma[pt][xf] = par[1];
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

double RHICfMassFitting::GetMassMean(int run, int type, int dle, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        if(GetMassTableFlag() == kNotExist){return mMassFitResults[run][type][dle].allMean;}
        return mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 0);
    }
    if(GetMassTableFlag() != kExistTable){return mMassFitResults[run][type][dle].mean[ptIdx][xfIdx];}
    return mTableMaker->GetTableData("Mass", run, type, dle, ptIdx, xfIdx, 0);
}

double RHICfMassFitting::GetMassSigma(int run, int type, int dle, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        if(GetMassTableFlag() == kNotExist){return mMassFitResults[run][type][dle].allSigma;}
        return mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 1);
    }
    if(GetMassTableFlag() != kExistTable){return mMassFitResults[run][type][dle].sigma[ptIdx][xfIdx];}
    return mTableMaker->GetTableData("Mass", run, type, dle, ptIdx, xfIdx, 1);
}

double RHICfMassFitting::GetMassLowerBoundary(int run, int type, int dle, int ptIdx, int xfIdx)
{
    return GetMassMean(run, type, dle, ptIdx, xfIdx) - 3. * GetMassSigma(run, type, dle, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassUpperBoundary(int run, int type, int dle, int ptIdx, int xfIdx)
{
    return GetMassMean(run, type, dle, ptIdx, xfIdx) + 3. * GetMassSigma(run, type, dle, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassBkgLowerBoundary(int run, int type, int dle, int ptIdx, int xfIdx)
{
    return GetMassMean(run, type, dle, ptIdx, xfIdx) - 5. * GetMassSigma(run, type, dle, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassBkgUpperBoundary(int run, int type, int dle, int ptIdx, int xfIdx)
{
    return GetMassMean(run, type, dle, ptIdx, xfIdx) + 5. * GetMassSigma(run, type, dle, ptIdx, xfIdx);  
}

double RHICfMassFitting::GetMassAllCounts(int run, int type, int dle, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        for(int i=0; i<kMassFitParNum; i++){
            double par = 0.;
            if(GetMassTableFlag() == kNotExist){par = mMassFitResults[run][type][dle].allPar[i];}
            else{par = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 5+i);}
            mMassFitter -> SetParameter(i, par);
        }
        return mMassFitter -> Integral(GetMassLowerBoundary(run, type, dle), GetMassUpperBoundary(run, type, dle));
    }

    for(int i=0; i<kMassFitParNum; i++){
        double par = 0.;
        if(GetMassTableFlag() == kExistPartOfTable){par = mMassFitResults[run][type][dle].par[i][ptIdx][xfIdx];}
        else{par = mTableMaker->GetTableData("Mass", run, type, dle, ptIdx, xfIdx, 5+i);}
        mMassFitter -> SetParameter(i, par);
    }
    return mMassFitter -> Integral(GetMassLowerBoundary(run, type, dle, ptIdx, xfIdx), GetMassUpperBoundary(run, type, dle, ptIdx, xfIdx));
}

double RHICfMassFitting::GetMassSignalCounts(int run, int type, int dle, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        for(int i=0; i<3; i++){
            double par = 0.;
            if(GetMassTableFlag() == kNotExist){par = mMassFitResults[run][type][dle].allPar[i];}
            else{par = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 5+i);}
            mSignalFitter -> SetParameter(i, par);
        }
        return mSignalFitter -> Integral(GetMassLowerBoundary(run, type, dle), GetMassUpperBoundary(run, type, dle));
    }

    for(int i=0; i<3; i++){
        double par = 0.;
        if(GetMassTableFlag() == kExistPartOfTable){par = mMassFitResults[run][type][dle].par[i][ptIdx][xfIdx];}
        else{par = mTableMaker->GetTableData("Mass", run, type, dle, ptIdx, xfIdx, 5+i);}
        mSignalFitter -> SetParameter(i, par);
    }
    return mSignalFitter -> Integral(GetMassLowerBoundary(run, type, dle, ptIdx, xfIdx), GetMassUpperBoundary(run, type, dle, ptIdx, xfIdx));
}

double RHICfMassFitting::GetMassBkgCounts(int run, int type, int dle, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        for(int i=3; i<kMassFitParNum; i++){
            double par = 0.;
            if(GetMassTableFlag() == kNotExist){par = mMassFitResults[run][type][dle].allPar[i];}
            else{par = mTableMaker->GetTableData("Mass", run, type, dle, -1, -1, 5+i);}
            mBkgFitter -> SetParameter(i-3, par);
        }
        return mBkgFitter -> Integral(GetMassLowerBoundary(run, type, dle), GetMassUpperBoundary(run, type, dle));
    }

    for(int i=3; i<kMassFitParNum; i++){
        double par = 0.;
        if(GetMassTableFlag() == kExistPartOfTable){par = mMassFitResults[run][type][dle].par[i][ptIdx][xfIdx];}
        else{par = mTableMaker->GetTableData("Mass", run, type, dle, ptIdx, xfIdx, 5+i);}
        mBkgFitter -> SetParameter(i-3, par);
    }
    return mBkgFitter -> Integral(GetMassLowerBoundary(run, type, dle, ptIdx, xfIdx), GetMassUpperBoundary(run, type, dle, ptIdx, xfIdx));
}

TH1D* RHICfMassFitting::GetMassHist(int run, int type, int dle, int ptIdx, int xfIdx)
{
    if(ptIdx == -1 && xfIdx == -1){
        return mMassHistAll[run][type][dle];
    }
    return mMassHistKinematic[run][type][dle][ptIdx][xfIdx];
}

double* RHICfMassFitting::FitMassHist(TH1D* hist)
{
    double par[10];
    memset(par, 0., sizeof(par));

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
