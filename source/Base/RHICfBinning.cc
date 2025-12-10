#include "RHICfBinning.hh"

RHICfBinning::RHICfBinning(TString tableName) : mTableName(tableName) 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfBinning::~RHICfBinning()
{
}

void RHICfBinning::Init()
{
    mLocalBinNum = 50;
    mPtBins = 5;
    mXfBins = 4;

    if(mTableName.Sizeof() == 1){mTableName = mOptContainer->GetTableSubName();}
    mTableMaker -> InitTable("Binning"+mTableName);

    int particleRunType = mOptContainer -> GetParticleRunIdx();
    if(particleRunType == kPi0Run){Pi0GlobalBinning();}
    if(particleRunType == kNeutronRun){ NeutronGlobalBinning();}

    DefaultBinning();
   
    cout << "RHICfBinning::Init() -- Done." << endl;
}

void RHICfBinning::InitBinningData()
{
    if(GetBinningTableFlag() == kExistTable){
        for(int run=0; run<kRunNum; run++){
            for(int type=0; type<kTypeNum; type++){
                for(int dle=0; dle<kDLENum; dle++){
                    int ptBoundaryNum = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, 1, -1, -1);
                    int xfBoundaryNum = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, -1, 1, -1);
                    if(ptBoundaryNum > 0 && xfBoundaryNum > 0){    
                        mPtBoundary[run][type][dle].resize(ptBoundaryNum);
                        mXfBoundary[run][type][dle].resize(xfBoundaryNum);
                        for(int pt=0; pt<ptBoundaryNum; pt++){
                            mPtBoundary[run][type][dle][pt] = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, 1, -1, pt);
                        }
                        for(int xf=0; xf<xfBoundaryNum; xf++){
                            mXfBoundary[run][type][dle][xf] = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, -1, 1, xf);
                        }

                        int ptBinNum = ptBoundaryNum-1;
                        int xfBinNum = xfBoundaryNum-1;
                        mPtMean[run][type][dle].resize(ptBinNum, vector<double>(xfBinNum));
                        mXfMean[run][type][dle].resize(ptBinNum, vector<double>(xfBinNum));
                        for(int pt=0; pt<ptBinNum; pt++){
                            for(int xf=0; xf<xfBinNum; xf++){
                                mPtMean[run][type][dle][pt][xf] = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, pt, xf, 0);
                                mXfMean[run][type][dle][pt][xf] = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, pt, xf, 1);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RHICfBinning::InitHist()
{
    TString particleName = mOptContainer -> GetParticleRunName();
    for(int run=0; run<kRunNum; run++){
        TString runName = mOptContainer -> GetRunTypeName(run);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer -> GetDLEName(dle);
                
                mPtHist[run][type][dle] = new TH1D(Form("pTHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1, dleName.Data()), "", mLocalBinNum, 0., 1.);
                mXfHist[run][type][dle] = new TH1D(Form("xFHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1,dleName.Data()), "", mLocalBinNum, 0., 1.);
                mKinematicsHist[run][type][dle] = new TH2D(Form("KinematicsHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1, dleName.Data()), "", 100, 0., 1., 100, 0., 1.);

                mPtHist[run][type][dle] -> SetStats(0);
                mXfHist[run][type][dle] -> SetStats(0);
                mKinematicsHist[run][type][dle] -> SetStats(0);

                mPtHist[run][type][dle] -> SetTitle("p_{T} distribution; p_{T} [GeV/c]; N/x_{F}");
                mXfHist[run][type][dle] -> SetTitle("x_{F} distribution; x_{F}; N/p_{F}");
                mKinematicsHist[run][type][dle] -> SetTitle("Kinematics; x_{F}; p_{T} [GeV/c]");

                mPtHist[run][type][dle] -> Sumw2();
                mXfHist[run][type][dle] -> Sumw2();


                for(int ptBin=0; ptBin<GetPtBinNum(run, type, dle); ptBin++){
                    for(int xfBin=0; xfBin<GetXfBinNum(run, type, dle); xfBin++){
                        mPtMean[run][type][dle][ptBin][xfBin] = 0.;
                        mXfMean[run][type][dle][ptBin][xfBin] = 0.;
                        mPtMeanNum[run][type][dle][ptBin][xfBin] = 0.;
                        mXfMeanNum[run][type][dle][ptBin][xfBin] = 0.;
                    }
                }
            }
        }
    }
}

void RHICfBinning::FillKinematics(int runIdx, int typeIdx, int dleIdx, double pt, double xf)
{
    mPtHist[runIdx][typeIdx][dleIdx] -> Fill(pt);
    mXfHist[runIdx][typeIdx][dleIdx] -> Fill(xf);
    mKinematicsHist[runIdx][typeIdx][dleIdx] -> Fill(xf, pt);

    for(int ptBin=0; ptBin<GetPtBinNum(runIdx, typeIdx, dleIdx); ptBin++){
        double ptLowerBoundary = GetPtBinBoundary(runIdx, typeIdx, dleIdx, ptBin);
        double ptUpperBoundary = GetPtBinBoundary(runIdx, typeIdx, dleIdx, ptBin+1);

        for(int xfBin=0; xfBin<GetXfBinNum(runIdx, typeIdx, dleIdx); xfBin++){
            double xfLowerBoundary = GetXfBinBoundary(runIdx, typeIdx, dleIdx, xfBin);
            double xfUpperBoundary = GetXfBinBoundary(runIdx, typeIdx, dleIdx, xfBin+1);
            
            if(ptLowerBoundary <= pt && pt < ptUpperBoundary){
                if(xfLowerBoundary <= xf && xf < xfUpperBoundary){
                    mPtMean[runIdx][typeIdx][dleIdx][ptBin][xfBin] += pt;
                    mXfMean[runIdx][typeIdx][dleIdx][ptBin][xfBin] += xf;
                    mPtMeanNum[runIdx][typeIdx][dleIdx][ptBin][xfBin] += 1.;
                    mXfMeanNum[runIdx][typeIdx][dleIdx][ptBin][xfBin] += 1.;
                }
            }
        }
    }
}

void RHICfBinning::Binning()
{
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                for(int ptBin=0; ptBin<GetPtBinNum(run, type, dle); ptBin++){
                    for(int xfBin=0; xfBin<GetPtBinNum(run, type, dle); xfBin++){
                        mPtMean[run][type][dle][ptBin][xfBin] /= mPtMeanNum[run][type][dle][ptBin][xfBin];
                        mXfMean[run][type][dle][ptBin][xfBin] /= mXfMeanNum[run][type][dle][ptBin][xfBin];  
                        if(mPtMeanNum[run][type][dle][ptBin][xfBin] <= 1){mPtMean[run][type][dle][ptBin][xfBin] = 0.;}                      
                        if(mXfMeanNum[run][type][dle][ptBin][xfBin] <= 1){mXfMean[run][type][dle][ptBin][xfBin] = 0.;}                      
                    }
                }
            }
        }
    }
}

void RHICfBinning::SaveBinningData()
{
    vector<RHICfTableMaker::TableData> table;
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                for(int ptxf=0; ptxf<2; ptxf++){
                    int binBoundaryNum = (ptxf==0)? GetPtBinNum(run, type, dle)+1 : GetXfBinNum(run, type, dle)+1;
                    int ptIdx = (ptxf==0)? 1 : -1;
                    int xfIdx = (ptxf==0)? -1 : 1;

                    RHICfTableMaker::TableData data;
                    data.runIdx = run;
                    data.typeIdx = type;
                    data.dleIdx = dle;
                    data.ptIdx = ptIdx;
                    data.xfIdx = xfIdx;

                    for(int i=0; i<binBoundaryNum; i++){
                        double boundary = (ptxf==0)? GetPtBinBoundary(run, type, dle, i) : GetXfBinBoundary(run, type, dle, i);
                        data.values.push_back(boundary);
                    }
                    table.push_back(data);
                } 

                int ptNum = GetPtBinNum(run, type, dle);
                int xfNum = GetXfBinNum(run, type, dle);
                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        RHICfTableMaker::TableData data;
                        data.runIdx = run;
                        data.typeIdx = type;
                        data.dleIdx = dle;
                        data.ptIdx = pt;
                        data.xfIdx = xf;

                        double ptMean = GetPtBinMean(run, type, dle, pt, xf);
                        double xfMean = GetXfBinMean(run, type, dle, pt, xf);

                        data.values.push_back(ptMean);
                        data.values.push_back(xfMean);
                        table.push_back(data);   
                    }
                }
            }
        }
    }

    mTableMaker ->SaveTable("Binning"+mTableName, table);
}

int RHICfBinning::GetBinningTableFlag()
{
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Binning"+mTableName, runIdx, 0, 0, 1, -1, 1) > 0.){
        return kExistTable;
    }
    return kNotExist;
}

int RHICfBinning::GetPtBinNum(int runIdx, int typeIdx, int dleIdx)
{
    return int(mPtBoundary[runIdx][typeIdx][dleIdx].size()-1);
}

int RHICfBinning::GetXfBinNum(int runIdx, int typeIdx, int dleIdx)
{
    return int(mXfBoundary[runIdx][typeIdx][dleIdx].size()-1);
}

double RHICfBinning::GetPtBinBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx)
{
    if(GetPtBinNum(runIdx, typeIdx, dleIdx) < 1){return -999.;}
    return mPtBoundary[runIdx][typeIdx][dleIdx][ptIdx];
}

double RHICfBinning::GetXfBinBoundary(int runIdx, int typeIdx, int dleIdx, int xfIdx)
{
    if(GetXfBinNum(runIdx, typeIdx, dleIdx) < 1){return -999.;}
    return mXfBoundary[runIdx][typeIdx][dleIdx][xfIdx];
}

double RHICfBinning::GetPtBinMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(mPtMean[runIdx][typeIdx][dleIdx].size() < 1){return -999.;}
    if(mPtMean[runIdx][typeIdx][dleIdx][ptIdx].size() < 1){return -999.;}
    return mPtMean[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

double RHICfBinning::GetXfBinMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(mXfMean[runIdx][typeIdx][dleIdx].size() < 1){return -999.;}
    if(mXfMean[runIdx][typeIdx][dleIdx][ptIdx].size() < 1){return -999.;}
    return mXfMean[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

int RHICfBinning::GetGlobalPtBinNum()
{
    return int(mGlobalPtBoundary.size()-1);
}

int RHICfBinning::GetGlobalXfBinNum()
{
    return int(mGlobalXfBoundary.size()-1);
}

double RHICfBinning::GetGlobalPtBinBoundary(int ptIdx)
{
    if(GetGlobalPtBinNum() < 0){return -999.;}
    return mGlobalPtBoundary[ptIdx];
}

double RHICfBinning::GetGlobalXfBinBoundary(int xfIdx)
{
    if(GetGlobalXfBinNum() < 0){return -999.;}
    return mGlobalXfBoundary[xfIdx];
}

TH1D* RHICfBinning::GetPtHist(int runIdx, int typeIdx, int dleIdx){return mPtHist[runIdx][typeIdx][dleIdx];}
TH1D* RHICfBinning::GetXfHist(int runIdx, int typeIdx, int dleIdx){return mXfHist[runIdx][typeIdx][dleIdx];}
TH2D* RHICfBinning::GetKinematicsHist(int runIdx, int typeIdx, int dleIdx){return mKinematicsHist[runIdx][typeIdx][dleIdx];}

void RHICfBinning::Pi0GlobalBinning()
{
    mGlobalPtBoundary.clear();
    mGlobalXfBoundary.clear();

    // RHICf Collaboration, PRL 124, 252501 (2020)
    mGlobalXfBoundary.push_back(0.25);
    mGlobalXfBoundary.push_back(0.34);
    mGlobalXfBoundary.push_back(0.44);
    mGlobalXfBoundary.push_back(0.58);
    mGlobalXfBoundary.push_back(1.00);

    mGlobalPtBoundary.push_back(0.00);
    mGlobalPtBoundary.push_back(0.07);
    mGlobalPtBoundary.push_back(0.19);
    mGlobalPtBoundary.push_back(0.30);
    mGlobalPtBoundary.push_back(0.50);
    mGlobalPtBoundary.push_back(0.69);
    mGlobalPtBoundary.push_back(1.00);
}

void RHICfBinning::NeutronGlobalBinning()
{
    mGlobalPtBoundary.clear();
    mGlobalXfBoundary.clear();

    // RHICf Collaboration, PRD 109, 012003 (2024)
    mGlobalXfBoundary.push_back(0.2);
    mGlobalXfBoundary.push_back(0.4);
    mGlobalXfBoundary.push_back(0.6);
    mGlobalXfBoundary.push_back(0.8);
    mGlobalXfBoundary.push_back(1.0);

    mGlobalPtBoundary.push_back(0.00);
    mGlobalPtBoundary.push_back(0.1);
    mGlobalPtBoundary.push_back(0.2);
    mGlobalPtBoundary.push_back(0.35);
    mGlobalPtBoundary.push_back(0.55);
    mGlobalPtBoundary.push_back(0.85);
}

void RHICfBinning::DefaultBinning()
{
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){

                int globalPtNum = GetGlobalPtBinNum()+1;
                for(int pt=0; pt<globalPtNum; pt++){
                    mPtBoundary[run][type][dle].push_back(GetGlobalPtBinBoundary(pt));
                }

                int globalXfNum = GetGlobalXfBinNum()+1;
                for(int xf=0; xf<globalXfNum; xf++){
                    mXfBoundary[run][type][dle].push_back(GetGlobalXfBinBoundary(xf));
                }

                mPtMean[run][type][dle].resize(globalPtNum, vector<double>(globalXfNum));
                mXfMean[run][type][dle].resize(globalPtNum, vector<double>(globalXfNum));

                mPtMeanNum[run][type][dle].resize(globalPtNum, vector<double>(globalXfNum));
                mXfMeanNum[run][type][dle].resize(globalPtNum, vector<double>(globalXfNum));
            }
        }
    }
}