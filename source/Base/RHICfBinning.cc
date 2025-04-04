#include "RHICfBinning.hh"

RHICfBinning::RHICfBinning() 
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

    Pi0GlobalBinning();

    TString particleName = mOptContainer -> GetParticleRunName();
    for(int run=0; run<kRunNum; run++){
        TString runName = mOptContainer -> GetRunTypeName(run);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer -> GetDLEName(dle);
                
                mPtHist[run][type][dle] = new TH1D(Form("pTHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1, dleName.Data()), "", mLocalBinNum, 0., 1.);
                mXfHist[run][type][dle] = new TH1D(Form("xFHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1,dleName.Data()), "", mLocalBinNum, 0., 1.);
                mPtXfBinHist[run][type][dle] = new TH2D(Form("ptxfBinHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1,dleName.Data()), "", mLocalBinNum, 0., 1., mLocalBinNum, 0., 1.);
                mKinematicsHist[run][type][dle] = new TH2D(Form("KinematicsHist_%s_%s_type%i_%s", particleName.Data(), runName.Data(), type+1, dleName.Data()), "", 100, 0., 1., 100, 0., 1.);

                mPtHist[run][type][dle] -> SetStats(0);
                mXfHist[run][type][dle] -> SetStats(0);
                mKinematicsHist[run][type][dle] -> SetStats(0);

                mPtHist[run][type][dle] -> SetTitle("p_{T} distribution; p_{T} [GeV/c]; N/x_{F}");
                mXfHist[run][type][dle] -> SetTitle("x_{F} distribution; x_{F}; N/p_{F}");
                mKinematicsHist[run][type][dle] -> SetTitle("Kinematics; x_{F}; p_{T} [GeV/c]");

                mPtHist[run][type][dle] -> Sumw2();
                mXfHist[run][type][dle] -> Sumw2();

                mPtBoundary[run][type][dle].clear();
                mXfBoundary[run][type][dle].clear();

                mPtMean[run][type][dle].clear();
                mXfMean[run][type][dle].clear();
            }
        }
    }

    mTableMaker -> InitTable("Binning");
}

void RHICfBinning::FillKinematics(int runIdx, int typeIdx, int dleIdx, double pt, double xf)
{
    mPtHist[runIdx][typeIdx][dleIdx] -> Fill(pt);
    mXfHist[runIdx][typeIdx][dleIdx] -> Fill(xf);
    mPtXfBinHist[runIdx][typeIdx][dleIdx] -> Fill(xf, pt);
    mKinematicsHist[runIdx][typeIdx][dleIdx] -> Fill(xf, pt);
}

void RHICfBinning::Binning()
{
    TH1D* projectionXf;
    TH1D* projectionPt;
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                vector<int> mXfBoundBinIdx;
                vector<int> mPtBoundBinIdx;
                mXfBoundBinIdx.clear();
                mPtBoundBinIdx.clear();

                // Xf Binning 
                double xfTotalEntry = 0.;
                for(int xf=0; xf<mLocalBinNum; xf++){
                    double entry = mXfHist[run][type][dle]->GetBinContent(xf+1);
                    if(entry < 1){continue;}

                    projectionPt = (TH1D*)mPtXfBinHist[run][type][dle]->ProjectionY(Form("project_%i_%i_%i_xf%i", run, type, dle, xf), xf+1, xf+1);
                    double ptMean = projectionPt -> GetMean();
                    if(ptMean < 0.00001){continue;}

                    double normEntry = entry/ptMean;
                    mXfHist[run][type][dle] -> SetBinContent(xf+1, normEntry);
                    xfTotalEntry += normEntry;
                }

                double xfNumByBin = xfTotalEntry/double(mXfBins);
                double tmpXfEntrySum = 0.;

                mXfBoundBinIdx.push_back(1);
                mXfBoundary[run][type][dle].push_back(0.);
                for(int xf=0; xf<mLocalBinNum; xf++){
                    double entry = mXfHist[run][type][dle] -> GetBinContent(xf+1);
                    tmpXfEntrySum += entry;

                    if(tmpXfEntrySum > xfNumByBin){
                        double binCenter = mXfHist[run][type][dle] -> GetBinCenter(xf+1);
                        double binWidth = mXfHist[run][type][dle] -> GetBinWidth(xf+1);
                        double binBondary = binCenter + binWidth/2.;

                        mXfBoundBinIdx.push_back(xf+1);
                        mXfBoundary[run][type][dle].push_back(binBondary);
                        tmpXfEntrySum = 0;
                    }
                }
                mXfBoundBinIdx.push_back(mLocalBinNum);
                mXfBoundary[run][type][dle].push_back(1.);


                // Pt Binning 
                double ptTotalEntry = 0.;
                for(int pt=0; pt<mLocalBinNum; pt++){
                    double entry = mPtHist[run][type][dle]->GetBinContent(pt+1);
                    if(entry < 1){continue;}

                    projectionXf = (TH1D*)mPtXfBinHist[run][type][dle]->ProjectionX(Form("project_%i_%i_%i_pt%i", run, type, dle, pt), pt+1, pt+1);
                    double xfMean = projectionXf -> GetMean();
                    if(xfMean < 0.00001){continue;}

                    double normEntry = entry/xfMean;
                    mPtHist[run][type][dle] -> SetBinContent(pt+1, normEntry);
                    ptTotalEntry += normEntry;
                }

                double ptNumByBin = ptTotalEntry/double(mPtBins);
                double tmpPtEntrySum = 0.;

                mPtBoundBinIdx.push_back(1);
                mPtBoundary[run][type][dle].push_back(0.);
                for(int pt=0; pt<mLocalBinNum; pt++){
                    double entry = mPtHist[run][type][dle] -> GetBinContent(pt+1);
                    tmpPtEntrySum += entry;

                    if(tmpPtEntrySum > ptNumByBin){
                        double binCenter = mPtHist[run][type][dle] -> GetBinCenter(pt+1);
                        double binWidth = mPtHist[run][type][dle] -> GetBinWidth(pt+1);
                        double binBondary = binCenter + binWidth/2.;

                        mPtBoundBinIdx.push_back(pt+1);
                        mPtBoundary[run][type][dle].push_back(binBondary);
                        tmpPtEntrySum = 0;
                    }
                }
                mPtBoundBinIdx.push_back(mLocalBinNum);
                mPtBoundary[run][type][dle].push_back(1.);


                // Find mean value in binning
                int ptNum = GetPtBinNum(run, type, dle);
                int xfNum = GetXfBinNum(run, type, dle);

                mPtMean[run][type][dle].resize(ptNum, vector<double>(xfNum));
                mXfMean[run][type][dle].resize(ptNum, vector<double>(xfNum));

                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        int xfBinIdx_lower = mXfBoundBinIdx[xf];
                        int ptBinIdx_lower = mPtBoundBinIdx[pt];
                        int xfBinIdx_upper = mXfBoundBinIdx[xf+1];
                        int ptBinIdx_upper = mPtBoundBinIdx[pt+1];

                        double entrySum = 0.;
                        double xfWeightSum = 0.;
                        double ptWeightSum = 0.;

                        for(int ptBin=ptBinIdx_lower; ptBin<ptBinIdx_upper; ptBin++){
                            for(int xfBin=xfBinIdx_lower; xfBin<xfBinIdx_upper; xfBin++){
                                double entry = mPtXfBinHist[run][type][dle] -> GetBinContent(xfBin, ptBin);
                                double ptValue = double(ptBin)/double(mLocalBinNum)+(1./double(mLocalBinNum*2));
                                double xfValue = double(xfBin)/double(mLocalBinNum)+(1./double(mLocalBinNum*2));

                                entrySum += entry;
                                ptWeightSum += (ptValue * entry);
                                xfWeightSum += (xfValue * entry);
                            }
                        }

                        ptWeightSum = ptWeightSum/entrySum;
                        xfWeightSum = xfWeightSum/entrySum;
                        if(entrySum < 10){
                            ptWeightSum = -1.;
                            xfWeightSum = -1.;
                        }

                        mPtMean[run][type][dle][pt][xf] = ptWeightSum;
                        mXfMean[run][type][dle][pt][xf] = xfWeightSum;
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

    mTableMaker ->SaveTable("Binning", table);
}

int RHICfBinning::GetBinningTableFlag()
{
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Binning", runIdx, 0, 0, 1, -1, 1) > 0.){
        return kExistTable;
    }
    return kNotExist;
}

int RHICfBinning::GetPtBinNum(int runIdx, int typeIdx, int dleIdx)
{
    if(GetBinningTableFlag() == kNotExist){return int(mPtBoundary[runIdx][typeIdx][dleIdx].size()-1);}
    return int(mTableMaker->GetTableData("Binning", runIdx, typeIdx, dleIdx, 1, -1, -1)-1);
}

int RHICfBinning::GetXfBinNum(int runIdx, int typeIdx, int dleIdx)
{
    if(GetBinningTableFlag() == kNotExist){return int(mXfBoundary[runIdx][typeIdx][dleIdx].size()-1);}
    return int(mTableMaker->GetTableData("Binning", runIdx, typeIdx, dleIdx, -1, 1, -1)-1);
}

double RHICfBinning::GetPtBinBoundary(int runIdx, int typeIdx, int dleIdx, int ptIdx)
{
    if(GetBinningTableFlag() == kNotExist){return mPtBoundary[runIdx][typeIdx][dleIdx][ptIdx];}
    return mTableMaker->GetTableData("Binning", runIdx, typeIdx, dleIdx, 1, -1, ptIdx);
}

double RHICfBinning::GetXfBinBoundary(int runIdx, int typeIdx, int dleIdx, int xfIdx)
{
    if(GetBinningTableFlag() == kNotExist){return mXfBoundary[runIdx][typeIdx][dleIdx][xfIdx];}
    return mTableMaker->GetTableData("Binning", runIdx, typeIdx, dleIdx, -1, 1, xfIdx);
}

double RHICfBinning::GetPtBinMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(GetBinningTableFlag() == kNotExist){return mPtMean[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];}
    return mTableMaker->GetTableData("Binning", runIdx, typeIdx, dleIdx, ptIdx, xfIdx, 0);
}

double RHICfBinning::GetXfBinMean(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    if(GetBinningTableFlag() == kNotExist){return mXfMean[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];}
    return mTableMaker->GetTableData("Binning", runIdx, typeIdx, dleIdx, ptIdx, xfIdx, 1);
}

int RHICfBinning::GetGlobalPtBinNum()
{
    if(GetBinningTableFlag() == kNotExist){return int(mGlobalPtBoundary.size()-1);}
    return int(mTableMaker->GetTableData("Binning", 0, 0, 0, 1, -1, -1)-1);
}

int RHICfBinning::GetGlobalXfBinNum()
{
    if(GetBinningTableFlag() == kNotExist){return int(mGlobalXfBoundary.size()-1);}
    return int(mTableMaker->GetTableData("Binning", 0, 0, 0, -1, 1, -1)-1);
}

double RHICfBinning::GetGlobalPtBinBoundary(int ptIdx)
{
    if(GetBinningTableFlag() == kNotExist){return mGlobalPtBoundary[ptIdx];}
    return mTableMaker->GetTableData("Binning", 0, 0, 0, 1, -1, ptIdx);
}

double RHICfBinning::GetGlobalXfBinBoundary(int xfIdx)
{
    if(GetBinningTableFlag() == kNotExist){return mGlobalXfBoundary[xfIdx];}
    return mTableMaker->GetTableData("Binning", 0, 0, 0, -1, 1, xfIdx);
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
