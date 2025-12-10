#include "RHICfDilutionFactor.hh"

RHICfDilutionFactor::RHICfDilutionFactor(TString tableName) : mTableName(tableName) 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfDilutionFactor::~RHICfDilutionFactor()
{
}

void RHICfDilutionFactor::Init()
{
    if(mTableName.Sizeof() == 1){mTableName = mOptContainer->GetTableSubName();}

    mTableMaker -> InitTable("Dilution"+mTableName);
    cout << "RHICfDilutionFactor::Init() -- Done." << endl;
}

void RHICfDilutionFactor::InitDilutionData()
{
    int flag = GetDilutionTableFlag();
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptBinNum = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, 1, -1, -1)-1;
                int xfBinNum = mTableMaker->GetTableData("Binning"+mTableName, run, type, dle, -1, 1, -1)-1;

                if(ptBinNum > 0 && xfBinNum > 0){
                    mDilutionFactor[run][type][dle][0].resize(ptBinNum, vector<double>(xfBinNum));
                    mDilutionFactor[run][type][dle][1].resize(ptBinNum, vector<double>(xfBinNum));
                }

                if(flag == kExistTable){
                    for(int pt=0; pt<ptBinNum; pt++){
                        for(int xf=0; xf<xfBinNum; xf++){
                            mDilutionFactor[run][type][dle][0][pt][xf] = mTableMaker->GetTableData("Dilution"+mTableName, run, type, dle, pt, xf, 0);
                            mDilutionFactor[run][type][dle][1][pt][xf] = mTableMaker->GetTableData("Dilution"+mTableName, run, type, dle, pt, xf, 1);
                        }
                    }
                }
            }
        }
    }
}

void RHICfDilutionFactor::InitHist()
{
    for(int run=0; run<kRunNum; run++){
        TString runName = mOptContainer -> GetRunTypeName(run);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer -> GetDLEName(dle);
            
                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                if(ptNum > 0 && xfNum > 0){
                    mDilutionHist[run][type][dle].resize(ptNum, vector<TH1D*>(xfNum));
                    for(int pt=0; pt<ptNum; pt++){
                        for(int xf=0; xf<xfNum; xf++){
                            mDilutionHist[run][type][dle][pt][xf] = new TH1D(Form("Dilution_%s_type%i_%s_pt%i_xf%i", runName.Data(), type, dleName.Data(), pt, xf), "", 180, 0., 180.);
                            mDilutionHist[run][type][dle][pt][xf] -> SetStats(0);
                            mDilutionHist[run][type][dle][pt][xf] -> SetTitle("particle #phi distribution; #phi [degree]; Counts");
                            mDilutionHist[run][type][dle][pt][xf] -> Sumw2();
                        }
                    }

                    if(mDilutionFactor[run][type][dle][0].size() == 0){
                        mDilutionFactor[run][type][dle][0].resize(ptNum, vector<double>(xfNum));
                        mDilutionFactor[run][type][dle][1].resize(ptNum, vector<double>(xfNum));
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                mDilutionFactor[run][type][dle][0][pt][xf] = 0.;
                                mDilutionFactor[run][type][dle][1][pt][xf] = 0.;
                            }
                        }
                    }
                }
            }
        }
    }
}

void RHICfDilutionFactor::SetBinning(RHICfBinning* binning)
{
    mBinning = binning;
}

void RHICfDilutionFactor::FillAngle(int run, int type, int dle, int ptIdx, int xfIdx, double angle)
{
    angle = fabs(angle);
    if(mDilutionHist[run][type][dle][ptIdx][xfIdx]){
        mDilutionHist[run][type][dle][ptIdx][xfIdx] -> Fill(angle);
    }
}

void RHICfDilutionFactor::Calculate()
{
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);

                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        int binNum = mDilutionHist[run][type][dle][pt][xf]->GetNbinsX();
                        double SinPhi = 0.;
                        double sumEntry = 0.;
                        for(int bin=0; bin<binNum; bin++){
                            double binCenter = mDilutionHist[run][type][dle][pt][xf]->GetBinCenter(bin+1);
                            double entry = mDilutionHist[run][type][dle][pt][xf]->GetBinContent(bin+1);
                            SinPhi += (TMath::Sin(binCenter*TMath::Pi()/180.)*entry);
                            sumEntry += entry;
                        }
                        SinPhi /= sumEntry;

                        double sigma = 0.;
                        double sumEntry2 = 0.;
                        for(int bin=0; bin<binNum; bin++){
                            double binCenter = mDilutionHist[run][type][dle][pt][xf]->GetBinCenter(bin+1);
                            double entry = mDilutionHist[run][type][dle][pt][xf]->GetBinContent(bin+1);
                            double sinPhi = TMath::Sin(binCenter*TMath::Pi()/180.);
                            sigma += (pow(sinPhi - SinPhi, 2.)*entry);
                            sumEntry2 += (entry*entry);
                        }
                        sigma /= (sumEntry - sumEntry2/sumEntry);
                        double SinPhiErr = sigma/sqrt(sumEntry*sumEntry/sumEntry2);
                        if(sumEntry < 10.){
                            SinPhi = 0.;
                            SinPhiErr = 0.;
                        }
                        mDilutionFactor[run][type][dle][0][pt][xf] = SinPhi;
                        mDilutionFactor[run][type][dle][1][pt][xf] = SinPhiErr;
                    }
                }
            }
        }
    }
    cout << "RHICfDilutionFactor::Calculate() -- Dilution factor calculation has done. " << endl;
}

void RHICfDilutionFactor::SaveDilutionData()
{
    vector<RHICfTableMaker::TableData> table;
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);

                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        RHICfTableMaker::TableData data;
                        data.runIdx = run;
                        data.typeIdx = type;
                        data.dleIdx = dle;
                        data.ptIdx = pt;
                        data.xfIdx = xf;
                        data.values.push_back(GetDilutionFactor(run, type, dle, pt, xf));
                        data.values.push_back(GetDilutionFactorErr(run, type, dle, pt, xf));
                        table.push_back(data);
                    }
                }       
            }
        }
    }
    mTableMaker -> SaveTable("Dilution"+mTableName, table);
}

int RHICfDilutionFactor::GetDilutionTableFlag()
{    
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Dilution"+mTableName, runIdx, 0, 3, 1, 1, 0) >= 0.){
        return kExistTable;
    }
    return kNotExist;
}

double RHICfDilutionFactor::GetDilutionFactor(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mDilutionFactor[runIdx][typeIdx][dleIdx][0][ptIdx][xfIdx];
}

double RHICfDilutionFactor::GetDilutionFactorErr(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mDilutionFactor[runIdx][typeIdx][dleIdx][1][ptIdx][xfIdx];
}

TH1D* RHICfDilutionFactor::GetDilutionHist(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mDilutionHist[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}
