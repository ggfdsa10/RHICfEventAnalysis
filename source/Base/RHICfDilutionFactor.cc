#include "RHICfDilutionFactor.hh"

RHICfDilutionFactor::RHICfDilutionFactor() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfDilutionFactor::~RHICfDilutionFactor()
{
}

void RHICfDilutionFactor::Init()
{
    mTableMaker -> InitTable("Dilution");
    cout << "RHICfDilutionFactor::Init() -- Done." << endl;
}

void RHICfDilutionFactor::InitDilutionData()
{
    int flag = GetDilutionTableFlag();
    for(int run=0; run<kRunNum; run++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptBinNum = mTableMaker->GetTableData("Binning", run, type, dle, 1, -1, -1)-1;
                int xfBinNum = mTableMaker->GetTableData("Binning", run, type, dle, -1, 1, -1)-1;

                if(ptBinNum > 0 && xfBinNum > 0){
                    mDilutionFactor[run][type][dle][0].resize(ptBinNum, vector<double>(xfBinNum));
                    mDilutionFactor[run][type][dle][1].resize(ptBinNum, vector<double>(xfBinNum));
                }

                if(flag == kExistTable){
                    for(int pt=0; pt<ptBinNum; pt++){
                        for(int xf=0; xf<xfBinNum; xf++){
                            mDilutionFactor[run][type][dle][0][pt][xf] = mTableMaker->GetTableData("Dilution", run, type, dle, pt, xf, 0);
                            mDilutionFactor[run][type][dle][1][pt][xf] = 0.;
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
        mDilutionFactor[run][type][dle][0][ptIdx][xfIdx] += TMath::Sin(angle*TMath::Pi()/180.);
        mDilutionFactor[run][type][dle][1][ptIdx][xfIdx] += 1.;
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
                        mDilutionFactor[run][type][dle][0][pt][xf] = mDilutionFactor[run][type][dle][0][pt][xf]/mDilutionFactor[run][type][dle][1][pt][xf];
                        if(mDilutionFactor[run][type][dle][1][pt][xf] < 10){mDilutionFactor[run][type][dle][0][pt][xf] = 0.;}
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
                        table.push_back(data);
                    }
                }       
            }
        }
    }
    mTableMaker -> SaveTable("Dilution", table);
}

int RHICfDilutionFactor::GetDilutionTableFlag()
{    
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Dilution", runIdx, 0, 3, 1, 1, 0) >= 0.){
        return kExistTable;
    }
    return kNotExist;
}

double RHICfDilutionFactor::GetDilutionFactor(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mDilutionFactor[runIdx][typeIdx][dleIdx][0][ptIdx][xfIdx];
}

TH1D* RHICfDilutionFactor::GetDilutionHist(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    return mDilutionHist[runIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}
