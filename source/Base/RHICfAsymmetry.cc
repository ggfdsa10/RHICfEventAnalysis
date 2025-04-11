#include "RHICfAsymmetry.hh"

RHICfAsymmetry::RHICfAsymmetry() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfAsymmetry::~RHICfAsymmetry()
{
}

void RHICfAsymmetry::Init()
{
    mTableMaker -> InitTable("Asymmetry");
    cout << "RHICfAsymmetry::Init() -- Done." << endl;
}

void RHICfAsymmetry::InitAsymmetryData()
{
    int flag = GetAsymmetryTableFlag();
    for(int fill=0; fill<kFillNum; fill++){
        int runIdx = mOptContainer->GetFillToRunIdx(fill);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptNum = mTableMaker->GetTableData("Binning", runIdx, type, dle, 1, -1, -1)-1;
                int xfNum = mTableMaker->GetTableData("Binning", runIdx, type, dle, -1, 1, -1)-1;

                if(ptNum > 0 && xfNum > 0){
                    mRightAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mLeftAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mBkgRightAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mBkgLeftAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));   

                    for(int pt=0; pt<ptNum; pt++){
                        for(int xf=0; xf<xfNum; xf++){
                            mRightAsymmetry[fill][type][dle][pt][xf] = 0.;
                            mLeftAsymmetry[fill][type][dle][pt][xf] = 0.;
                            mBkgRightAsymmetry[fill][type][dle][pt][xf] = 0.;
                            mBkgLeftAsymmetry[fill][type][dle][pt][xf] = 0.;

                            if(flag == kExistTable && !mOptContainer->IsForceCalculateAsymmetry()){
                                mRightAsymmetry[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Asymmetry", runIdx, type, dle, pt, xf, fill, 0);
                                mLeftAsymmetry[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Asymmetry", runIdx, type, dle, pt, xf, fill, 1);
                                mBkgRightAsymmetry[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Asymmetry", runIdx, type, dle, pt, xf, fill, 2);
                                mBkgLeftAsymmetry[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Asymmetry", runIdx, type, dle, pt, xf, fill, 3);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RHICfAsymmetry::SetBinning(RHICfBinning* binning)
{
    mBinning = binning;
}

void RHICfAsymmetry::SetMassFitting(RHICfMassFitting* massFitting)
{
    mMassFitting = massFitting;
}

void RHICfAsymmetry::SetDilution(RHICfDilutionFactor* dilution)
{
    mDilution = dilution;
}

void RHICfAsymmetry::FillPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double angle, bool isSpinUp)
{
    // int runIdx = mOptContainer->GetFillToRunIdx(fillIdx);
    // if(isSpinUp == true){
    //     if(runIdx != kTLRun){
    //         if(angle > 0.){mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // right
    //         if(angle < 0.){mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //     }
    //     if(runIdx == kTLRun){
    //         if(angle > 0.){mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //         if(angle < 0.){mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // right
    //     }
    // }
    // if(isSpinUp == false){
    //     if(runIdx != kTLRun){
    //         if(angle > 0.){mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //         if(angle < 0.){mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // right
    //     }
    //     if(runIdx == kTLRun){
    //         if(angle > 0.){mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;}// right
    //         if(angle < 0.){mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //     }
    // }

    if(0. < angle && angle < 180.){
        if(isSpinUp == true){
            mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
        else{
            mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
    }
    else{
        if(isSpinUp == true){
            mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
        else{
            mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
    }
}

void RHICfAsymmetry::FillBkgPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double angle, bool isSpinUp)
{
    // int runIdx = mOptContainer->GetFillToRunIdx(fillIdx);
    // if(isSpinUp == true){
    //     if(runIdx != kTLRun){
    //         if(angle > 0.){mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // right
    //         if(angle < 0.){mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //     }
    //     if(runIdx == kTLRun){
    //         if(angle > 0.){mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //         if(angle < 0.){mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // right
    //     }
    // }
    // if(isSpinUp == false){
    //     if(runIdx != kTLRun){
    //         if(angle > 0.){mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //         if(angle < 0.){mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // right
    //     }
    //     if(runIdx == kTLRun){
    //         if(angle > 0.){mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;}// right
    //         if(angle < 0.){mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;} // left
    //     }
    // }

    if(0. < angle && angle < 180.){
        if(isSpinUp == true){
            mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
        else{
            mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
    }
    else{
        if(isSpinUp == true){
            mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
        else{
            mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
        }
    }
}

void RHICfAsymmetry::CalculateAN()
{
    cout << "RHICfAsymmetry::CalculateAN() -- start.." << endl;
    InitGraph();


    int particleRunIdx = mOptContainer -> GetParticleRunIdx();
    if(particleRunIdx == kPi0Run){AsymmetryPi0();}
    if(particleRunIdx == kNeutronRun){AsymmetryNeutron();}

    cout << "RHICfAsymmetry::CalculateAN() -- Done." << endl;
}

void RHICfAsymmetry::SaveAsymmetryData()
{
    int runType = mOptContainer -> GetRunType();
    vector<RHICfTableMaker::TableData> table;
    for(int fill=0; fill<kFillNum; fill++){
        int runIdx = mOptContainer -> GetFillToRunIdx(fill);

        if(runType != kALLRun && runType != runIdx){continue;}
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptNum = mBinning -> GetPtBinNum(runIdx, type, dle);
                int xfNum = mBinning -> GetXfBinNum(runIdx, type, dle);

                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        RHICfTableMaker::TableData data;
                        data.runIdx = mOptContainer->GetFillToRunIdx(fill);
                        data.typeIdx = type;
                        data.dleIdx = dle;
                        data.ptIdx = pt;
                        data.xfIdx = xf;

                        data.values.clear();
                        data.values.push_back(fill);
                        data.values.push_back(GetPolNum(fill, type, dle, pt, xf, true));
                        data.values.push_back(GetPolNum(fill, type, dle, pt, xf, false));
                        data.values.push_back(GetBkgPolNum(fill, type, dle, pt, xf, true));
                        data.values.push_back(GetBkgPolNum(fill, type, dle, pt, xf, false));

                        table.push_back(data);
                    }
                }       
            }
        }
    }
    mTableMaker -> SaveTable("Asymmetry", table);
}

int RHICfAsymmetry::GetAsymmetryTableFlag()
{
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Asymmetry", runIdx, 0, 0, 0, 0, 0) >= 0){
        return kExistTable;
    }
    return kNotExist;
}

double RHICfAsymmetry::GetPolNum(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand)
{
    if(isRightHand){return mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];}
    return mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

double RHICfAsymmetry::GetBkgPolNum(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand)
{
    if(isRightHand){return mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];}
    return mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

TGraphErrors* RHICfAsymmetry::GetANGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx)
{
    if(isPtGraph){return mGraphAN_Global_pT[runIdx][typeIdx][dleIdx][anType][binIdx];}
    return mGraphAN_Global_xF[runIdx][typeIdx][dleIdx][anType][binIdx];
}

TGraphErrors* RHICfAsymmetry::GetANGraph(bool isPtGraph, int anType, int dleIdx, int binIdx)
{
    if(isPtGraph){return mGraphAN_pT[dleIdx][anType][binIdx];}
    return mGraphAN_xF[dleIdx][anType][binIdx];
}

TGraphErrors* RHICfAsymmetry::GetANSummaryGraph(bool isPtGraph, int anType, int dleIdx)
{
    if(isPtGraph){return mGraphAN_pTSummary[dleIdx][anType];}
    return mGraphAN_xFSummary[dleIdx][anType];
}


void RHICfAsymmetry::AsymmetryPi0()
{
    int globalPtNum = mBinning -> GetGlobalPtBinNum();
    int globalXfNum = mBinning -> GetGlobalXfBinNum();
    for(int run=0; run<kRunNum; run++){
        // if(run == kTLRun){continue;}
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){

                for(int globalPt=0; globalPt<globalPtNum; globalPt++){
                    for(int globalXf=0; globalXf<globalXfNum; globalXf++){
                        double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                        double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                        double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                        double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                    if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){
                                        for(int i=0; i<3; i++){
                                            double sumAN = 0.;
                                            double sumANStatErr = 0.;
                                            double sumStatistic = 0.;

                                            for(int fill=0; fill<kFillNum; fill++){
                                                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                                if(run != runIdx){continue;}

                                                double an = -999.;
                                                double statErr = -999.;
                                                double counts = 0.;
                                                if(i == 0){
                                                    an = GetInclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 1){
                                                    an = GetBkgAN(fill, type, dle, pt, xf);
                                                    statErr = GetBkgANStatError(fill, type, dle, pt, xf);
                                                    counts = GetBkgANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 2){
                                                    an = GetExclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetExclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetExclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(an <= -900. || statErr <= -900. || counts <= 10.){continue;}
                                                
                                                sumAN += (an*counts);
                                                sumANStatErr += (counts/(statErr*statErr));
                                                sumStatistic += counts;
                                            }

                                            sumAN /= sumStatistic; 
                                            sumANStatErr = 1./sqrt(sumANStatErr/sumStatistic);

                                            int pTGraphNum = mGraphAN_Global_pT[run][type][dle][i][globalXf]->GetN();
                                            mGraphAN_Global_pT[run][type][dle][i][globalXf] -> SetPoint(pTGraphNum, ptMean, sumAN);
                                            mGraphAN_Global_pT[run][type][dle][i][globalXf] -> SetPointError(pTGraphNum, 0., sumANStatErr);

                                            int xFGraphNum = mGraphAN_Global_xF[run][type][dle][i][globalPt]->GetN();
                                            mGraphAN_Global_xF[run][type][dle][i][globalPt] -> SetPoint(xFGraphNum, xfMean, sumAN);
                                            mGraphAN_Global_xF[run][type][dle][i][globalPt] -> SetPointError(xFGraphNum, 0., sumANStatErr);
                                        }
                                    }
                                }
                            }
                        }           
                    }
                }
            }
        }
    }


    for(int dle=0; dle<kDLENum; dle++){
        for(int globalPt=0; globalPt<globalPtNum; globalPt++){
            for(int globalXf=0; globalXf<globalXfNum; globalXf++){
                double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                for(int i=0; i<3; i++){
                    double sumAN = 0.;
                    double sumPt = 0.;
                    double sumXf = 0.;
                    double sumANStatErr = 0.;
                    double sumStatistic = 0.;
                    for(int run=0; run<kRunNum; run++){
                        // if(run == kTLRun){continue;}
                        for(int type=0; type<kTypeNum; type++){
                            int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                            int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                            for(int pt=0; pt<ptNum; pt++){
                                for(int xf=0; xf<xfNum; xf++){
                                    double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                    double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                    if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                        if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){
                    
                                            for(int fill=0; fill<kFillNum; fill++){
                                                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                                if(run != runIdx){continue;}

                                                double an = -999.;
                                                double statErr = -999.;
                                                double counts = 0.;
                                                if(i == 0){
                                                    an = GetInclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 1){
                                                    an = GetBkgAN(fill, type, dle, pt, xf);
                                                    statErr = GetBkgANStatError(fill, type, dle, pt, xf);
                                                    counts = GetBkgANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 2){
                                                    an = GetExclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetExclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetExclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(an <= -900. || statErr <= -900. || counts <= 10.){continue;}
                                                
                                                sumAN += (an*counts);
                                                sumPt += (ptMean*counts);
                                                sumXf += (xfMean*counts);
                                                sumANStatErr += (1./(statErr*statErr));
                                                sumStatistic += counts;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    sumAN /= sumStatistic;
                    sumPt /= sumStatistic;
                    sumXf /= sumStatistic;
                    sumANStatErr = sqrt(1./(sumANStatErr));

                    int pTGraphNum = mGraphAN_pT[dle][i][globalXf]->GetN();
                    mGraphAN_pT[dle][i][globalXf] -> SetPoint(pTGraphNum, sumPt, sumAN);
                    mGraphAN_pT[dle][i][globalXf] -> SetPointError(pTGraphNum, 0., sumANStatErr);

                    int xFGraphNum = mGraphAN_xF[dle][i][globalPt]->GetN();
                    mGraphAN_xF[dle][i][globalPt] -> SetPoint(xFGraphNum, sumXf, sumAN);
                    mGraphAN_xF[dle][i][globalPt] -> SetPointError(xFGraphNum, 0., sumANStatErr);
                }
            }
        }
    }

    for(int dle=0; dle<kDLENum; dle++){
        for(int globalPt=0; globalPt<globalPtNum; globalPt++){

            for(int i=0; i<3; i++){
                double sumAN = 0.;
                double sumPt = 0.;
                double sumANStatErr = 0.;
                double sumStatistic = 0.;
                for(int globalXf=0; globalXf<globalXfNum; globalXf++){
                    double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                    double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                    double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                    double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                    for(int run=0; run<kRunNum; run++){
                        // if(run == kTLRun){continue;}
                        for(int type=0; type<kTypeNum; type++){
                            int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                            int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                            for(int pt=0; pt<ptNum; pt++){
                                for(int xf=0; xf<xfNum; xf++){
                                    double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                    double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                    if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                        if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){
                                            for(int fill=0; fill<kFillNum; fill++){
                                                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                                if(run != runIdx){continue;}

                                                double an = -999.;
                                                double statErr = -999.;
                                                double counts = 0.;
                                                if(i == 0){
                                                    an = GetInclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 1){
                                                    an = GetBkgAN(fill, type, dle, pt, xf);
                                                    statErr = GetBkgANStatError(fill, type, dle, pt, xf);
                                                    counts = GetBkgANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 2){
                                                    an = GetExclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetExclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetExclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(an <= -900. || statErr <= -900. || counts <= 10.){continue;}

                                                sumAN += (an*counts);
                                                sumPt += (ptMean*counts);
                                                sumANStatErr += (1./(statErr*statErr));
                                                sumStatistic += counts;
                                            }
                                        
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                sumAN /= sumStatistic;
                sumPt /= sumStatistic;
                sumANStatErr = sqrt(1./(sumANStatErr));

                int pTGraphNum = mGraphAN_pTSummary[dle][i]->GetN();
                mGraphAN_pTSummary[dle][i] -> SetPoint(pTGraphNum, sumPt, sumAN);
                mGraphAN_pTSummary[dle][i] -> SetPointError(pTGraphNum, 0., sumANStatErr);
            }
        }


        for(int globalXf=0; globalXf<globalXfNum; globalXf++){

            for(int i=0; i<3; i++){
                double sumAN = 0.;
                double sumXf = 0.;
                double sumANStatErr = 0.;
                double sumStatistic = 0.;
                for(int globalPt=0; globalPt<globalPtNum; globalPt++){
                    double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                    double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                    double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                    double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                    for(int run=0; run<kRunNum; run++){
                        // if(run == kTLRun){continue;}
                        for(int type=0; type<kTypeNum; type++){
                            int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                            int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                            for(int pt=0; pt<ptNum; pt++){
                                for(int xf=0; xf<xfNum; xf++){
                                    double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                    double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                    if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                        if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){
                                            for(int fill=0; fill<kFillNum; fill++){
                                                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                                if(run != runIdx){continue;}

                                                double an = -999.;
                                                double statErr = -999.;
                                                double counts = 0.;
                                                if(i == 0){
                                                    an = GetInclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 1){
                                                    an = GetBkgAN(fill, type, dle, pt, xf);
                                                    statErr = GetBkgANStatError(fill, type, dle, pt, xf);
                                                    counts = GetBkgANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(i == 2){
                                                    an = GetExclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetExclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetExclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(an <= -900. || statErr <= -900. || counts <= 10.){continue;}
                                                
                                                sumAN += (an*counts);
                                                sumXf += (xfMean*counts);
                                                sumANStatErr += (1./(statErr*statErr));
                                                sumStatistic += counts;
                                            }
                                        
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                sumAN /= sumStatistic;
                sumXf /= sumStatistic;
                sumANStatErr = sqrt(1./(sumANStatErr));

                int GraphNum = mGraphAN_xFSummary[dle][i]->GetN();
                mGraphAN_xFSummary[dle][i] -> SetPoint(GraphNum, sumXf, sumAN);
                mGraphAN_xFSummary[dle][i] -> SetPointError(GraphNum, 0., sumANStatErr);
            }
        }

    }

}

void RHICfAsymmetry::AsymmetryNeutron()
{
    int globalPtNum = mBinning -> GetGlobalPtBinNum();
    int globalXfNum = mBinning -> GetGlobalXfBinNum();
    for(int run=0; run<kRunNum; run++){
        
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                for(int globalPt=0; globalPt<globalPtNum; globalPt++){
                    for(int globalXf=0; globalXf<globalXfNum; globalXf++){
                        double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                        double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                        double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                        double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                    if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){

                                        double sumAN = 0.;
                                        double sumANStatErr = 0.;
                                        double sumStatistic = 0.;

                                        for(int fill=0; fill<kFillNum; fill++){
                                            int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                            if(run != runIdx){continue;}

                                            double an = GetInclusiveAN(fill, type, dle, pt, xf);
                                            double statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                            double counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                            if(an <= -900. || statErr <= -900. || counts <= 200.){continue;}
                                            
                                            sumAN += (an*counts);
                                            sumANStatErr += (1./(statErr*statErr));
                                            sumStatistic += counts;
                                        }

                                        sumAN /= sumStatistic;
                                        sumANStatErr = sqrt(1./(sumANStatErr));

                                        int pTGraphNum = mGraphAN_Global_pT[run][type][dle][0][globalXf]->GetN();
                                        mGraphAN_Global_pT[run][type][dle][0][globalXf] -> SetPoint(pTGraphNum, ptMean, sumAN);
                                        mGraphAN_Global_pT[run][type][dle][0][globalXf] -> SetPointError(pTGraphNum, 0., sumANStatErr);

                                        int xFGraphNum = mGraphAN_Global_xF[run][type][dle][0][globalPt]->GetN();
                                        mGraphAN_Global_xF[run][type][dle][0][globalPt] -> SetPoint(xFGraphNum, xfMean, sumAN);
                                        mGraphAN_Global_xF[run][type][dle][0][globalPt] -> SetPointError(xFGraphNum, 0., sumANStatErr);
                                    
                                    }
                                }
                            }
                        }           
                    }
                }
            }
        }
    }


    for(int dle=0; dle<kDLENum; dle++){
        for(int globalPt=0; globalPt<globalPtNum; globalPt++){
            for(int globalXf=0; globalXf<globalXfNum; globalXf++){
                double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                double sumAN = 0.;
                double sumPt = 0.;
                double sumXf = 0.;
                double sumANStatErr = 0.;
                double sumStatistic = 0.;
                for(int run=0; run<kRunNum; run++){
                    
                    for(int type=1; type<kTypeNum; type++){
                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                    if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){

                                        for(int fill=0; fill<kFillNum; fill++){
                                            int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                            if(run != runIdx){continue;}

                                            double an = GetInclusiveAN(fill, type, dle, pt, xf);
                                            double statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                            double counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                            if(an <= -900. || statErr <= -900. || counts <= 200.){continue;}
                                            
                                            sumAN += (an*counts);
                                            sumPt += (ptMean*counts);
                                            sumXf += (xfMean*counts);
                                            sumANStatErr += (1./(statErr*statErr));
                                            sumStatistic += counts;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                sumAN /= sumStatistic;
                sumPt /= sumStatistic;
                sumXf /= sumStatistic;
                sumANStatErr = sqrt(1./(sumANStatErr));
                
                int pTGraphNum = mGraphAN_pT[dle][0][globalXf]->GetN();
                mGraphAN_pT[dle][0][globalXf] -> SetPoint(pTGraphNum, sumPt, sumAN);
                mGraphAN_pT[dle][0][globalXf] -> SetPointError(pTGraphNum, 0., sumANStatErr);

                int xFGraphNum = mGraphAN_xF[dle][0][globalPt]->GetN();
                mGraphAN_xF[dle][0][globalPt] -> SetPoint(xFGraphNum, sumXf, sumAN);
                mGraphAN_xF[dle][0][globalPt] -> SetPointError(xFGraphNum, 0., sumANStatErr);
            }
        }
    }



    for(int dle=0; dle<kDLENum; dle++){
        for(int globalPt=0; globalPt<globalPtNum; globalPt++){
            double sumAN = 0.;
            double sumPt = 0.;
            double sumANStatErr = 0.;
            double sumStatistic = 0.;
            for(int globalXf=0; globalXf<globalXfNum; globalXf++){
                double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                for(int run=0; run<kRunNum; run++){
                    
                    for(int type=0; type<kTypeNum; type++){
                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                    if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){
                                        for(int i=0; i<1; i++){
                                            for(int fill=0; fill<kFillNum; fill++){
                                                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                                if(run != runIdx){continue;}

                                                double an = -999.;
                                                double statErr = -999.;
                                                double counts = 0.;
                                                if(i == 0){
                                                    an = GetInclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(an <= -900. || statErr <= -900. || counts <= 200.){continue;}
                                                
                                                sumAN += (an*counts);
                                                sumPt += (ptMean*counts);
                                                sumANStatErr += (1./(statErr*statErr));
                                                sumStatistic += counts;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            sumAN /= sumStatistic;
            sumPt /= sumStatistic;
            sumANStatErr = sqrt(1./(sumANStatErr));

            int pTGraphNum = mGraphAN_pTSummary[dle][0]->GetN();
            mGraphAN_pTSummary[dle][0] -> SetPoint(pTGraphNum, sumPt, sumAN);
            mGraphAN_pTSummary[dle][0] -> SetPointError(pTGraphNum, 0., sumANStatErr);
        }


        for(int globalXf=0; globalXf<globalXfNum; globalXf++){
            double sumAN = 0.;
            double sumXf = 0.;
            double sumANStatErr = 0.;
            double sumStatistic = 0.;
            for(int globalPt=0; globalPt<globalPtNum; globalPt++){
                double globalPtLowerBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt);
                double globalPtUpperBoundary = mBinning -> GetGlobalPtBinBoundary(globalPt+1);
                double globalXfLowerBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf);
                double globalXfUpperBoundary = mBinning -> GetGlobalXfBinBoundary(globalXf+1);

                for(int run=0; run<kRunNum; run++){
                    
                    for(int type=0; type<kTypeNum; type++){
                        int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                        int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                        for(int pt=0; pt<ptNum; pt++){
                            for(int xf=0; xf<xfNum; xf++){
                                double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                                double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

                                if(globalPtLowerBoundary <= ptMean && ptMean < globalPtUpperBoundary){
                                    if(globalXfLowerBoundary <= xfMean && xfMean < globalXfUpperBoundary){
                                        for(int i=0; i<1; i++){
                                            for(int fill=0; fill<kFillNum; fill++){
                                                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                                if(run != runIdx){continue;}

                                                double an = -999.;
                                                double statErr = -999.;
                                                double counts = 0.;
                                                if(i == 0){
                                                    an = GetInclusiveAN(fill, type, dle, pt, xf);
                                                    statErr = GetInclusiveANStatError(fill, type, dle, pt, xf);
                                                    counts = GetInclusiveANCounts(fill, type, dle, pt, xf);   
                                                }
                                                if(an <= -999. || statErr <= -999. || counts <= 200.){continue;}
        
                                                sumAN += (an*counts);
                                                sumXf += (xfMean*counts);
                                                sumANStatErr += (1./(statErr*statErr));
                                                sumStatistic += counts;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            sumAN /= sumStatistic;
            sumXf /= sumStatistic;
            sumANStatErr = sqrt(1./(sumANStatErr));

            int GraphNum = mGraphAN_xFSummary[dle][0]->GetN();
            mGraphAN_xFSummary[dle][0] -> SetPoint(GraphNum, sumXf, sumAN);
            mGraphAN_xFSummary[dle][0] -> SetPointError(GraphNum, 0., sumANStatErr);
        }

    }

}

double RHICfAsymmetry::GetInclusiveAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double dilutionFactor = mDilution -> GetDilutionFactor(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(dilutionFactor <= 0.1){return -999.;}

    double beamPol = BeamPolarization(fillIdx);
    double luminosity = RelativeLuminosity(fillIdx);

    double rightAllNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double leftAllNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(rightAllNum < 10 || leftAllNum < 10){return -999.;}

    double allAN = -1./(dilutionFactor*beamPol) * (leftAllNum - luminosity*rightAllNum)/(leftAllNum + luminosity*rightAllNum);
    return allAN;
}

double RHICfAsymmetry::GetInclusiveANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double rightAllNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double leftAllNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(rightAllNum < 10 || leftAllNum < 10){return 0.;}
    return rightAllNum+leftAllNum;
}

double RHICfAsymmetry::GetInclusiveANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double counts = GetInclusiveANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(counts <= 10.){return -999.;}
    return 1./sqrt(counts);
}

double RHICfAsymmetry::GetBkgAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double dilutionFactor = mDilution -> GetDilutionFactor(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(dilutionFactor <= 0.1){return -999.;}

    double beamPol = BeamPolarization(fillIdx);
    double luminosity = RelativeLuminosity(fillIdx);

    double rightBkgNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double leftBkgNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(rightBkgNum < 1 || leftBkgNum < 1){return -999.;}

    double bkgAN = -1./(dilutionFactor*beamPol) * (leftBkgNum - luminosity*rightBkgNum)/(leftBkgNum + luminosity*rightBkgNum);
    return bkgAN;
}

double RHICfAsymmetry::GetBkgANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double rightBkgNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double leftBkgNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(rightBkgNum < 1 || leftBkgNum < 1){return 0.;}
    return rightBkgNum+leftBkgNum;
}

double RHICfAsymmetry::GetBkgANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double counts = GetBkgANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    return 1./sqrt(counts);
}

double RHICfAsymmetry::GetExclusiveAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double inclusiveAN = GetInclusiveAN(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgAN = GetBkgAN(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);

    double inclusiveANCounts = GetInclusiveANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgANCounts = GetBkgANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(inclusiveANCounts < 0 || bkgANCounts < 0){return -999.;}

    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double allMassCounts = mMassFitting -> GetMassAllCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double signalMassCounts = mMassFitting -> GetMassSignalCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgMassCounts = mMassFitting -> GetMassBkgCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(allMassCounts <= 0. || signalMassCounts <= 0. || bkgMassCounts <= 0.){return -999.;}
    if(signalMassCounts > allMassCounts){return -999.;}

    double exclusiveAN = (allMassCounts/signalMassCounts) * inclusiveAN - (bkgMassCounts/signalMassCounts) * bkgAN;
    return exclusiveAN;
}

double RHICfAsymmetry::GetExclusiveANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double inclusiveANCounts = GetInclusiveANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgANCounts = GetBkgANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(inclusiveANCounts < 1 || bkgANCounts < 1){return 0;}

    double allMassCounts = mMassFitting -> GetMassAllCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double signalMassCounts = mMassFitting -> GetMassSignalCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgMassCounts = mMassFitting -> GetMassBkgCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(allMassCounts <= 0. || signalMassCounts <= 0. || bkgMassCounts <= 0.){return -999.;}
    if(signalMassCounts > allMassCounts){return -999.;}

    double exclusiveAN = (signalMassCounts/allMassCounts) * inclusiveANCounts;
    return exclusiveAN;
}

double RHICfAsymmetry::GetExclusiveANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double inclusiveANError = GetInclusiveANStatError(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgANError = GetBkgANStatError(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);

    double allMassCounts = mMassFitting -> GetMassAllCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double signalMassCounts = mMassFitting -> GetMassSignalCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgMassCounts = mMassFitting -> GetMassBkgCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(allMassCounts <= 0. || signalMassCounts <= 0. || bkgMassCounts <= 0.){return -999.;}
    if(signalMassCounts > allMassCounts){return -999.;}
    
    double exclusiveError = (allMassCounts/signalMassCounts)*(allMassCounts/signalMassCounts) * (inclusiveANError*inclusiveANError) + (bkgMassCounts/signalMassCounts)*(bkgMassCounts/signalMassCounts) * (bkgANError*bkgANError);
    return sqrt(exclusiveError);
}

bool RHICfAsymmetry::RejectPoint(int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double allCount = mMassFitting -> GetMassAllCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double sigCount = mMassFitting -> GetMassSignalCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double bkgCount = mMassFitting -> GetMassBkgCounts(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double fitChi2 = mMassFitting -> GetMassFitChi2(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double sigma = mMassFitting -> GetMassSigma(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double peakDiff= mMassFitting -> GetMassPeakDifference(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    double fitPeak = mMassFitting -> GetMassFitPeak(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);

    if(allCount <= 1.){return true;}
    if(sigCount <= 1.){return true;}
    if(bkgCount <= 1.){return true;}
    // if(fitChi2 >= 400.){return true;}
    if(sigma >= 14.){return true;}
    // if(peakDiff/fitPeak >= 9.){return true;}
    return false;
}


void RHICfAsymmetry::InitGraph()
{
    int ptNum = mBinning -> GetGlobalPtBinNum();
    int xfNum = mBinning -> GetGlobalXfBinNum();
    for(int dle=0; dle<kDLENum; dle++){
        TString dleName = mOptContainer -> GetDLEName(dle);
        for(int i=0; i<3; i++){

            mGraphAN_pT[dle][i].resize(xfNum);
            for(int xf=0; xf<xfNum; xf++){
                mGraphAN_pT[dle][i][xf] = new TGraphErrors();
                mGraphAN_pT[dle][i][xf] -> SetName(Form("AN_pT_%s_xfBin%i", dleName.Data(), xf));
                mGraphAN_pT[dle][i][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
            }
            
            mGraphAN_xF[dle][i].resize(ptNum);
            for(int pt=0; pt<ptNum; pt++){
                mGraphAN_xF[dle][i][pt] = new TGraphErrors();
                mGraphAN_xF[dle][i][pt] -> SetName(Form("AN_xF_%s_ptBin%i", dleName.Data(), pt));
                mGraphAN_xF[dle][i][pt] -> SetTitle("; x_{F}; A_{N}");
            }

            for(int run=0; run<kRunNum; run++){
                TString runName = mOptContainer -> GetRunTypeName(run);
                for(int type=0; type<kTypeNum; type++){

                    mGraphAN_Global_pT[run][type][dle][i].resize(xfNum);
                    for(int xf=0; xf<xfNum; xf++){
                        mGraphAN_Global_pT[run][type][dle][i][xf] = new TGraphErrors();
                        mGraphAN_Global_pT[run][type][dle][i][xf] -> SetName(Form("AN_pT_%s_type%i_%s_xfBin%i", runName.Data(), type, dleName.Data(), xf));
                        mGraphAN_Global_pT[run][type][dle][i][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
                    }

                    mGraphAN_Global_xF[run][type][dle][i].resize(ptNum);
                    for(int pt=0; pt<ptNum; pt++){
                        mGraphAN_Global_xF[run][type][dle][i][pt] = new TGraphErrors();
                        mGraphAN_Global_xF[run][type][dle][i][pt] -> SetName(Form("AN_xF_%s_type%i_%s_ptBin%i", runName.Data(), type, dleName.Data(), pt));
                        mGraphAN_Global_xF[run][type][dle][i][pt] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
                    }
                }
            }
            mGraphAN_pTSummary[dle][i] = new TGraphErrors();
            mGraphAN_pTSummary[dle][i]  -> SetName(Form("AN_pTSummary_%s_%i", dleName.Data(), i));
            mGraphAN_pTSummary[dle][i]  -> SetTitle("; p_{T} [GeV/c]; A_{N}");

            mGraphAN_xFSummary[dle][i] = new TGraphErrors();
            mGraphAN_xFSummary[dle][i]  -> SetName(Form("AN_xFSummary_%s_%i", dleName.Data(), i));
            mGraphAN_xFSummary[dle][i]  -> SetTitle("; x_{F}; A_{N}");
        }
    }
}

double RHICfAsymmetry::RelativeLuminosity(int fillIdx)
{
    if(fillIdx == 0){return 0.9581;}
    if(fillIdx == 1){return 0.9623;}
    if(fillIdx == 2){return 0.9924;}
    if(fillIdx == 3){return 0.9949;}
    if(fillIdx == 4){return 0.9774;}
    return 0.;
}

double RHICfAsymmetry::BeamPolarization(int fillIdx)
{
    if(fillIdx = 0){return 0.536;}
    if(fillIdx = 1){return 0.554;}
    if(fillIdx = 2){return 0.590;}
    if(fillIdx = 3){return 0.566;}
    if(fillIdx = 4){return 0.592;}
    return 0.;
}

