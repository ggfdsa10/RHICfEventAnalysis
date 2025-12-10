#include "RHICfPolarization.hh"

RHICfPolarization::RHICfPolarization(TString tableName) : mTableName(tableName) 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfPolarization::~RHICfPolarization()
{
}

void RHICfPolarization::Init()
{
    if(mTableName.Sizeof() == 1){mTableName = mOptContainer->GetTableSubName();}
    mTableMaker -> InitTable("Polarization"+mTableName);
    cout << "RHICfPolarization::Init() -- Done." << endl;
}

void RHICfPolarization::InitPolarizationData()
{
    int flag = GetPolarizationTableFlag();
    for(int fill=0; fill<kFillNum; fill++){
        int runIdx = mOptContainer->GetFillToRunIdx(fill);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int ptNum = mTableMaker->GetTableData("Binning"+mTableName, runIdx, type, dle, 1, -1, -1)-1;
                int xfNum = mTableMaker->GetTableData("Binning"+mTableName, runIdx, type, dle, -1, 1, -1)-1;

                if(ptNum > 0 && xfNum > 0){
                    mUpCounts[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mDownCounts[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mBkgUpCounts[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mBkgDownCounts[fill][type][dle].resize(ptNum, vector<double>(xfNum));   

                    for(int pt=0; pt<ptNum; pt++){
                        for(int xf=0; xf<xfNum; xf++){
                            mUpCounts[fill][type][dle][pt][xf] = 0.;
                            mDownCounts[fill][type][dle][pt][xf] = 0.;
                            mBkgUpCounts[fill][type][dle][pt][xf] = 0.;
                            mBkgDownCounts[fill][type][dle][pt][xf] = 0.;

                            if(flag == kExistTable && !mOptContainer->IsForceCalculatePolarization()){
                                mUpCounts[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Polarization"+mTableName, runIdx, type, dle, pt, xf, fill, 0);
                                mDownCounts[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Polarization"+mTableName, runIdx, type, dle, pt, xf, fill, 1);
                                mBkgUpCounts[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Polarization"+mTableName, runIdx, type, dle, pt, xf, fill, 2);
                                mBkgDownCounts[fill][type][dle][pt][xf] = mTableMaker->GetTableData("Polarization"+mTableName, runIdx, type, dle, pt, xf, fill, 3);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RHICfPolarization::SetBinning(RHICfBinning* binning)
{
    mBinning = binning;
}

void RHICfPolarization::SetMassFitting(RHICfMassFitting* massFitting)
{
    mMassFitting = massFitting;
}

void RHICfPolarization::SetDilution(RHICfDilutionFactor* dilution)
{
    mDilution = dilution;
}

void RHICfPolarization::FillPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double angle, bool isSpinUp)
{
    int runIdx = mOptContainer->GetFillToRunIdx(fillIdx);
    if(runIdx != kTLRun){
        if(0. < angle){
            if(isSpinUp == true){
                mUpCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
            else{
                mDownCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
        }
    }
    else{
        if(0. > angle){
            if(isSpinUp == true){
                mUpCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
            else{
                mDownCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
        }
    }
}

void RHICfPolarization::FillBkgPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, double angle, bool isSpinUp)
{
    int runIdx = mOptContainer->GetFillToRunIdx(fillIdx);
    if(runIdx != kTLRun){
        if(0. < angle){
            if(isSpinUp == true){
                mBkgUpCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
            else{
                mBkgDownCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
        }
    }
    else{
        if(0. > angle){
            if(isSpinUp == true){
                mBkgUpCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
            else{
                mBkgDownCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
            }
        }
    }
}

void RHICfPolarization::CalculateAN()
{
    cout << "RHICfPolarization::CalculateAN() -- start.." << endl;
    InitGraph();


    int particleRunIdx = mOptContainer -> GetParticleRunIdx();
    // if(particleRunIdx == kPi0Run){PolarizationPi0();}
    // if(particleRunIdx == kNeutronRun){PolarizationNeutron();}

    cout << "RHICfPolarization::CalculateAN() -- Done." << endl;
}

void RHICfPolarization::SavePolarizationData()
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
    mTableMaker -> SaveTable("Polarization"+mTableName, table);
}

int RHICfPolarization::GetPolarizationTableFlag()
{
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("Polarization"+mTableName, runIdx, 0, 0, 0, 0, 0) >= 0){
        return kExistTable;
    }
    return kNotExist;
}

double RHICfPolarization::GetPolNum(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isUP)
{
    if(isUP){return mUpCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];}
    return mDownCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

double RHICfPolarization::GetBkgPolNum(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isUP)
{
    if(isUP){return mBkgUpCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];}
    return mBkgDownCounts[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx];
}

// TGraphErrors* RHICfPolarization::GetANGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx)
// {
//     if(isPtGraph){return mGraphAN_Global_pT[runIdx][typeIdx][dleIdx][anType][binIdx];}
//     return mGraphAN_Global_xF[runIdx][typeIdx][dleIdx][anType][binIdx];
// }

// TGraphErrors* RHICfPolarization::GetANGraph(bool isPtGraph, int anType, int dleIdx, int binIdx)
// {
//     if(isPtGraph){return mGraphAN_pT[dleIdx][anType][binIdx];}
//     return mGraphAN_xF[dleIdx][anType][binIdx];
// }

// TGraphErrors* RHICfPolarization::GetANSummaryGraph(bool isPtGraph, int anType, int dleIdx)
// {
//     if(isPtGraph){return mGraphAN_pTSummary[dleIdx][anType];}
//     return mGraphAN_xFSummary[dleIdx][anType];
// }


// void RHICfPolarization::PolarizationPi0()
// {
//     const int ptNum = mBinning -> GetGlobalPtBinNum();
//     const int xfNum = mBinning -> GetGlobalXfBinNum();

//     for(int run=0; run<kRunNum; run++){
//         for(int type=0; type<kTypeNum; type++){
//             for(int dle=0; dle<kDLENum; dle++){
//                 for(int pt=0; pt<ptNum; pt++){
//                     for(int xf=0; xf<xfNum; xf++){
//                         double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
//                         double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);

//                         double dilutionFactor = mDilution -> GetDilutionFactor(run, type, dle, pt, xf);
//                         if(dilutionFactor <= 0.1){continue;}

//                         double anWeight = 0.;
//                         double errWeight = 0.;
//                         for(int fill=0; fill<kFillNum; fill++){
//                             int runIdx = mOptContainer -> GetFillToRunIdx(fill);
//                             if(run != runIdx){continue;}
//                             double pol = BeamPolarization(fill);
//                             double luminosity = RelativeLuminosity(fill);
//                             double up = GetPolNum(fill, type, dle, pt, xf, true);
//                             double down = GetPolNum(fill, type, dle, pt, xf, false);

//                             double an = -1./(pol*dilutionFactor) * (up - luminosity*down)/(up + luminosity*down);
//                             if(run == kTLRun){an = -an;}

//                             double err = 2.*luminosity*sqrt(up*down*(up+down))/(pol*dilutionFactor*pow((up + luminosity*down), 2.));
                            
//                             double w = 1./(err*err);
//                             anWeight += (an*w);
//                             errWeight += w;
//                         }
                        
//                         anWeight /= errWeight;
//                         errWeight = 1./sqrt(errWeight);

//                         int pTGraphNum = mGraphAN_Global_pT[run][type][dle][0][xf]->GetN();
//                         mGraphAN_Global_pT[run][type][dle][0][xf] -> SetPoint(pTGraphNum, ptMean, anWeight);
//                         mGraphAN_Global_pT[run][type][dle][0][xf] -> SetPointError(pTGraphNum, 0., errWeight);

//                         int xFGraphNum = mGraphAN_Global_xF[run][type][dle][0][pt]->GetN();
//                         mGraphAN_Global_xF[run][type][dle][0][pt] -> SetPoint(xFGraphNum, xfMean, anWeight);
//                         mGraphAN_Global_xF[run][type][dle][0][pt] -> SetPointError(xFGraphNum, 0., errWeight);
//                     }
//                 }
//             }
//         }
//     }

// }

double RHICfPolarization::GetInclusiveAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double dilutionFactor = mDilution -> GetDilutionFactor(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    // double dilutionFactor = 1.;
    if(dilutionFactor <= 0.1){return -999.;}

    double beamPol = BeamPolarization(fillIdx);
    double luminosity = RelativeLuminosity(fillIdx);

    double upNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double downNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(upNum < 10 || downNum < 10){return -999.;}

    double an = -1./(dilutionFactor*beamPol) * (upNum - luminosity*downNum)/(upNum + luminosity*downNum);
    if(runIdx == kTLRun){an = -an;}
    return an;
}

double RHICfPolarization::GetInclusiveANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double upNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double downNum = GetPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(upNum < 10 || downNum < 10){return 0.;}
    return upNum+downNum;
}

double RHICfPolarization::GetInclusiveANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double counts = GetInclusiveANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    if(counts <= 10.){return -999.;}
    return 1./sqrt(counts);
}

double RHICfPolarization::GetBkgAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    int runIdx = mOptContainer -> GetFillToRunIdx(fillIdx);
    double dilutionFactor = mDilution -> GetDilutionFactor(runIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    // double dilutionFactor = 1.;
    if(dilutionFactor <= 0.1){return -999.;}

    double beamPol = BeamPolarization(fillIdx);
    double luminosity = RelativeLuminosity(fillIdx);

    double upBkgNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double downBkgNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(upBkgNum < 1 || downBkgNum < 1){return -999.;}

    double an = -1./(dilutionFactor*beamPol) * (upBkgNum - luminosity*downBkgNum)/(upBkgNum + luminosity*downBkgNum);
    if(runIdx == kTLRun){an = -an;}
    return an;
}

double RHICfPolarization::GetBkgANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double upNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, true);
    double downNum = GetBkgPolNum(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx, false);
    if(upNum < 1 || downNum < 1){return 0.;}
    return upNum+downNum;
}

double RHICfPolarization::GetBkgANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
{
    double counts = GetBkgANCounts(fillIdx, typeIdx, dleIdx, ptIdx, xfIdx);
    return 1./sqrt(counts);
}

double RHICfPolarization::GetExclusiveAN(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
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

double RHICfPolarization::GetExclusiveANCounts(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
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

double RHICfPolarization::GetExclusiveANStatError(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx)
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

void RHICfPolarization::InitGraph()
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

double RHICfPolarization::RelativeLuminosity(int fillIdx)
{
    if(fillIdx == 0){return 0.9581;}
    if(fillIdx == 1){return 0.9623;}
    if(fillIdx == 2){return 0.9924;}
    if(fillIdx == 3){return 0.9949;}
    if(fillIdx == 4){return 0.9774;}
    return 0.;
}

double RHICfPolarization::BeamPolarization(int fillIdx)
{
    if(fillIdx = 0){return 0.536;}
    if(fillIdx = 1){return 0.554;}
    if(fillIdx = 2){return 0.590;}
    if(fillIdx = 3){return 0.566;}
    if(fillIdx = 4){return 0.592;}
    return 0.;
}

