#include "RHICfSystematicError.hh"

RHICfSystematicError::RHICfSystematicError() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfSystematicError::~RHICfSystematicError()
{
}

void RHICfSystematicError::Init()
{
    mBTofMult_thrBoundary[0] = 1;
    mBBCSmallEast_thrBoundary[0] = 40;
    mBBCSmallWest_thrBoundary[0] = 40;
    mBBCLargeEast_thrBoundary[0] = 6;
    mBBCLargeWest_thrBoundary[0] = 4;

    mBTofMult_thrBoundary[1] = 6;
    mBBCSmallEast_thrBoundary[1] = 300;
    mBBCSmallWest_thrBoundary[1] = 300;
    mBBCLargeEast_thrBoundary[1] = 250;
    mBBCLargeWest_thrBoundary[1] = 150;

    mBTofMult_thrBoundary[2] = 12;
    mBBCSmallEast_thrBoundary[2] = 600;
    mBBCSmallWest_thrBoundary[2] = 600;
    mBBCLargeEast_thrBoundary[2] = 500;
    mBBCLargeWest_thrBoundary[2] = 250;

    mBTofMult_thrBoundary[3] = 500;
    mBBCSmallEast_thrBoundary[3] = 100000;
    mBBCSmallWest_thrBoundary[3] = 100000;
    mBBCLargeEast_thrBoundary[3] = 100000;
    mBBCLargeWest_thrBoundary[3] = 100000;

    mTableMaker -> InitTable("SDLESystematic");
    mTableMaker -> InitTable("DDLESystematic");
    mTableMaker -> InitTable("NDLESystematic");

    cout << "RHICfSystematicError::Init() -- Done." << endl;
}

void RHICfSystematicError::InitDLESystematicData()
{
    int flag = GetDLESystematicTableFlag();
    for(int pol=0; pol<2; pol++){
        for(int fill=0; fill<kFillNum; fill++){
            int runIdx = mOptContainer->GetFillToRunIdx(fill);
            for(int type=0; type<kTypeNum; type++){
                int ptNumSDLE = mBinning -> GetPtBinNum(runIdx, type, kSDLE);
                int xfNumSDLE = mBinning -> GetXfBinNum(runIdx, type, kSDLE);
                if(ptNumSDLE > 0 && xfNumSDLE > 0){
                    if(mBeamSpinSDLE[pol][fill][type].size() == 0){
                        mBeamSpinSDLE[pol][fill][type].resize(ptNumSDLE, vector<double>(xfNumSDLE));
                        mBeamSpinSDLE[pol][fill][type].resize(ptNumSDLE, vector<double>(xfNumSDLE));
                        for(int pt=0; pt<ptNumSDLE; pt++){
                            for(int xf=0; xf<xfNumSDLE; xf++){
                                mBeamSpinSDLE[pol][fill][type][pt][xf] = 0.;

                                if(flag == kExistTable && !mOptContainer->IsForceCalculateSystematicError()){
                                    mBeamSpinSDLE[pol][fill][type][pt][xf] = mTableMaker->GetTableData("SDLESystematic", runIdx, type, fill, pt, xf, pol, 0);
                                }
                            }
                        }
                    }
                }

                int ptNumDDLE = mBinning -> GetPtBinNum(runIdx, type, kDDLE);
                int xfNumDDLE = mBinning -> GetXfBinNum(runIdx, type, kDDLE);
                if(ptNumDDLE > 0 && xfNumDDLE > 0){
                    for(int ddle=0; ddle<mDDLEBins; ddle++){
                        if(mBeamSpinDDLE[pol][ddle][fill][type].size() == 0){
                            mBeamSpinDDLE[pol][ddle][fill][type].resize(ptNumDDLE, vector<double>(xfNumDDLE));
                            mBeamSpinDDLE[pol][ddle][fill][type].resize(ptNumDDLE, vector<double>(xfNumDDLE));
                            for(int pt=0; pt<ptNumDDLE; pt++){
                                for(int xf=0; xf<xfNumDDLE; xf++){
                                    mBeamSpinDDLE[pol][ddle][fill][type][pt][xf] = 0.;

                                    if(flag == kExistTable && !mOptContainer->IsForceCalculateSystematicError()){
                                        mBeamSpinDDLE[pol][ddle][fill][type][pt][xf] = mTableMaker->GetTableData("DDLESystematic", runIdx, type, fill, pt, xf, pol, ddle);
                                    }
                                }
                            }
                        }   
                    }
                }

                int ptNumNDLE = mBinning -> GetPtBinNum(runIdx, type, kNDLE);
                int xfNumNDLE = mBinning -> GetXfBinNum(runIdx, type, kNDLE);
                if(ptNumNDLE > 0 && xfNumNDLE > 0){
                    int tmpNDLEIdx = 0;
                    for(int ndle=0; ndle<mNDLEBins; ndle++){
                        for(int ndle2=0; ndle2<mNDLEBins; ndle2++){
                            for(int ndle3=0; ndle3<mNDLEBins; ndle3++){
                                if(mBeamSpinNDLE[pol][ndle][ndle2][ndle3][fill][type].size() == 0){
                                    mBeamSpinNDLE[pol][ndle][ndle2][ndle3][fill][type].resize(ptNumNDLE, vector<double>(xfNumNDLE));
                                    mBeamSpinNDLE[pol][ndle][ndle2][ndle3][fill][type].resize(ptNumNDLE, vector<double>(xfNumNDLE));
                                    for(int pt=0; pt<ptNumNDLE; pt++){
                                        for(int xf=0; xf<xfNumNDLE; xf++){
                                            mBeamSpinNDLE[pol][ndle][ndle2][ndle3][fill][type][pt][xf] = 0.;

                                            if(flag == kExistTable && !mOptContainer->IsForceCalculateSystematicError()){
                                                mBeamSpinNDLE[pol][ndle][ndle2][ndle3][fill][type][pt][xf] = mTableMaker->GetTableData("NDLESystematic", runIdx, type, fill, pt, xf, pol, tmpNDLEIdx);
                                            }
                                        }
                                    }
                                }   
                                tmpNDLEIdx++;
                            }
                        }
                    }   
                }
            }
        }
    }
}

void RHICfSystematicError::SetBinning(RHICfBinning* binning)
{
    mBinning = binning;
}

void RHICfSystematicError::CalculateDLESystematic()
{
    InitGraph();

    // for(int fill=0; fill<kFillNum; fill++){
    //     int runIdx = mOptContainer->GetFillToRunIdx(fill);
    //     for(int type=0; type<kTypeNum; type++){
    //         int ptNumSDLE = mBinning -> GetPtBinNum(runIdx, type, kSDLE);
    //         int xfNumSDLE = mBinning -> GetXfBinNum(runIdx, type, kSDLE);
    //         for(int pt=0; pt<ptNumSDLE; pt++){
    //             for(int xf=0; xf<xfNumSDLE; xf++){
    //                 double asym = GetRawAsymmetrySDLE(fill, type, pt, xf);
    //                 mGraphAsymSDLE[fill][type][pt][xf] -> SetPoint(0, );
    //             }
    //         }
            

    //     //             double GetRawAsymmetrySDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx);
    //     // double GetRawAsymmetryDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ddleIdx);
    //     // double GetRawAsymmetryNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ndleIdx1, int ndleIdx2, int ndleIdx3);


    //         int ptNumDDLE = mBinning -> GetPtBinNum(runIdx, type, kDDLE);
    //         int xfNumDDLE = mBinning -> GetXfBinNum(runIdx, type, kDDLE);
    //         if(ptNumDDLE > 0 && xfNumDDLE > 0){
    //             for(int pt=0; pt<ptNumDDLE; pt++){
    //                 for(int xf=0; xf<xfNumDDLE; xf++){
    //                     mGraphAsymDDLE[ddle][fill][type][pt][xf];

    //                 }
    //             }
    //         }

    //         int ptNumNDLE = mBinning -> GetPtBinNum(runIdx, type, kNDLE);
    //         int xfNumNDLE = mBinning -> GetXfBinNum(runIdx, type, kNDLE);
    //         if(ptNumNDLE > 0 && xfNumNDLE > 0){
    //             for(int ndle=0; ndle<mNDLEBins; ndle++){
    //                 for(int ndle2=0; ndle2<mNDLEBins; ndle2++){
    //                     for(int ndle3=0; ndle3<mNDLEBins; ndle3++){
    //                         for(int pt=0; pt<ptNumNDLE; pt++){
    //                             for(int xf=0; xf<xfNumNDLE; xf++){
    //                                 mGraphAsymNDLE[ndle][ndle2][ndle3][fill][type][pt][xf] = new TGraphErrors();

    //                             }
    //                         }
    //                     }
    //                 }
    //             }   
    //         }
    //     }
    // }
}

void RHICfSystematicError::SaveDLESystematicData()
{
    int runType = mOptContainer -> GetRunType();
    vector<RHICfTableMaker::TableData> tableSDLE;
    vector<RHICfTableMaker::TableData> tableDDLE;
    vector<RHICfTableMaker::TableData> tableNDLE;
    for(int pol=0; pol<2; pol++){
        for(int fill=0; fill<kFillNum; fill++){
            int runIdx = mOptContainer -> GetFillToRunIdx(fill);
            if(runType != kALLRun && runType != runIdx){continue;}
            
            for(int type=0; type<kTypeNum; type++){
                int ptNumSDLE = mBinning -> GetPtBinNum(runIdx, type, kSDLE);
                int xfNumSDLE = mBinning -> GetXfBinNum(runIdx, type, kSDLE);
                for(int pt=0; pt<ptNumSDLE; pt++){
                    for(int xf=0; xf<xfNumSDLE; xf++){
                        RHICfTableMaker::TableData data;
                        data.runIdx = runIdx;
                        data.typeIdx = type;
                        data.dleIdx = fill;
                        data.ptIdx = pt;
                        data.xfIdx = xf;

                        data.values.clear();
                        data.values.push_back(pol);
                        data.values.push_back(GetPolSDLE(fill, type, pt, xf, pol));
                        tableSDLE.push_back(data);
                    }
                }

                int ptNumDDLE = mBinning -> GetPtBinNum(runIdx, type, kDDLE);
                int xfNumDDLE = mBinning -> GetXfBinNum(runIdx, type, kDDLE);
                for(int pt=0; pt<ptNumDDLE; pt++){
                    for(int xf=0; xf<xfNumDDLE; xf++){
                        RHICfTableMaker::TableData data;
                        data.runIdx = runIdx;
                        data.typeIdx = type;
                        data.dleIdx = fill;
                        data.ptIdx = pt;
                        data.xfIdx = xf;

                        data.values.clear();
                        data.values.push_back(pol);
                        for(int ddle=0; ddle<mDDLEBins; ddle++){
                            data.values.push_back(GetPolDDLE(fill, type, pt, xf, ddle, pol));
                        }
                        tableDDLE.push_back(data);
                    }
                }

                int ptNumNDLE = mBinning -> GetPtBinNum(runIdx, type, kNDLE);
                int xfNumNDLE = mBinning -> GetXfBinNum(runIdx, type, kNDLE);
                for(int pt=0; pt<ptNumNDLE; pt++){
                    for(int xf=0; xf<xfNumNDLE; xf++){
                        RHICfTableMaker::TableData data;
                        data.runIdx = runIdx;
                        data.typeIdx = type;
                        data.dleIdx = fill;
                        data.ptIdx = pt;
                        data.xfIdx = xf;

                        data.values.clear();
                        data.values.push_back(pol);
                        for(int ndle=0; ndle<mNDLEBins; ndle++){
                            for(int ndle2=0; ndle2<mNDLEBins; ndle2++){
                                for(int ndle3=0; ndle3<mNDLEBins; ndle3++){
                                    data.values.push_back(GetPolNDLE(fill, type, pt, xf, ndle, ndle2, ndle3, pol));
                                }
                            }
                        }
                        tableNDLE.push_back(data);
                    }
                }
            }
        }
    }
    mTableMaker -> SaveTable("SDLESystematic", tableSDLE);
    mTableMaker -> SaveTable("DDLESystematic", tableDDLE);
    mTableMaker -> SaveTable("NDLESystematic", tableNDLE);
}

int RHICfSystematicError::GetDLESystematicTableFlag()
{
    int runIdx = mOptContainer->GetRunType();
    if(runIdx == kALLRun){runIdx = 0;}
    if(mTableMaker -> GetTableData("SDLESystematic", runIdx, 0, 0, 0, 0, 0) >= 0){
        if(mTableMaker -> GetTableData("DDLESystematic", runIdx, 0, 0, 0, 0, 0) >= 0){
            if(mTableMaker -> GetTableData("NDLESystematic", runIdx, 0, 0, 0, 0, 0) >= 0){
                return kExistTable;
            }
        }
    }
    return kNotExist;
}

void RHICfSystematicError::FillPolSDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, bool blueSpinUp)
{
    if(blueSpinUp == true){mBeamSpinSDLE[0][fillIdx][typeIdx][ptIdx][xfIdx] += 1.;}
    if(blueSpinUp == false){mBeamSpinSDLE[1][fillIdx][typeIdx][ptIdx][xfIdx] += 1.;}
}

void RHICfSystematicError::FillPolDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int bbcSE, int bbcLE, bool blueSpinUp)
{
    for(int i=0; i<mDDLEBins; i++){
        if(mBBCSmallEast_thrBoundary[i] <= bbcSE && bbcSE < mBBCSmallEast_thrBoundary[i+1]){
            if(mBBCLargeEast_thrBoundary[i] <= bbcLE && bbcLE < mBBCLargeEast_thrBoundary[i+1]){
                if(blueSpinUp == true){mBeamSpinDDLE[0][i][fillIdx][typeIdx][ptIdx][xfIdx] += 1.;}
                if(blueSpinUp == false){mBeamSpinDDLE[1][i][fillIdx][typeIdx][ptIdx][xfIdx] += 1.;}
            }
        }
    }
}

void RHICfSystematicError::FillPolNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int btof, int bbcSE, int bbcLE, int bbcSW, int bbcLW, bool blueSpinUp)
{
    for(int i=0; i<mNDLEBins; i++){
        for(int j=0; j<mNDLEBins; j++){
            for(int k=0; k<mNDLEBins; k++){
                if(mBTofMult_thrBoundary[i] < btof && btof <= mBTofMult_thrBoundary[i+1]){
                    
                    if(mBBCSmallEast_thrBoundary[j] <= bbcSE && bbcSE < mBBCSmallEast_thrBoundary[j+1]){
                        if(mBBCLargeEast_thrBoundary[j] <= bbcLE && bbcLE < mBBCLargeEast_thrBoundary[j+1]){
                        
                            if(mBBCSmallWest_thrBoundary[k] <= bbcSW && bbcSW < mBBCSmallWest_thrBoundary[k+1]){
                                if(mBBCLargeWest_thrBoundary[k] <= bbcLW && bbcLW < mBBCLargeWest_thrBoundary[k+1]){
                                    if(blueSpinUp == true){mBeamSpinNDLE[0][i][j][k][fillIdx][typeIdx][ptIdx][xfIdx] += 1.;}
                                    if(blueSpinUp == false){mBeamSpinNDLE[1][i][j][k][fillIdx][typeIdx][ptIdx][xfIdx] += 1.;}
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void RHICfSystematicError::InitGraph()
{
    for(int fill=0; fill<kFillNum; fill++){
        int runIdx = mOptContainer->GetFillToRunIdx(fill);
        for(int type=0; type<kTypeNum; type++){
            int ptNumSDLE = mBinning -> GetPtBinNum(runIdx, type, kSDLE);
            int xfNumSDLE = mBinning -> GetXfBinNum(runIdx, type, kSDLE);
            if(xfNumSDLE > 0){
                if(mGraphAsymSDLEPt[fill][type].size() == 0){
                    mGraphAsymSDLEPt[fill][type].resize(xfNumSDLE);
                    for(int xf=0; xf<xfNumSDLE; xf++){
                        mGraphAsymSDLEPt[fill][type][xf] = new TGraphErrors();
                        mGraphAsymSDLEPt[fill][type][xf] -> SetName(Form("asymSDLE_fill%i_type%i_xf%i", fill, type, xf));
                        mGraphAsymSDLEPt[fill][type][xf] -> SetTitle("; conditions; Raw Asymmetry");
                    }
                }
            }
            if(ptNumSDLE > 0){
                if(mGraphAsymSDLEXf[fill][type].size() == 0){
                    mGraphAsymSDLEXf[fill][type].resize(ptNumSDLE);
                    for(int pt=0; pt<ptNumSDLE; pt++){
                        mGraphAsymSDLEXf[fill][type][pt] = new TGraphErrors();
                        mGraphAsymSDLEXf[fill][type][pt] -> SetName(Form("asymSDLE_fill%i_type%i_pt%i", fill, type, pt));
                        mGraphAsymSDLEXf[fill][type][pt] -> SetTitle("; conditions; Raw Asymmetry");
                    }
                }
            }

            

            int ptNumDDLE = mBinning -> GetPtBinNum(runIdx, type, kDDLE);
            int xfNumDDLE = mBinning -> GetXfBinNum(runIdx, type, kDDLE);
            for(int ddle=0; ddle<mDDLEBins; ddle++){
                if(xfNumDDLE > 0){
                    if(mGraphAsymDDLEPt[ddle][fill][type].size() == 0){
                        mGraphAsymDDLEPt[ddle][fill][type].resize(xfNumDDLE);
                        for(int xf=0; xf<xfNumDDLE; xf++){
                            mGraphAsymDDLEPt[ddle][fill][type][xf] = new TGraphErrors();
                            mGraphAsymDDLEPt[ddle][fill][type][xf] -> SetName(Form("asymDDLE%i_fill%s_type%i_xf%i", ddle, fill, type, xf));
                            mGraphAsymDDLEPt[ddle][fill][type][xf] -> SetTitle("; conditions; Raw Asymmetry");
                        }
                    }
                }   
                if(ptNumDDLE > 0){
                    if(mGraphAsymDDLEXf[ddle][fill][type].size() == 0){
                        mGraphAsymDDLEXf[ddle][fill][type].resize(ptNumDDLE);
                        for(int pt=0; pt<ptNumDDLE; pt++){
                            mGraphAsymDDLEXf[ddle][fill][type][pt] = new TGraphErrors();
                            mGraphAsymDDLEXf[ddle][fill][type][pt] -> SetName(Form("asymDDLE%i_fill%s_type%i_pt%i", ddle, fill, type, pt));
                            mGraphAsymDDLEXf[ddle][fill][type][pt] -> SetTitle("; conditions; Raw Asymmetry");
                        }
                    }
                }   
            }


            // int ptNumNDLE = mBinning -> GetPtBinNum(runIdx, type, kNDLE);
            // int xfNumNDLE = mBinning -> GetXfBinNum(runIdx, type, kNDLE);

            // for(int ndle=0; ndle<mNDLEBins; ndle++){
            //     for(int ndle2=0; ndle2<mNDLEBins; ndle2++){
            //         for(int ndle3=0; ndle3<mNDLEBins; ndle3++){
            //             if(xfNumDDLE > 0){
            //             if(mGraphAsymNDLE[ndle][ndle2][ndle3][fill][type].size() == 0){
            //                 mGraphAsymNDLE[ndle][ndle2][ndle3][fill][type].resize(xfNumNDLE);
            //                 for(int pt=0; pt<ptNumNDLE; pt++){
            //                     for(int xf=0; xf<xfNumNDLE; xf++){
            //                         mGraphAsymNDLE[ndle][ndle2][ndle3][fill][type][pt][xf] = new TGraphErrors();
            //                         mGraphAsymNDLE[ndle][ndle2][ndle3][fill][type][pt][xf] -> SetName(Form("asymDDLE%i_%i_%i_fill%s_type%i_pt%i_xf%i", ndle, ndle2, ndle3, fill, type, pt, xf));
            //                         mGraphAsymNDLE[ndle][ndle2][ndle3][fill][type][pt][xf] -> SetTitle("; conditions; Raw Asymmetry");
            //                     }
            //                 }
            //             }   
            //         }
            //     }
            // }   
        
        }
    }
}

double RHICfSystematicError::GetPolSDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int polIdx)
{   
    return mBeamSpinSDLE[polIdx][fillIdx][typeIdx][ptIdx][xfIdx];
}

double RHICfSystematicError::GetPolDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ddleIdx, int polIdx)
{   
    return mBeamSpinDDLE[polIdx][ddleIdx][fillIdx][typeIdx][ptIdx][xfIdx];
}

double RHICfSystematicError::GetPolNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ndleIdx1, int ndleIdx2, int ndleIdx3, int polIdx)
{   
    return mBeamSpinNDLE[polIdx][ndleIdx1][ndleIdx2][ndleIdx3][fillIdx][typeIdx][ptIdx][xfIdx];
}

double RHICfSystematicError::GetRawAsymmetrySDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx)
{
    double upPolNum = GetPolSDLE(fillIdx, typeIdx, ptIdx, xfIdx, 0);
    double downPolNum = GetPolSDLE(fillIdx, typeIdx, ptIdx, xfIdx, 1);
    double luminosity = RelativeLuminosity(fillIdx);

    double asymmetry = (luminosity*downPolNum - upPolNum)/(luminosity*downPolNum + upPolNum);
    return asymmetry;
}

double RHICfSystematicError::GetRawAsymmetryDDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ddleIdx)
{
    double upPolNum = GetPolDDLE(fillIdx, typeIdx, ptIdx, xfIdx, ddleIdx, 0);
    double downPolNum = GetPolDDLE(fillIdx, typeIdx, ptIdx, xfIdx, ddleIdx, 1);
    double luminosity = RelativeLuminosity(fillIdx);

    double asymmetry = (luminosity*downPolNum - upPolNum)/(luminosity*downPolNum + upPolNum);
    return asymmetry;
}

double RHICfSystematicError::GetRawAsymmetryNDLE(int fillIdx, int typeIdx, int ptIdx, int xfIdx, int ndleIdx1, int ndleIdx2, int ndleIdx3)
{
    double upPolNum = GetPolNDLE(fillIdx, typeIdx, ptIdx, xfIdx, ndleIdx1, ndleIdx2, ndleIdx3, 0);
    double downPolNum = GetPolNDLE(fillIdx, typeIdx, ptIdx, xfIdx, ndleIdx1, ndleIdx2, ndleIdx3, 1);
    double luminosity = RelativeLuminosity(fillIdx);

    double asymmetry = (luminosity*downPolNum - upPolNum)/(luminosity*downPolNum + upPolNum);
    return asymmetry;
}

double RHICfSystematicError::RelativeLuminosity(int fillIdx)
{
    if(fillIdx == 0){return 0.9581;}
    if(fillIdx == 1){return 0.9623;}
    if(fillIdx == 2){return 0.9924;}
    if(fillIdx == 3){return 0.9949;}
    if(fillIdx == 4){return 0.9774;}
    return 0.;
}
