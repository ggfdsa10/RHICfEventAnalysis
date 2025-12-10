#include "RHICfAsymmetry.hh"

RHICfAsymmetry::RHICfAsymmetry()
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mTableMaker = RHICfTableMaker::GetTableMaker();
}

RHICfAsymmetry::~RHICfAsymmetry()
{
}

int RHICfAsymmetry::Init()
{
    RHICfOptContainer::SetOffTPCTrack();
    RHICfOptContainer::SetOffBTof();
    RHICfOptContainer::SetOffVPD();
    RHICfOptContainer::SetOffZDC();
    RHICfOptContainer::SetOffFMS();

    GetOptContainer()->SetConditionName("");
    GetOptContainer()->CalculatePi0(); 

    if(!RHICfOptContainer::Init()){return 0;}
    if(!RHICfTableMaker::Init()){return 0;}

    for(int run=0; run<kRunNum; run++){
        TString runName = GetOptContainer()->GetRunTypeName(run);
        GetOptContainer()->SetRunType(runName);
        int beamTOPRefNum = (run == kTOPRun)? kBeamTOPRefNum : 1;

        for(int i=0; i<kBeamMetNum; i++){
            for(int j=0; j<beamTOPRefNum; j++){
                GetOptContainer()->SetBeamCenterMethod(i+1, j+1);
                TString tableSubName = GetOptContainer()->GetTableSubName();
            
                mMassFitting[run][i][j] = new RHICfMassFitting(tableSubName);
                mBinning[run][i][j] = new RHICfBinning(tableSubName);
                mDilution[run][i][j] = new RHICfDilutionFactor(tableSubName);
                mPolarization[run][i][j] = new RHICfPolarization(tableSubName);

                mMassFitting[run][i][j] -> Init();
                mBinning[run][i][j] -> Init();
                mDilution[run][i][j] -> Init();
                mPolarization[run][i][j] -> Init();

                mMassFitting[run][i][j] -> InitMassData();
                mBinning[run][i][j] -> InitBinningData();
                mDilution[run][i][j] -> InitDilutionData();
                mPolarization[run][i][j] -> InitPolarizationData();
            }
        }
    }

    cout << "RHICfAsymmetry::Init() -- Done." << endl;
}

void RHICfAsymmetry::Calculate()
{
    cout << "RHICfAsymmetry::Calculate() -- start.." << endl;

    InitGraph();
    SystematicErrorCalculatePi0();
    AsymmetryPi0();

    DrawAsymmetryGraph();

    cout << "RHICfAsymmetry::Calculate() -- Done." << endl;
}

TGraphErrors* RHICfAsymmetry::GetANGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx)
{
    if(isPtGraph){return mGraphAN_pT[runIdx][typeIdx][dleIdx][anType][binIdx];}
    return mGraphAN_xF[runIdx][typeIdx][dleIdx][anType][binIdx];
}

TGraphErrors* RHICfAsymmetry::GetANBkgGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx)
{
    if(isPtGraph){return mGraphAN_Bkg_pT[runIdx][typeIdx][dleIdx][anType][binIdx];}
    return mGraphAN_Bkg_xF[runIdx][typeIdx][dleIdx][anType][binIdx];
}

TGraphErrors* RHICfAsymmetry::GetANSubtGraph(bool isPtGraph, int anType, int runIdx, int typeIdx, int dleIdx, int binIdx)
{
    if(isPtGraph){return mGraphAN_Subt_pT[runIdx][typeIdx][dleIdx][anType][binIdx];}
    return mGraphAN_Subt_xF[runIdx][typeIdx][dleIdx][anType][binIdx];
}

void RHICfAsymmetry::AsymmetryPi0()
{
    int ptNum = 0;
    int xfNum = 0;
    for(int run=0; run<kRunNum; run++){
        TString runName = GetOptContainer()->GetRunTypeName(run);
        GetOptContainer()->SetRunType(runName);
        if(GetOptContainer()->GetRunType() != run && GetOptContainer()->GetRunType() != kALLRun){continue;}
        int beamRefNum = (run == kTOPRun)? kBeamTOPRefNum : 1;

        if(mBinning[run][0][0] -> GetBinningTableFlag() == kNotExist){continue;}
        ptNum = mBinning[run][0][0] -> GetGlobalPtBinNum();
        xfNum = mBinning[run][0][0] -> GetGlobalXfBinNum();

        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){

                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        if(mBinning[run][0][0] -> GetBinningTableFlag() == kNotExist){continue;}
                        double ptMean = mBinning[run][0][0] -> GetPtBinMean(run, type, dle, pt, xf);
                        double xfMean = mBinning[run][0][0] -> GetXfBinMean(run, type, dle, pt, xf);

                        double sumAN = 0.;
                        double sumErr = 0.;

                        double sumBkgAN = 0.;
                        double sumBkgErr = 0.;

                        double sumSubtAN = 0.;
                        double sumSubtErr = 0.;

                        int graphIdx = -1;
                        for(int i=0; i<kBeamMetNum; i++){
                            for(int j=0; j<beamRefNum; j++){
                                graphIdx++;
                                if(mMassFitting[run][i][j] -> GetMassTableFlag() == kNotExist){continue;}
                                if(mDilution[run][i][j] -> GetDilutionTableFlag() == kNotExist){continue;}
                                if(mPolarization[run][i][j] -> GetPolarizationTableFlag() == kNotExist){continue;}

                                double dilutionFactor = mDilution[run][i][j] -> GetDilutionFactor(run, type, dle, pt, xf);
                                if(dilutionFactor <= 0.1){continue;}
                                double allCounts = mMassFitting[run][i][j] -> GetMassAllCounts(false, run, type, dle, pt, xf);
                                double signalCounts = mMassFitting[run][i][j] -> GetMassSignalCounts(false, run, type, dle, pt, xf);
                                double bkgCounts = mMassFitting[run][i][j] -> GetMassBkgCounts(false, run, type, dle, pt, xf);

                                double anWeight = 0.;
                                double errWeight = 0.;

                                double anBkgWeight = 0.;
                                double errBkgWeight = 0.;

                                double anSubtWeight = 0.;
                                double errSubtWeight = 0.;

                                int tmp = 0;
                                for(int fill=0; fill<kFillNum; fill++){
                                    int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                    if(run != runIdx){continue;}
                                    double pol = BeamPolarization(fill);
                                    double luminosity = RelativeLuminosity(fill);

                                    double up = mPolarization[run][i][j] -> GetPolNum(fill, type, dle, pt, xf, true);
                                    double down = mPolarization[run][i][j] -> GetPolNum(fill, type, dle, pt, xf, false);
                                    if(up <= 10 || down <= 10){continue;}

                                    double an = -1./(pol*dilutionFactor) * (up - luminosity*down)/(up + luminosity*down);
                                    if(run == kTLRun){an = -an;}

                                    double err = 2.*luminosity*sqrt(up*down*(up+down))/(pol*dilutionFactor*pow((up + luminosity*down), 2.));
                                    double w = 1./(err*err);
                                    anWeight += (an*w);
                                    errWeight += w;

                                    double upBkg = mPolarization[run][i][j] -> GetBkgPolNum(fill, type, dle, pt, xf, true);
                                    double downBkg = mPolarization[run][i][j] -> GetBkgPolNum(fill, type, dle, pt, xf, false);

                                    double anBkg = -1./(pol*dilutionFactor) * (upBkg - luminosity*downBkg)/(upBkg + luminosity*downBkg);
                                    if(run == kTLRun){anBkg = -anBkg;}

                                    double errBkg = 2.*luminosity*sqrt(upBkg*downBkg*(upBkg+downBkg))/(pol*dilutionFactor*pow((upBkg + luminosity*downBkg), 2.));
                                    double wBkg = 1./(errBkg*errBkg);
                                    anBkgWeight += (anBkg*wBkg);
                                    errBkgWeight += wBkg;

                                    double massMean = mMassFitting[run][i][j] -> GetMassMean(run, type, dle, pt, xf);
                                    if(massMean <= 10.){continue;}
                                    double BSratio = mMassFitting[run][i][j] -> GetBSRatio_GPR(run, type, dle, pt, xf);
                                    double ASratio = 1. + BSratio;

                                    if(upBkg <= 1 || downBkg <= 1){continue;}

                                    double subtAN = ASratio*an - BSratio*anBkg;
                                    double subtErr = sqrt(ASratio*ASratio*err*err + BSratio*BSratio*errBkg*errBkg);

                                    double wSubt = 1./(subtErr*subtErr);
                                    anSubtWeight += (subtAN*wSubt);
                                    errSubtWeight += wSubt;
                                
                                }

                                if(errWeight > 0.1){
                                    anWeight /= errWeight;
                                    errWeight = 1./sqrt(errWeight);
                                    if(errWeight > 0.1){continue;}

                                    sumAN += anWeight;
                                    sumErr += (errWeight*errWeight);

                                    int pTGraphNum = mGraphAN_pT[run][type][dle][graphIdx][xf]->GetN();
                                    mGraphAN_pT[run][type][dle][graphIdx][xf] -> SetPoint(pTGraphNum, ptMean, anWeight);
                                    mGraphAN_pT[run][type][dle][graphIdx][xf] -> SetPointError(pTGraphNum, 0., errWeight);

                                    int xFGraphNum = mGraphAN_xF[run][type][dle][graphIdx][pt]->GetN();
                                    mGraphAN_xF[run][type][dle][graphIdx][pt] -> SetPoint(xFGraphNum, xfMean, anWeight);
                                    mGraphAN_xF[run][type][dle][graphIdx][pt] -> SetPointError(xFGraphNum, 0., errWeight);
                                }
                                if(errBkgWeight > 0.1){
                                    anBkgWeight /= errBkgWeight;
                                    errBkgWeight = 1./sqrt(errBkgWeight);

                                    sumBkgAN += anBkgWeight;
                                    sumBkgErr += (errBkgWeight*errBkgWeight);
                                    
                                    int pTBkgGraphNum = mGraphAN_Bkg_pT[run][type][dle][graphIdx][xf]->GetN();
                                    mGraphAN_Bkg_pT[run][type][dle][graphIdx][xf] -> SetPoint(pTBkgGraphNum, ptMean, anBkgWeight);
                                    mGraphAN_Bkg_pT[run][type][dle][graphIdx][xf] -> SetPointError(pTBkgGraphNum, 0., errBkgWeight);

                                    int xFBkgGraphNum = mGraphAN_Bkg_xF[run][type][dle][graphIdx][pt]->GetN();
                                    mGraphAN_Bkg_xF[run][type][dle][graphIdx][pt] -> SetPoint(xFBkgGraphNum, xfMean, anBkgWeight);
                                    mGraphAN_Bkg_xF[run][type][dle][graphIdx][pt] -> SetPointError(xFBkgGraphNum, 0., errBkgWeight);
                                }
                                if(errSubtWeight > 0.1){
                                    anSubtWeight /= errSubtWeight;
                                    errSubtWeight = 1./sqrt(errSubtWeight);
                                    if(errSubtWeight > 0.1){continue;}
                                    sumSubtAN += anSubtWeight;
                                    sumSubtErr += (errSubtWeight*errSubtWeight);
                                    
                                    int pTSubtGraphNum = mGraphAN_Subt_pT[run][type][dle][graphIdx][xf]->GetN();
                                    mGraphAN_Subt_pT[run][type][dle][graphIdx][xf] -> SetPoint(pTSubtGraphNum, ptMean, anSubtWeight);
                                    mGraphAN_Subt_pT[run][type][dle][graphIdx][xf] -> SetPointError(pTSubtGraphNum, 0., errSubtWeight);

                                    int xFSubtGraphNum = mGraphAN_Subt_xF[run][type][dle][graphIdx][pt]->GetN();
                                    mGraphAN_Subt_xF[run][type][dle][graphIdx][pt] -> SetPoint(xFSubtGraphNum, xfMean, anSubtWeight);
                                    mGraphAN_Subt_xF[run][type][dle][graphIdx][pt] -> SetPointError(xFSubtGraphNum, 0., errSubtWeight);
                                }
                            }
                        }

                        if(sumErr > 0.000001){
                            sumAN /= double(graphIdx+1);
                            sumErr = sqrt(sumErr)/double(graphIdx+1);
                            int pTGraphNum = mGraphAN_pT[run][type][dle][4][xf]->GetN();
                            mGraphAN_pT[run][type][dle][4][xf] -> SetPoint(pTGraphNum, ptMean, sumAN);
                            mGraphAN_pT[run][type][dle][4][xf] -> SetPointError(pTGraphNum, 0., sumErr);

                            int xFGraphNum = mGraphAN_xF[run][type][dle][4][pt]->GetN();
                            mGraphAN_xF[run][type][dle][4][pt] -> SetPoint(xFGraphNum, xfMean, sumAN);
                            mGraphAN_xF[run][type][dle][4][pt] -> SetPointError(xFGraphNum, 0., sumErr);
                        }
                        if(sumBkgErr > 0.000001){
                            sumBkgAN /= double(graphIdx+1);
                            sumBkgErr = sqrt(sumBkgErr)/double(graphIdx+1);
                            int pTBkgGraphNum = mGraphAN_Bkg_pT[run][type][dle][4][xf]->GetN();
                            mGraphAN_Bkg_pT[run][type][dle][4][xf] -> SetPoint(pTBkgGraphNum, ptMean, sumBkgAN);
                            mGraphAN_Bkg_pT[run][type][dle][4][xf] -> SetPointError(pTBkgGraphNum, 0., sumBkgErr);

                            int xFBkgGraphNum = mGraphAN_Bkg_xF[run][type][dle][4][pt]->GetN();
                            mGraphAN_Bkg_xF[run][type][dle][4][pt] -> SetPoint(xFBkgGraphNum, xfMean, sumBkgAN);
                            mGraphAN_Bkg_xF[run][type][dle][4][pt] -> SetPointError(xFBkgGraphNum, 0., sumBkgErr);
                        }
                        if(sumSubtErr > 0.000001){
                            sumSubtAN /= double(graphIdx+1);
                            sumSubtErr = sqrt(sumSubtErr)/double(graphIdx+1);

                            int pTSubtGraphNum = mGraphAN_Subt_pT[run][type][dle][4][xf]->GetN();
                            mGraphAN_Subt_pT[run][type][dle][4][xf] -> SetPoint(pTSubtGraphNum, ptMean, sumSubtAN);
                            mGraphAN_Subt_pT[run][type][dle][4][xf] -> SetPointError(pTSubtGraphNum, 0., sumSubtErr);

                            int xFSubtGraphNum = mGraphAN_Subt_xF[run][type][dle][4][pt]->GetN();
                            mGraphAN_Subt_xF[run][type][dle][4][pt] -> SetPoint(xFSubtGraphNum, xfMean, sumSubtAN);
                            mGraphAN_Subt_xF[run][type][dle][4][pt] -> SetPointError(xFSubtGraphNum, 0., sumSubtErr);

                            double sysErrPt = 0.;
                            double sysErrXf = 0.;
                            for(int s=0; s<4; s++){
                                if(s==1){continue;}
                                double minPointPt = 999.;
                                double maxPointPt = -999.;
                                double minPointXf = 999.;
                                double maxPointXf = -999.;

                                if(s < 3){
                                    for(int k=0; k<2; k++){
                                        for(int met=0; met<4; met++){
                                            double Ptx, Pty; 
                                            mGraphANSysErr_Calcul_subt_pT[run][type][dle][met][s][k][xf] -> GetPoint(pTSubtGraphNum, Ptx, Pty);
                                            if(minPointPt > Pty){
                                                minPointPt = Pty;
                                            }
                                            if(maxPointPt < Pty){
                                                maxPointPt = Pty;
                                            }

                                            double Xfx, Xfy; 
                                            mGraphANSysErr_Calcul_subt_xF[run][type][dle][met][s][k][pt] -> GetPoint(xFSubtGraphNum, Xfx, Xfy);
                                            if(minPointXf > Xfy){
                                                minPointXf = Xfy;
                                            }
                                            if(maxPointXf < Xfy){
                                                maxPointXf = Xfy;
                                            }
                                        }
                                    }
                                }

                                if(s == 3){
                                    for(int g=0; g<4; g++){
                                        int pTGraphNum = mGraphAN_Subt_pT[run][type][dle][g][xf]->GetN();
                                        for(int p=0; p<pTGraphNum; p++){
                                            double x, y;
                                            mGraphAN_Subt_pT[run][type][dle][g][xf] -> GetPoint(p, x, y);
                                            if(fabs(x - ptMean) < 0.0001){
                                                if(minPointPt > y){
                                                    minPointPt = y;
                                                }
                                                if(maxPointPt < y){
                                                    maxPointPt = y;
                                                }
                                            }
                                        }

                                        int xFGraphNum = mGraphAN_Subt_xF[run][type][dle][g][pt]->GetN();
                                        for(int p=0; p<xFGraphNum; p++){
                                            double x, y;
                                            mGraphAN_Subt_xF[run][type][dle][g][pt] -> GetPoint(p, x, y);
                                            if(fabs(x - xfMean) < 0.0001){
                                                if(minPointXf > y){
                                                    minPointXf = y;
                                                }
                                                if(maxPointXf < y){
                                                    maxPointXf = y;
                                                }
                                            }
                                        }
                                    }
                                }
                                double errPt = fabs(maxPointPt - minPointPt)/2.;
                                sysErrPt += (errPt*errPt);

                                double errXf = fabs(maxPointXf - minPointXf)/2.;
                                sysErrXf += (errXf*errXf);

                                mGraphANSysErr_Subt_pT[run][type][dle][s][xf] -> SetPoint(pTSubtGraphNum, ptMean, sumSubtAN); 
                                mGraphANSysErr_Subt_pT[run][type][dle][s][xf] -> SetPointError(pTSubtGraphNum, 0.01, errPt);
                                mGraphANSysErr_Subt_xF[run][type][dle][s][pt] -> SetPoint(xFSubtGraphNum, xfMean, sumSubtAN);
                                mGraphANSysErr_Subt_xF[run][type][dle][s][pt] -> SetPointError(xFSubtGraphNum, 0.01, errXf);
                            }
                            sysErrPt = sqrt(sysErrPt);
                            sysErrXf = sqrt(sysErrXf);

                            mGraphANSysErr_Subt_pT[run][type][dle][4][xf] -> SetPoint(pTSubtGraphNum, ptMean, sumSubtAN); 
                            mGraphANSysErr_Subt_pT[run][type][dle][4][xf] -> SetPointError(pTSubtGraphNum, 0.01, sysErrPt);
                            mGraphANSysErr_Subt_xF[run][type][dle][4][pt] -> SetPoint(xFSubtGraphNum, xfMean, sumSubtAN);
                            mGraphANSysErr_Subt_xF[run][type][dle][4][pt] -> SetPointError(xFSubtGraphNum, 0.01, sysErrXf);
                        }
                    }
                }
            }
        }
    }

    for(int dle=0; dle<kDLENum; dle++){
        for(int xf=0; xf<xfNum; xf++){
            for(int pt=0; pt<ptNum; pt++){
                double ptBoundL = mBinning[0][0][0] -> GetGlobalPtBinBoundary(pt);
                double ptBoundU = mBinning[0][0][0] -> GetGlobalPtBinBoundary(pt+1);
                double xfBoundL = mBinning[0][0][0] -> GetGlobalXfBinBoundary(xf);
                double xfBoundU = mBinning[0][0][0] -> GetGlobalXfBinBoundary(xf+1);

                double sumPt = 0.;
                double sumAN = 0.;
                double sumErr = 0.;
                double sumSysErr[5];
                memset(sumSysErr, 0., sizeof(sumSysErr));

                for(int run=0; run<kRunNum; run++){
                    for(int type=0; type<kTypeNum; type++){
                        int gNum = mGraphAN_Subt_pT[run][type][dle][4][xf]->GetN();
                        for(int p=0; p<gNum; p++){
                            double pt, an;
                            mGraphAN_Subt_pT[run][type][dle][4][xf] -> GetPoint(p, pt, an);
                            double err = mGraphAN_Subt_pT[run][type][dle][4][xf] -> GetErrorY(p);
                            if(ptBoundL <= pt && pt < ptBoundU){
                                double w = 1./(err*err);

                                sumPt += (pt*w);
                                sumAN += (an*w);
                                sumErr += w;

                                for(int s=0; s<5; s++){
                                    if(s==1){continue;}
                                    double sysErr = mGraphANSysErr_Subt_pT[run][type][dle][s][xf] -> GetErrorY(p);
                                    sumSysErr[s] += (sysErr*sysErr*w*w);
                                }
                            }
                        }
                    }
                }

                if(sumErr > 0.1){
                    sumPt /= sumErr;
                    sumAN /= sumErr;
                    double statErr = 1./sqrt(sumErr);

                    int pTGraphNum = mGraphAN_pTSummary[dle][xf]->GetN();
                    mGraphAN_pTSummary[dle][xf] -> SetPoint(pTGraphNum, sumPt, sumAN);
                    mGraphAN_pTSummary[dle][xf] -> SetPointError(pTGraphNum, 0., statErr);

                    for(int s=0; s<5; s++){
                        if(s==1){continue;}
                        double sysErr = sqrt(sumSysErr[s])/sumErr;
                        mGraphAN_pTSummarySysErr[dle][s][xf] -> SetPoint(pTGraphNum, sumPt, sumAN);
                        mGraphAN_pTSummarySysErr[dle][s][xf] -> SetPointError(pTGraphNum, 0.01, sysErr);
                    }
                }

                double sumXf = 0.;
                sumAN = 0.;
                sumErr = 0.;
                memset(sumSysErr, 0., sizeof(sumSysErr));

                for(int run=0; run<kRunNum; run++){
                    for(int type=0; type<kTypeNum; type++){
                        int gNum = mGraphAN_Subt_xF[run][type][dle][4][pt]->GetN();
                        for(int p=0; p<gNum; p++){
                            double xf, an;
                            mGraphAN_Subt_xF[run][type][dle][4][pt] -> GetPoint(p, xf, an);
                            double err = mGraphAN_Subt_xF[run][type][dle][4][pt] -> GetErrorY(p);
                            if(xfBoundL <= xf && xf < xfBoundU){
                                double w = 1./(err*err);

                                sumXf += (xf*w);
                                sumAN += (an*w);
                                sumErr += w;

                                for(int s=0; s<5; s++){
                                    if(s==1){continue;}
                                    double sysErr = mGraphANSysErr_Subt_xF[run][type][dle][s][pt] -> GetErrorY(p);
                                    sumSysErr[s] += (sysErr*sysErr*w*w);
                                }
                            }
                        }
                    }
                }

                if(sumErr > 0.1){
                    sumXf /= sumErr;
                    sumAN /= sumErr;
                    double statErr = 1./sqrt(sumErr);

                    int xFGraphNum = mGraphAN_xFSummary[dle][pt]->GetN();
                    mGraphAN_xFSummary[dle][pt] -> SetPoint(xFGraphNum, sumXf, sumAN);
                    mGraphAN_xFSummary[dle][pt] -> SetPointError(xFGraphNum, 0., statErr);

                    for(int s=0; s<5; s++){
                        if(s==1){continue;}
                        double sysErr = sqrt(sumSysErr[s])/sumErr;
                        mGraphAN_xFSummarySysErr[dle][s][pt] -> SetPoint(xFGraphNum, sumXf, sumAN);
                        mGraphAN_xFSummarySysErr[dle][s][pt] -> SetPointError(xFGraphNum, 0.01, sysErr);
                    }
                }
            }
        }
    }



    for(int dle=0; dle<kDLENum; dle++){
        for(int i=0; i<2; i++){
            for(int pt=0; pt<4; pt++){
                double ptBoundL = mBinning[0][0][0] -> GetGlobalPtBinBoundary(pt);
                double ptBoundU = mBinning[0][0][0] -> GetGlobalPtBinBoundary(pt+1);

                double sumPt = 0.;
                double sumAN = 0.;
                double sumErr = 0.;
                double sumSysErr[5];
                memset(sumSysErr, 0., sizeof(sumSysErr));
            
                for(int xf=(i*2); xf<(i*2+2); xf++){
                    int gNum = mGraphAN_pTSummary[dle][xf]->GetN();
                    for(int p=0; p<gNum; p++){
                        double pt, an;
                        mGraphAN_pTSummary[dle][xf] -> GetPoint(p, pt, an);
                        double err = mGraphAN_pTSummary[dle][xf] -> GetErrorY(p);
                        if(ptBoundL <= pt && pt < ptBoundU){
                            double w = 1./(err*err);

                            sumPt += (pt*w);
                            sumAN += (an*w);
                            sumErr += w;
                            for(int s=0; s<5; s++){
                                if(s==1){continue;}
                                double sysErr = mGraphAN_pTSummarySysErr[dle][s][xf] -> GetErrorY(p);
                                sumSysErr[s] += (sysErr*sysErr*w*w);
                            }
                        }
                    }
                }

                if(sumErr > 0.1){
                    sumPt /= sumErr;
                    sumAN /= sumErr;
                    double statErr = 1./sqrt(sumErr);

                    int gNum = mGraphAN_pTSummary_merge[dle][i]->GetN();
                    mGraphAN_pTSummary_merge[dle][i] -> SetPoint(gNum, sumPt, sumAN);
                    mGraphAN_pTSummary_merge[dle][i] -> SetPointError(gNum, 0., statErr);

                    for(int s=0; s<5; s++){
                        if(s==1){continue;}
                        double sysErr = sqrt(sumSysErr[s])/sumErr;
                        sumSysErr[s] = sysErr;
                        mGraphAN_pTSummarySysErr_merge[dle][s][i] -> SetPoint(gNum, sumPt, sumAN);
                        mGraphAN_pTSummarySysErr_merge[dle][s][i] -> SetPointError(gNum, 0.01, sysErr);
                    }
                    if(i == 1){
                        cout << mOptContainer->GetDLEName(dle) << " PT " << endl;
                        cout << fixed;     
                        cout.precision(5);  
                        cout << "& " << sumPt << " & " << sumSysErr[3] << " & " << sumSysErr[1] << " & " << sumSysErr[0] << " & " << sumSysErr[2] << " & " << sumSysErr[4] << endl;
                    }
                }
            }
        }
    }


    for(int dle=0; dle<kDLENum; dle++){
        for(int i=0; i<2; i++){
            for(int xf=0; xf<xfNum; xf++){
                double xfBoundL = mBinning[0][0][0] -> GetGlobalXfBinBoundary(xf);
                double xfBoundU = mBinning[0][0][0] -> GetGlobalXfBinBoundary(xf+1);

                double sumXf = 0.;
                double sumAN = 0.;
                double sumErr = 0.;
                double sumSysErr[5];
                memset(sumSysErr, 0., sizeof(sumSysErr));
            
                for(int pt=(i*2); pt<(i*2+2); pt++){
                    int gNum = mGraphAN_xFSummary[dle][pt]->GetN();
                    for(int p=0; p<gNum; p++){
                        double xf, an;
                        mGraphAN_xFSummary[dle][pt] -> GetPoint(p, xf, an);
                        double err = mGraphAN_xFSummary[dle][pt] -> GetErrorY(p);
                        if(xfBoundL <= xf && xf < xfBoundU){
                            double w = 1./(err*err);

                            sumXf += (xf*w);
                            sumAN += (an*w);
                            sumErr += w;
                            for(int s=0; s<5; s++){
                                if(s==1){continue;}
                                double sysErr = mGraphAN_xFSummarySysErr[dle][s][pt] -> GetErrorY(p);
                                sumSysErr[s] += (sysErr*sysErr*w*w);
                            }
                        }
                    }
                }

                if(sumErr > 0.1){
                    sumXf /= sumErr;
                    sumAN /= sumErr;
                    double statErr = 1./sqrt(sumErr);

                    int gNum = mGraphAN_xFSummary_merge[dle][i]->GetN();
                    mGraphAN_xFSummary_merge[dle][i] -> SetPoint(gNum, sumXf, sumAN);
                    mGraphAN_xFSummary_merge[dle][i] -> SetPointError(gNum, 0., statErr);

                    for(int s=0; s<5; s++){
                        if(s==1){continue;}
                        double sysErr = sqrt(sumSysErr[s])/sumErr;
                        sumSysErr[s] = sysErr;
                        mGraphAN_xFSummarySysErr_merge[dle][s][i] -> SetPoint(gNum, sumXf, sumAN);
                        mGraphAN_xFSummarySysErr_merge[dle][s][i] -> SetPointError(gNum, 0.01, sysErr);
                    }

                    if(i == 1){
                        cout << mOptContainer->GetDLEName(dle) << " XF " << endl;
                        cout << fixed;     
                        cout.precision(5);  
                        cout << "& " << sumXf << " & " << sumSysErr[3] << " & " << sumSysErr[1] << " & " << sumSysErr[0] << " & " << sumSysErr[2] << " & " << sumSysErr[4] << endl;
                    }
                }
            }
        }
    }



}

void RHICfAsymmetry::SystematicErrorCalculatePi0()
{
    int ptNum = 0;
    int xfNum = 0;
    for(int run=0; run<kRunNum; run++){
        TString runName = GetOptContainer()->GetRunTypeName(run);
        GetOptContainer()->SetRunType(runName);
        if(GetOptContainer()->GetRunType() != run && GetOptContainer()->GetRunType() != kALLRun){continue;}
        int beamRefNum = (run == kTOPRun)? kBeamTOPRefNum : 1;

        if(mBinning[run][0][0] -> GetBinningTableFlag() == kNotExist){continue;}
        ptNum = mBinning[run][0][0] -> GetGlobalPtBinNum();
        xfNum = mBinning[run][0][0] -> GetGlobalXfBinNum();

        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        if(mBinning[run][0][0] -> GetBinningTableFlag() == kNotExist){continue;}
                        double ptMean = mBinning[run][0][0] -> GetPtBinMean(run, type, dle, pt, xf);
                        double xfMean = mBinning[run][0][0] -> GetXfBinMean(run, type, dle, pt, xf);

                        double sumAN = 0.;
                        double sumErr = 0.;

                        double sumBkgAN = 0.;
                        double sumBkgErr = 0.;

                        double sumSubtAN = 0.;
                        double sumSubtErr = 0.;

                        int graphIdx = -1;
                        for(int i=0; i<kBeamMetNum; i++){
                            for(int j=0; j<beamRefNum; j++){
                                graphIdx++;
                                if(mMassFitting[run][i][j] -> GetMassTableFlag() == kNotExist){continue;}
                                if(mDilution[run][i][j] -> GetDilutionTableFlag() == kNotExist){continue;}
                                if(mPolarization[run][i][j] -> GetPolarizationTableFlag() == kNotExist){continue;}

                                for(int s=0; s<3; s++){
                                    if(s==1){continue;}
                                    for(int k=0; k<2; k++){     
                                        double sign = (k==0)? +1. : -1.;
                                                
                                        double dilutionFactor = mDilution[run][i][j] -> GetDilutionFactor(run, type, dle, pt, xf);
                                        double dilutionError = mDilution[run][i][j] -> GetDilutionFactorErr(run, type, dle, pt, xf);
                                        if(dilutionFactor <= 0.1){continue;}
                                        if(s==0){dilutionFactor = dilutionFactor + sign*dilutionError;}

                                        double allCounts = mMassFitting[run][i][j] -> GetMassAllCounts(false, run, type, dle, pt, xf);
                                        double signalCounts = mMassFitting[run][i][j] -> GetMassSignalCounts(false, run, type, dle, pt, xf);
                                        double bkgCounts = mMassFitting[run][i][j] -> GetMassBkgCounts(false, run, type, dle, pt, xf);

                                        double anSubtWeight = 0.;
                                        double errSubtWeight = 0.;

                                        for(int fill=0; fill<kFillNum; fill++){
                                            int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                                            if(run != runIdx){continue;}
                                            double pol =  BeamPolarization(fill);
                                            if(s==1){pol += (sign*BeamPolarizationError(fill));}
                                        
                                            double luminosity = RelativeLuminosity(fill);

                                            double up = mPolarization[run][i][j] -> GetPolNum(fill, type, dle, pt, xf, true);
                                            double down = mPolarization[run][i][j] -> GetPolNum(fill, type, dle, pt, xf, false);
                                            if(up <= 10 || down <= 10){continue;}

                                            double an = -1./(pol*dilutionFactor) * (up - luminosity*down)/(up + luminosity*down);
                                            if(run == kTLRun){an = -an;}
                                            double err = 2.*luminosity*sqrt(up*down*(up+down))/(pol*dilutionFactor*pow((up + luminosity*down), 2.));
                                            
                                            double upBkg = mPolarization[run][i][j] -> GetBkgPolNum(fill, type, dle, pt, xf, true);
                                            double downBkg = mPolarization[run][i][j] -> GetBkgPolNum(fill, type, dle, pt, xf, false);

                                            double anBkg = -1./(pol*dilutionFactor) * (upBkg - luminosity*downBkg)/(upBkg + luminosity*downBkg);
                                            if(run == kTLRun){anBkg = -anBkg;}
                                            double errBkg = 2.*luminosity*sqrt(upBkg*downBkg*(upBkg+downBkg))/(pol*dilutionFactor*pow((upBkg + luminosity*downBkg), 2.));

                                            double massMean = mMassFitting[run][i][j] -> GetMassMean(run, type, dle, pt, xf);
                                            if(massMean <= 10.){continue;}
                                            if(upBkg <= 1 || downBkg <= 1){continue;}

                                            double BSratio = mMassFitting[run][i][j] -> GetBSRatio_GPR(run, type, dle, pt, xf);

                                            if(s==2){
                                                if(k==0){BSratio = mMassFitting[run][i][j] -> GetBSRatioUpper_GPR(run, type, dle, pt, xf);}
                                                if(k==1){BSratio = mMassFitting[run][i][j] -> GetBSRatioLower_GPR(run, type, dle, pt, xf);}
                                            }

                                            double ASratio = 1. + BSratio;
                                            double subtAN = ASratio*an - BSratio*anBkg;
                                            double subtErr = sqrt(ASratio*ASratio*err*err + BSratio*BSratio*errBkg*errBkg);

                                            if(err < 0.2){
                                                double wSubt = 1./(subtErr*subtErr);
                                                anSubtWeight += (subtAN*wSubt);
                                                errSubtWeight += wSubt;
                                            }
                                        }

                                        if(errSubtWeight > 0.1){
                                            anSubtWeight /= errSubtWeight;
                                            errSubtWeight = 1./sqrt(errSubtWeight);
                                            if(errSubtWeight > 0.15){continue;}

                                            int pTSubtGraphNum = mGraphANSysErr_Calcul_subt_pT[run][type][dle][graphIdx][s][k][xf]->GetN();
                                            mGraphANSysErr_Calcul_subt_pT[run][type][dle][graphIdx][s][k][xf] -> SetPoint(pTSubtGraphNum, ptMean, anSubtWeight);
                                            mGraphANSysErr_Calcul_subt_pT[run][type][dle][graphIdx][s][k][xf] -> SetPointError(pTSubtGraphNum, 0., errSubtWeight);

                                            int xFSubtGraphNum = mGraphANSysErr_Calcul_subt_xF[run][type][dle][graphIdx][s][k][pt]->GetN();
                                            mGraphANSysErr_Calcul_subt_xF[run][type][dle][graphIdx][s][k][pt] -> SetPoint(xFSubtGraphNum, xfMean, anSubtWeight);
                                            mGraphANSysErr_Calcul_subt_xF[run][type][dle][graphIdx][s][k][pt] -> SetPointError(xFSubtGraphNum, 0., errSubtWeight);
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
}

void RHICfAsymmetry::DrawAsymmetryGraph()
{
    cout << "RHICfAsymmetry::DrawAsymmetryGraph() " << endl;
    TString figurePath = mOptContainer -> GetFigurePath() + "/Asymmetry";
    TString ParticleTypeName = mOptContainer -> GetParticleRunName();
    int particleTypeIdx = mOptContainer -> GetParticleRunIdx();

    int anTypeNum = 0;
    const int runColor[kRunNum] = {3, 1, 2};
    const int typeStyle[kTypeNum] = {20, 24};
    const int binColor[7] = {1, 2, 4, 8, 51, 90, 93}; 
    const int markerStyleDLE[kDLENum]= {22, 23, 21, 20};
    const int markerColorDLE[kDLENum]= {8, 4, 2, 1};

    double legendPosY[2];
    double drawANRange[2]; 
    double latexPosY[2];

    legendPosY[0] = 0.23;
    legendPosY[1] = 0.485;
    drawANRange[0] = -0.3;
    drawANRange[1] = 0.3;
    latexPosY[0] = 0.13;
    latexPosY[1] = 0.06;
    anTypeNum = 3;


    TCanvas* cAsymmetry = new TCanvas("cAsymmetry", "", 600.*4., 600.*3.);

    TLatex* latex = new TLatex();
    TLegend* legend = new TLegend(0.65, legendPosY[0], 0.89, legendPosY[1]);
    legend -> SetBorderSize(0);

    TLine* line = new TLine(0., 0., 1., 0.);
    line -> SetLineStyle(2);
    line -> SetLineColor(kBlack);

    int cIdxArrXf[6] = {1, 2, 5, 6, 9, 10};
    int cIdxArrPt[4] = {3, 4, 7, 8};

    int availableRunIdx = -1;
    for(int run=0; run<kRunNum; run++){
        if(mBinning[run][0][0] -> GetBinningTableFlag() == kNotExist){continue;}
        availableRunIdx = run;
    }
    if(availableRunIdx == -1){return;}

    const int ptNum = mBinning[availableRunIdx][0][0] -> GetGlobalPtBinNum();
    const int xfNum = mBinning[availableRunIdx][0][0] -> GetGlobalXfBinNum();

    TGraph* base[2];
    for(int ptxf=0; ptxf<2; ptxf++){
        base[ptxf] = new TGraph();
        base[ptxf] -> SetPoint(0, 0., -1.);
        base[ptxf] -> SetPoint(1, 1.5, 1.);
        base[ptxf] -> SetMarkerSize(0);
        base[ptxf] -> SetMarkerStyle(20);
        base[ptxf] -> SetLineWidth(0);
        base[ptxf] -> GetYaxis() -> SetRangeUser(drawANRange[0]-0.05, drawANRange[1]); 
        base[ptxf] -> GetXaxis() -> SetRangeUser(0., 1.);

        TString xTitle = (ptxf==0)? "p_{T} [GeV/c]" : "x_{F}";
        base[ptxf] -> SetTitle(Form("; %s; A_{N}", xTitle.Data()));
    }

    for(int dle=0; dle<kDLENum; dle++){
        TString dleName = mOptContainer->GetDLEName(dle);
        if(dle == kDLENum-1){dleName = "Inclusive";}
        // if(dle != kALLDLE){continue;}

        for(int step=0; step<3; step++){
            TString stepName = "";
            if(step == 0){stepName = "Inc";}
            if(step == 1){stepName = "Bkg";}
            if(step == 2){stepName = "Subt";}

            cAsymmetry -> Clear();
            cAsymmetry -> Divide(4, 3);
        
            for(int ptxf=0; ptxf<2; ptxf++){
                bool isPt = (ptxf==0)? true : false;
                
                int binNum = (isPt)? xfNum : ptNum;

                for(int bin=0; bin<binNum; bin++){
                    const int cIdx = (isPt)? cIdxArrPt[bin] : cIdxArrXf[bin];
                    cAsymmetry -> cd(cIdx);
                    gPad -> SetGrid(1,1);
                    base[ptxf] -> Draw("ap");
                    legend -> Clear("ICESM");
                    for(int run=0; run<kRunNum; run++){
                        TString runName = mOptContainer -> GetRunTypeName(run);
                        for(int type=0; type<kTypeNum; type++){                        
                            TGraphErrors* graph = 0;
                            if(step == 0){graph = GetANGraph(isPt, 4, run, type, dle, bin);}
                            if(step == 1){graph = GetANBkgGraph(isPt, 4, run, type, dle, bin);}
                            if(step == 2){graph = GetANSubtGraph(isPt, 4, run, type, dle, bin);}
                            graph -> SetMarkerColor(runColor[run]);
                            graph -> SetMarkerStyle(typeStyle[type]);
                            graph -> SetMarkerSize(1.2);
                            graph -> SetLineWidth(1.5);
                            graph -> SetLineColor(runColor[run]);
                        
                            if(step == 2){
                                if(isPt){
                                    mGraphANSysErr_Subt_pT[run][type][dle][4][bin] -> SetLineColor(runColor[run]);
                                    mGraphANSysErr_Subt_pT[run][type][dle][4][bin] -> SetLineColor(runColor[run]);
                                    mGraphANSysErr_Subt_pT[run][type][dle][4][bin] -> SetFillStyle(0);
                                    mGraphANSysErr_Subt_pT[run][type][dle][4][bin] -> Draw("5, same");
                                }
                                else{
                                    mGraphANSysErr_Subt_xF[run][type][dle][4][bin] -> SetLineColor(runColor[run]);
                                    mGraphANSysErr_Subt_xF[run][type][dle][4][bin] -> SetLineColor(runColor[run]);
                                    mGraphANSysErr_Subt_xF[run][type][dle][4][bin] -> SetFillStyle(0);
                                    mGraphANSysErr_Subt_xF[run][type][dle][4][bin] -> Draw("5, same");    
                                }
                            }
                            graph -> Draw("p, 0, same");
                            legend -> AddEntry(graph, Form("%s Type%i", runName.Data(), type+1), "p");
                        }                    
                    }
                    
                    legend -> Draw("same");
                    line -> Draw("same");
                    latex -> DrawLatexNDC(0.13, latexPosY[0], GetSystemString());
                    latex -> DrawLatexNDC(0.13, latexPosY[0]+latexPosY[1], Form("%s", dleName.Data()));
                    latex -> DrawLatexNDC(0.33, latexPosY[0]+latexPosY[1], "#eta > 6");
                    TString binName = "";
                    if(isPt){
                        binName = "x_{F}";
                        latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(bin), binName.Data(), mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(bin+1) ));
                    }
                    else{
                        binName = "p_{T}";
                        latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(bin), binName.Data(), mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(bin+1) ));
                    }
                    
                }
            }
            cAsymmetry -> Draw();
            cAsymmetry -> SaveAs(Form("%s/AN_Summed_%s_%s.pdf", figurePath.Data(), stepName.Data(), dleName.Data()));
        }
    }

    // ======================== Summary plot ==========================

    TCanvas* cAsymSummary = new TCanvas("cAsymSummary", "", 1200., 600.);
    TLegend* sumLegend[2];
    TLegend* sumLegendDLE[2];
    for(int i=0; i<2; i++){
        sumLegend[i] = new TLegend(0.13, 0.65, 0.52, 0.89);
        sumLegendDLE[i] = new TLegend(0.50, 0.65, 0.89, 0.89);
        sumLegend[i] -> SetBorderSize(0);
        sumLegendDLE[i] -> SetBorderSize(0);
    } 

    int summaryColor[6] = {8, 66, 2, 9, 6, 1};

    for(int dle=0; dle<kDLENum; dle++){
        TString dleName = mOptContainer->GetDLEName(dle);
        if(dle == kDLENum-1){dleName = "Inclusive";}

        cAsymSummary -> Clear();
        cAsymSummary -> Divide(2,1);

        cAsymSummary -> cd(1);
        sumLegend[0] -> Clear();
        sumLegendDLE[0] -> Clear();
        base[0] -> GetYaxis()->SetTitleOffset(1.47);
        base[0] -> GetYaxis() -> SetRangeUser(-0.13, drawANRange[1]+0.02); 
        if(dle == kDDLE){base[0] -> GetYaxis() -> SetRangeUser(-0.155, drawANRange[1]+0.02); }
        base[0] -> GetXaxis() -> SetRangeUser(0., 1.);
        base[0] -> Draw("ap");
        for(int xf=0; xf<xfNum; xf++){
            mGraphAN_pTSummary[dle][xf] -> SetMarkerColor(summaryColor[xf]);
            mGraphAN_pTSummary[dle][xf] -> SetMarkerStyle(20);
            mGraphAN_pTSummary[dle][xf] -> SetMarkerSize(1.2);
            mGraphAN_pTSummary[dle][xf] -> SetLineWidth(1.2);
            mGraphAN_pTSummary[dle][xf] -> SetLineColor(summaryColor[xf]);

                mGraphAN_pTSummarySysErr[dle][4][xf] -> SetLineColor(summaryColor[xf]);
                mGraphAN_pTSummarySysErr[dle][4][xf]  -> SetFillStyle(0);
                mGraphAN_pTSummarySysErr[dle][4][xf] -> Draw("5, same");

            mGraphAN_pTSummary[dle][xf] -> Draw("p, same");
            

            sumLegend[0] -> AddEntry(mGraphAN_pTSummary[dle][xf], Form("%.2f < x_{F} < %.2f", mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(xf), mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(xf+1)), "p");
            sumLegendDLE[0] -> AddEntry(mGraphAN_pTSummary[dle][xf], Form("%.2f < x_{F} < %.2f", mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(xf), mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(xf+1)), "p");
        }

        line -> Draw("same");
        if(dle == kALLDLE || dle == kNDLE){sumLegend[0] -> Draw("same");}
        if(dle == kSDLE || dle == kDDLE){sumLegendDLE[0] -> Draw("same");}
        latex -> DrawLatexNDC(0.16, latexPosY[0], GetSystemString());
        latex -> DrawLatexNDC(0.16, latexPosY[0]+latexPosY[1], Form("%s", dleName.Data()));
        latex -> DrawLatexNDC(0.36, latexPosY[0]+latexPosY[1], "#eta > 6");

        cAsymSummary -> cd(2);
        sumLegend[1] -> Clear();
        sumLegendDLE[1] -> Clear();
        base[1] -> GetYaxis()->SetTitleOffset(1.47);
        base[1] -> GetYaxis() -> SetRangeUser(-0.1, drawANRange[1]+0.02); 
        if(dle == kSDLE || dle == kDDLE){base[1] -> GetYaxis() -> SetRangeUser(-0.155, drawANRange[1]+0.02); }
        base[1] -> GetXaxis() -> SetRangeUser(0., 1.);
        base[1] -> Draw("ap");
        for(int pt=0; pt<ptNum; pt++){
            mGraphAN_xFSummary[dle][pt] -> SetMarkerColor(summaryColor[pt]);
            mGraphAN_xFSummary[dle][pt] -> SetMarkerStyle(20);
            mGraphAN_xFSummary[dle][pt] -> SetMarkerSize(1.2);
            mGraphAN_xFSummary[dle][pt] -> SetLineWidth(1.5);
            mGraphAN_xFSummary[dle][pt] -> SetLineColor(summaryColor[pt]);


            mGraphAN_xFSummarySysErr[dle][4][pt] -> SetLineColor(summaryColor[pt]);
            mGraphAN_xFSummarySysErr[dle][4][pt]  -> SetFillStyle(0);
            mGraphAN_xFSummarySysErr[dle][4][pt] -> Draw("5, same");

            mGraphAN_xFSummary[dle][pt] -> Draw("p, same");

            sumLegend[1] -> AddEntry(mGraphAN_xFSummary[dle][pt], Form("%.2f < p_{T} < %.2f GeV/c", mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(pt), mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(pt+1)), "p");
            sumLegendDLE[1] -> AddEntry(mGraphAN_xFSummary[dle][pt], Form("%.2f < p_{T} < %.2f GeV/c", mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(pt), mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(pt+1)), "p");
        }

        line -> Draw("same");
        if(dle == kALLDLE || dle == kNDLE || dle == kDDLE){sumLegend[1] -> Draw("same");}
        if(dle == kSDLE){sumLegendDLE[1] -> Draw("same");}
        latex -> DrawLatexNDC(0.16, latexPosY[0], GetSystemString());
        latex -> DrawLatexNDC(0.16, latexPosY[0]+latexPosY[1], Form("%s", dleName.Data()));
        latex -> DrawLatexNDC(0.36, latexPosY[0]+latexPosY[1], "#eta > 6");

        cAsymSummary -> Draw();
        cAsymSummary -> SaveAs(Form("%s/AN_Summary_%s.pdf", figurePath.Data(), dleName.Data()));
    }


    TCanvas* cAsymDLE = new TCanvas("cAsymDLE", "", 1200., 1200.);
    cAsymDLE -> Divide(2,2);

    int dleColor[4] = {8, 66, 2, 9};
    int dleStyle[4] = {22, 23, 33, 20};

    TLegend* legDLE = new TLegend(0.13, 0.64, 0.44, 0.87);
    legDLE -> SetBorderSize(0);

    TLine* lineDLE0 = new TLine(0., 0., 0.44, 0.);
    // TLine* lineDLE0 = new TLine(0., 0., 1., 0.);
    lineDLE0 -> SetLineColor(kBlack);
    lineDLE0 -> SetLineStyle(2);

    TLine* lineDLE1 = new TLine(0.2, 0., 0.75, 0.);
    lineDLE1 -> SetLineColor(kBlack);
    lineDLE1 -> SetLineStyle(2);

    TGraph* clone = (TGraph*)base[0] ->Clone("testtest");

    for(int i=0; i<2; i++){
        const int cIdx = i*2+1;
        cAsymDLE -> cd(cIdx);
        if(i==0){clone -> GetYaxis() -> SetRangeUser(-0.135, 0.29);}
        else{clone -> GetYaxis() -> SetRangeUser(-0.135, 0.21);}
        clone -> GetXaxis() -> SetRangeUser(0., 0.44);
        clone -> Draw("ap");

        gPad->SetTicks(1, 1);

        for(int dle=0; dle<kDLENum; dle++){
            TString dleName = mOptContainer->GetDLEName(dle);
            if(dle == kDLENum-1){dleName = "Inclusive";}

            mGraphAN_pTSummary_merge[dle][i] -> SetMarkerColor(dleColor[dle]);
            mGraphAN_pTSummary_merge[dle][i] -> SetMarkerStyle(dleStyle[dle]);
            mGraphAN_pTSummary_merge[dle][i] -> SetMarkerSize(1.5);
            if(dle == kNDLE){mGraphAN_pTSummary_merge[dle][i] -> SetMarkerSize(2.2);}
            mGraphAN_pTSummary_merge[dle][i] -> SetLineWidth(1.5);
            mGraphAN_pTSummary_merge[dle][i] -> SetLineColor(dleColor[dle]);
            mGraphAN_pTSummary_merge[dle][i] -> Draw("p, same");

            mGraphAN_pTSummarySysErr_merge[dle][4][i] -> SetFillStyle(0);
            mGraphAN_pTSummarySysErr_merge[dle][4][i] -> SetLineColor(dleColor[dle]);
            mGraphAN_pTSummarySysErr_merge[dle][4][i] -> Draw("5, same");

            if(i == 0){
                legDLE -> AddEntry(mGraphAN_pTSummary_merge[dle][0], Form("%s", dleName.Data()), "p");
            }
        }
        legDLE -> Draw("same");
        latex -> SetTextColor(kRed);
        latex -> SetTextSize(0.04);
        latex -> DrawLatexNDC(0.36, 0.81, "RHICf and STAR Preliminary");

        latex -> SetTextSize(0.038);
        latex -> SetTextColor(kBlack);
        latex -> DrawLatexNDC(0.36, 0.76, "#font[12]{#bf{3.8% beam pol. uncer. not shown}}");

        latex -> SetTextSize(0.045);
        latex -> SetTextColor(kBlack);
        latex -> DrawLatexNDC(0.16, latexPosY[0]+0.005, "p^{#uparrow}+p #rightarrow #pi^{0} + X @ #sqrt{s} = 510 GeV");
        latex -> DrawLatexNDC(0.16, latexPosY[0]+latexPosY[1]+0.01, Form("6 < #eta ,   %.2f < x_{F} < %.2f",  mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(i*2), mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary((i*2)+2)));
        // latex -> DrawLatexNDC(0.16, latexPosY[0]+latexPosY[1], Form("6 < #eta , %.2f < x_{F} < %.2f",  mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(i*2), mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(xfNum)));
        lineDLE0 -> Draw("same");
    }

    for(int i=0; i<2; i++){
        const int cIdx = i*2+2;
        cAsymDLE -> cd(cIdx);
        base[1] ->  GetYaxis() -> SetRangeUser(-0.085, 0.21); 
        base[1] -> GetXaxis() -> SetRangeUser(0.2, 0.75);
        base[1] -> Draw("ap");

        gPad->SetTicks(1, 1);

        for(int dle=0; dle<kDLENum; dle++){
            TString dleName = mOptContainer->GetDLEName(dle);
            if(dle == kDLENum-1){dleName = "Inclusive";}

            mGraphAN_xFSummary_merge[dle][i] -> SetMarkerColor(dleColor[dle]);
            mGraphAN_xFSummary_merge[dle][i] -> SetMarkerStyle(dleStyle[dle]);
            mGraphAN_xFSummary_merge[dle][i] -> SetMarkerSize(1.5);
            if(dle == kNDLE){mGraphAN_xFSummary_merge[dle][i] -> SetMarkerSize(2.2);}
            mGraphAN_xFSummary_merge[dle][i] -> SetLineWidth(1.5);
            mGraphAN_xFSummary_merge[dle][i] -> SetLineColor(dleColor[dle]);
            mGraphAN_xFSummary_merge[dle][i] -> Draw("p, same");

            mGraphAN_xFSummarySysErr_merge[dle][4][i] -> SetFillStyle(0);
            mGraphAN_xFSummarySysErr_merge[dle][4][i] -> SetLineColor(dleColor[dle]);
            mGraphAN_xFSummarySysErr_merge[dle][4][i] -> Draw("5, same");
        }
        legDLE -> Draw("same");

        latex -> SetTextColor(kRed);
        latex -> SetTextSize(0.04);
        latex -> DrawLatexNDC(0.36, 0.81, "RHICf and STAR Preliminary");

        latex -> SetTextSize(0.038);
        latex -> SetTextColor(kBlack);
        latex -> DrawLatexNDC(0.36, 0.76, "#font[12]{#bf{3.8% beam pol. uncer. not shown}}");

        latex -> SetTextSize(0.045);
        latex -> SetTextColor(kBlack);
        latex -> DrawLatexNDC(0.16, latexPosY[0]+0.005, "p^{#uparrow}+p #rightarrow #pi^{0} + X @ #sqrt{s} = 510 GeV");
        latex -> DrawLatexNDC(0.16, latexPosY[0]+latexPosY[1]+0.01, Form("6 < #eta ,   %.2f < p_{T} < %.2f GeV/c",  mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(i*2), mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary((i*2)+2)));
        lineDLE1 -> Draw("same");
    }
    cAsymDLE -> Draw();
    cAsymDLE -> SaveAs(Form("%s/AN_DLE.pdf",figurePath.Data()));




    int metColor[2] = {2, 4};
    for(int step=2; step<3; step++){
        TString stepName = "";
        if(step == 0){stepName = "Inc";}
        if(step == 1){stepName = "Bkg";}
        if(step == 2){stepName = "Subt";}

        cAsymmetry -> Clear();
        cAsymmetry -> Divide(4, 3);
    
        for(int ptxf=0; ptxf<2; ptxf++){
            bool isPt = (ptxf==0)? true : false;
            
            int binNum = (isPt)? xfNum : ptNum;

            for(int bin=0; bin<binNum; bin++){
                const int cIdx = (isPt)? cIdxArrPt[bin] : cIdxArrXf[bin];
                cAsymmetry -> cd(cIdx);
                gPad -> SetGrid(1,1);
                base[ptxf] -> GetYaxis()->SetRangeUser(-0.2, 0.31);
                base[ptxf] -> Draw("ap");
                legend -> Clear("ICESM");
                for(int run=0; run<kRunNum; run++){
                    // if(run != kTSRun){continue;}
                    TString runName = mOptContainer -> GetRunTypeName(run);
                    for(int type=0; type<kTypeNum; type++){    
                        for(int met=0; met<2; met++){     
                            TString metName = (met==0)? "BeamHit" : "BeamScan";               
                            TGraphErrors* graph = 0;
                            if(step == 0){graph = GetANGraph(isPt, met, run, type, 3, bin);}
                            if(step == 1){graph = GetANBkgGraph(isPt, met, run, type, 3, bin);}
                            if(step == 2){graph = GetANSubtGraph(isPt, met, run, type, 3, bin);}

                            for(int k=0; k<2; k++){
                                if(isPt){
                                    // mGraphANSysErr_Subt_pT[run][type][dle][5][bin] -> SetLineColor(runColor[run]);
                                    // mGraphANSysErr_Subt_pT[run][type][dle][5][bin] -> SetLineColor(runColor[run]);
                                    // mGraphANSysErr_Subt_pT[run][type][dle][5][bin] -> SetFillStyle(0);
                                    // mGraphANSysErr_Subt_pT[run][type][dle][5][bin] -> Draw("p, same");

                                    mGraphANSysErr_Calcul_subt_pT[run][type][3][met][2][k][bin] -> SetMarkerColor(metColor[met]);
                                    mGraphANSysErr_Calcul_subt_pT[run][type][3][met][2][k][bin] -> SetMarkerStyle(typeStyle[type]);
                                    mGraphANSysErr_Calcul_subt_pT[run][type][3][met][2][k][bin] -> SetLineColor(metColor[met]);
                                    mGraphANSysErr_Calcul_subt_pT[run][type][3][met][2][k][bin] -> SetLineColor(metColor[met]);
                                    mGraphANSysErr_Calcul_subt_pT[run][type][3][met][2][k][bin] -> SetFillStyle(0);
                                    mGraphANSysErr_Calcul_subt_pT[run][type][3][met][2][k][bin] -> Draw("p, same");
                                }
                                else{
                                    // mGraphANSysErr_Subt_xF[run][type][dle][5][bin] -> SetLineColor(runColor[run]);
                                    // mGraphANSysErr_Subt_xF[run][type][dle][5][bin] -> SetLineColor(runColor[run]);
                                    // mGraphANSysErr_Subt_xF[run][type][dle][5][bin] -> SetFillStyle(0);
                                    // mGraphANSysErr_Subt_xF[run][type][dle][5][bin] -> Draw("p, same");  

                                    mGraphANSysErr_Calcul_subt_xF[run][type][3][met][2][k][bin] -> SetMarkerColor(metColor[met]);
                                    mGraphANSysErr_Calcul_subt_xF[run][type][3][met][2][k][bin] -> SetMarkerStyle(typeStyle[type]);
                                    mGraphANSysErr_Calcul_subt_xF[run][type][3][met][2][k][bin] -> SetLineColor(metColor[met]);
                                    mGraphANSysErr_Calcul_subt_xF[run][type][3][met][2][k][bin] -> SetLineColor(metColor[met]);
                                    mGraphANSysErr_Calcul_subt_xF[run][type][3][met][2][k][bin] -> SetFillStyle(0);
                                    mGraphANSysErr_Calcul_subt_xF[run][type][3][met][2][k][bin] -> Draw("p, same");  
                                }

                            }
                            graph -> SetMarkerColor(metColor[met]);
                            graph -> SetMarkerStyle(typeStyle[type]);
                            graph -> SetMarkerSize(1.2);
                            graph -> SetLineWidth(1.5);
                            graph -> SetLineColor(metColor[met]);
                        
                            // graph -> Draw("p, same");
                            legend -> AddEntry(graph, Form("%s Type%i", metName.Data(), type+1), "p");
                        }
                    }                    
                }
                
                legend -> Draw("same");
                line -> Draw("same");
                latex -> DrawLatexNDC(0.13, latexPosY[0], GetSystemString());
                latex -> DrawLatexNDC(0.13, latexPosY[0]+latexPosY[1], Form("%s", "Inclusive"));
                latex -> DrawLatexNDC(0.33, latexPosY[0]+latexPosY[1], "#eta > 6");
                TString binName = "";
                if(isPt){
                    binName = "x_{F}";
                    latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(bin), binName.Data(), mBinning[availableRunIdx][0][0]->GetGlobalXfBinBoundary(bin+1) ));
                }
                else{
                    binName = "p_{T}";
                    latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(bin), binName.Data(), mBinning[availableRunIdx][0][0]->GetGlobalPtBinBoundary(bin+1) ));
                }
                
            }
        }
        cAsymmetry -> Draw();
        cAsymmetry -> SaveAs(Form("%s/AN_BeamMet_BS_%s.pdf", figurePath.Data(), stepName.Data()));
    }
}

TString RHICfAsymmetry::GetSystemString()
{
    TString text = "";
    if(mOptContainer->GetParticleRunIdx() == kGammaRun){
        text = "p^{#uparrow}+p #rightarrow #gamma + X @ #sqrt{s} = 510 GeV";
    }
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){
        text = "p^{#uparrow}+p #rightarrow #pi^{0} + X @ #sqrt{s} = 510 GeV";
    }
    if(mOptContainer->GetParticleRunIdx() == kNeutronRun){
        text = "p^{#uparrow}+p #rightarrow n + X @ #sqrt{s} = 510 GeV";
    }
    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){
        text = "p^{#uparrow}+p #rightarrow #Lambda^{0} + X @ #sqrt{s} = 510 GeV";
    }
    return text;
}

void RHICfAsymmetry::InitGraph()
{
    int ptNum = mBinning[0][0][0] -> GetGlobalPtBinNum();
    int xfNum = mBinning[0][0][0] -> GetGlobalXfBinNum();

    for(int dle=0; dle<kDLENum; dle++){
        TString dleName = mOptContainer -> GetDLEName(dle);
        for(int i=0; i<5; i++){
            for(int run=0; run<kRunNum; run++){
                TString runName = mOptContainer -> GetRunTypeName(run);

                for(int type=0; type<kTypeNum; type++){
                    mGraphAN_pT[run][type][dle][i].resize(xfNum);
                    mGraphAN_Bkg_pT[run][type][dle][i].resize(xfNum);
                    mGraphAN_Subt_pT[run][type][dle][i].resize(xfNum);
                    for(int xf=0; xf<xfNum; xf++){
                        mGraphAN_pT[run][type][dle][i][xf] = new TGraphErrors();
                        mGraphAN_pT[run][type][dle][i][xf] -> SetName(Form("AN_pT_%s_type%i_%s_xfBin%i", runName.Data(), type, dleName.Data(), xf));
                        mGraphAN_pT[run][type][dle][i][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");

                        mGraphAN_Bkg_pT[run][type][dle][i][xf] = new TGraphErrors();
                        mGraphAN_Bkg_pT[run][type][dle][i][xf] -> SetName(Form("AN_Bkg_pT_%s_type%i_%s_xfBin%i", runName.Data(), type, dleName.Data(), xf));
                        mGraphAN_Bkg_pT[run][type][dle][i][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");


                        mGraphAN_Subt_pT[run][type][dle][i][xf] = new TGraphErrors();
                        mGraphAN_Subt_pT[run][type][dle][i][xf] -> SetName(Form("AN_Subt_pT_%s_type%i_%s_xfBin%i", runName.Data(), type, dleName.Data(), xf));
                        mGraphAN_Subt_pT[run][type][dle][i][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
                    }

                    mGraphAN_xF[run][type][dle][i].resize(ptNum);
                    mGraphAN_Bkg_xF[run][type][dle][i].resize(ptNum);
                    mGraphAN_Subt_xF[run][type][dle][i].resize(ptNum);
                    for(int pt=0; pt<ptNum; pt++){
                        mGraphAN_xF[run][type][dle][i][pt] = new TGraphErrors();
                        mGraphAN_xF[run][type][dle][i][pt] -> SetName(Form("AN_xF_%s_type%i_%s_ptBin%i", runName.Data(), type, dleName.Data(), pt));
                        mGraphAN_xF[run][type][dle][i][pt] -> SetTitle("; x_{F}; A_{N}");

                        mGraphAN_Bkg_xF[run][type][dle][i][pt] = new TGraphErrors();
                        mGraphAN_Bkg_xF[run][type][dle][i][pt] -> SetName(Form("AN_Bkg_xF_%s_type%i_%s_ptBin%i", runName.Data(), type, dleName.Data(), pt));
                        mGraphAN_Bkg_xF[run][type][dle][i][pt] -> SetTitle("; x_{F}; A_{N}");


                        mGraphAN_Subt_xF[run][type][dle][i][pt] = new TGraphErrors();
                        mGraphAN_Subt_xF[run][type][dle][i][pt] -> SetName(Form("AN_Sut_xF_%s_type%i_%s_ptBin%i", runName.Data(), type, dleName.Data(), pt));
                        mGraphAN_Subt_xF[run][type][dle][i][pt] -> SetTitle("; x_{F}; A_{N}");
                  
                    }
                }
            }
        }

        for(int run=0; run<kRunNum; run++){
            for(int type=0; type<kTypeNum; type++){
                for(int i=0; i<5; i++){
                    mGraphANSysErr_Subt_pT[run][type][dle][i].resize(xfNum);
                    for(int xf=0; xf<xfNum; xf++){
                        mGraphANSysErr_Subt_pT[run][type][dle][i][xf] = new TGraphErrors();
                        mGraphANSysErr_Subt_pT[run][type][dle][i][xf] -> SetName(Form("AN_SysErr_Subt_pT_%i_type%i_%s_%i_xfBin%i", run, type, dleName.Data(), i, xf));
                    }
                    mGraphANSysErr_Subt_xF[run][type][dle][i].resize(ptNum);
                    for(int pt=0; pt<ptNum; pt++){
                        mGraphANSysErr_Subt_xF[run][type][dle][i][pt] = new TGraphErrors();
                        mGraphANSysErr_Subt_xF[run][type][dle][i][pt] -> SetName(Form("AN_SysErr_Sut_xF_%i_type%i_%s_%i_ptBin%i", run, type, dleName.Data(), i, pt));  
                    }
                }

                for(int i=0; i<4; i++){
                    for(int j=0; j<4; j++){
                        for(int k=0; k<2; k++){
                            mGraphANSysErr_Calcul_subt_pT[run][type][dle][i][j][k].resize(xfNum);
                            for(int xf=0; xf<xfNum; xf++){
                                mGraphANSysErr_Calcul_subt_pT[run][type][dle][i][j][k][xf] = new TGraphErrors();
                                mGraphANSysErr_Calcul_subt_pT[run][type][dle][i][j][k][xf] -> SetName(Form("ANStsErrCalcul_pT_%i_type%i_%s_xfBin%i_%i_%i_%i", run, type, dleName.Data(), xf, i, j, k));
                            }
                            mGraphANSysErr_Calcul_subt_xF[run][type][dle][i][j][k].resize(ptNum);
                            for(int pt=0; pt<ptNum; pt++){
                                mGraphANSysErr_Calcul_subt_xF[run][type][dle][i][j][k][pt] = new TGraphErrors();
                                mGraphANSysErr_Calcul_subt_xF[run][type][dle][i][j][k][pt] -> SetName(Form("ANStsErrCalcul_xF_%i_type%i_%s_ptBin%i_%i_%i_%i", run, type, dleName.Data(), pt, i, j, k));
                            }
                        }
                    }
                }
            }
        }


        mGraphAN_pTSummary[dle].resize(xfNum);
        for(int xf=0; xf<xfNum; xf++){
            mGraphAN_pTSummary[dle][xf] = new TGraphErrors();
            mGraphAN_pTSummary[dle][xf] -> SetName(Form("AN_pTSummary_%s_xfBin%i", dleName.Data(), xf));
            mGraphAN_pTSummary[dle][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
        }
        mGraphAN_xFSummary[dle].resize(ptNum);
        for(int pt=0; pt<ptNum; pt++){
            mGraphAN_xFSummary[dle][pt] = new TGraphErrors();
            mGraphAN_xFSummary[dle][pt] -> SetName(Form("AN_xFSummary_%s_ptBin%i", dleName.Data(), pt));
            mGraphAN_xFSummary[dle][pt] -> SetTitle("; x_{F}; A_{N}");
        }
        for(int s=0; s<5; s++){
            mGraphAN_pTSummarySysErr[dle][s].resize(xfNum);
            for(int xf=0; xf<xfNum; xf++){
                mGraphAN_pTSummarySysErr[dle][s][xf] = new TGraphErrors();
                mGraphAN_pTSummarySysErr[dle][s][xf] -> SetName(Form("AN_pTSummary_SysErr_%s_xfBin%i", dleName.Data(), xf));
                mGraphAN_pTSummarySysErr[dle][s][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
            }
            mGraphAN_xFSummarySysErr[dle][s].resize(ptNum);
            for(int pt=0; pt<ptNum; pt++){
                mGraphAN_xFSummarySysErr[dle][s][pt] = new TGraphErrors();
                mGraphAN_xFSummarySysErr[dle][s][pt] -> SetName(Form("AN_xFSummary_SysErr_%s_ptBin%i", dleName.Data(), pt));
                mGraphAN_xFSummarySysErr[dle][s][pt] -> SetTitle("; x_{F}; A_{N}");
            }
        }




        mGraphAN_pTSummary_merge[dle].resize(2);
        for(int xf=0; xf<2; xf++){
            mGraphAN_pTSummary_merge[dle][xf] = new TGraphErrors();
            mGraphAN_pTSummary_merge[dle][xf] -> SetName(Form("AN_pTSummary_merge_%s_xfBin%i", dleName.Data(), xf));
            mGraphAN_pTSummary_merge[dle][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
        }
        mGraphAN_xFSummary_merge[dle].resize(4);
        for(int pt=0; pt<4; pt++){
            mGraphAN_xFSummary_merge[dle][pt] = new TGraphErrors();
            mGraphAN_xFSummary_merge[dle][pt] -> SetName(Form("AN_xFSummary_merge_%s_ptBin%i", dleName.Data(), pt));
            mGraphAN_xFSummary_merge[dle][pt] -> SetTitle("; x_{F}; A_{N}");
        }

        for(int s=0; s<5; s++){
            mGraphAN_pTSummarySysErr_merge[dle][s].resize(2);
            for(int xf=0; xf<2; xf++){
                mGraphAN_pTSummarySysErr_merge[dle][s][xf] = new TGraphErrors();
                mGraphAN_pTSummarySysErr_merge[dle][s][xf] -> SetName(Form("AN_pTSummary_SysErr_merge_%s_xfBin%i", dleName.Data(), xf));
                mGraphAN_pTSummarySysErr_merge[dle][s][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
            }
            mGraphAN_xFSummarySysErr_merge[dle][s].resize(4);
            for(int pt=0; pt<4; pt++){
                mGraphAN_xFSummarySysErr_merge[dle][s][pt] = new TGraphErrors();
                mGraphAN_xFSummarySysErr_merge[dle][s][pt] -> SetName(Form("AN_xFSummary_SysErr_merge_%s_ptBin%i", dleName.Data(), pt));
                mGraphAN_xFSummarySysErr_merge[dle][s][pt] -> SetTitle("; x_{F}; A_{N}");
            }
        }

        if(dle < 3){
            mGraphAN_pTSummary_DLEmerge[dle].resize(2);
            for(int xf=0; xf<2; xf++){
                mGraphAN_pTSummary_DLEmerge[dle][xf] = new TGraphErrors();
                mGraphAN_pTSummary_DLEmerge[dle][xf] -> SetName(Form("AN_pTSummary_DLEmerge_%s_xfBin%i", dleName.Data(), xf));
                mGraphAN_pTSummary_DLEmerge[dle][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
            }
            mGraphAN_xFSummary_DLEmerge[dle].resize(4);
            for(int pt=0; pt<4; pt++){
                mGraphAN_xFSummary_DLEmerge[dle][pt] = new TGraphErrors();
                mGraphAN_xFSummary_DLEmerge[dle][pt] -> SetName(Form("AN_xFSummary_DLEmerge_%s_ptBin%i", dleName.Data(), pt));
                mGraphAN_xFSummary_DLEmerge[dle][pt] -> SetTitle("; x_{F}; A_{N}");
            }

            for(int s=0; s<5; s++){
                mGraphAN_pTSummarySysErr_DLEmerge[dle][s].resize(2);
                for(int xf=0; xf<2; xf++){
                    mGraphAN_pTSummarySysErr_DLEmerge[dle][s][xf] = new TGraphErrors();
                    mGraphAN_pTSummarySysErr_DLEmerge[dle][s][xf] -> SetName(Form("AN_pTSummary_SysErr_DLEmerge_%s_xfBin%i", dleName.Data(), xf));
                    mGraphAN_pTSummarySysErr_DLEmerge[dle][s][xf] -> SetTitle("; p_{T} [GeV/c]; A_{N}");
                }
                mGraphAN_xFSummarySysErr_DLEmerge[dle][s].resize(4);
                for(int pt=0; pt<4; pt++){
                    mGraphAN_xFSummarySysErr_DLEmerge[dle][s][pt] = new TGraphErrors();
                    mGraphAN_xFSummarySysErr_DLEmerge[dle][s][pt] -> SetName(Form("AN_xFSummary_SysErr_DLEmerge_%s_ptBin%i", dleName.Data(), pt));
                    mGraphAN_xFSummarySysErr_DLEmerge[dle][s][pt] -> SetTitle("; x_{F}; A_{N}");
                }
            }
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
    // if(fillIdx = 0){return 0.536;}
    // if(fillIdx = 1){return 0.554;}
    // if(fillIdx = 2){return 0.590;}
    // if(fillIdx = 3){return 0.566;}
    // if(fillIdx = 4){return 0.592;}
    if(fillIdx = 0){return 0.53978;}
    if(fillIdx = 1){return 0.60523;}
    if(fillIdx = 2){return 0.60171;}
    if(fillIdx = 3){return 0.61295;}
    if(fillIdx = 4){return 0.61571;}
    return 0.;
}

double RHICfAsymmetry::BeamPolarizationError(int fillIdx)
{
    if(fillIdx = 0){return 0.02;}
    if(fillIdx = 1){return 0.017;}
    if(fillIdx = 2){return 0.016;}
    if(fillIdx = 3){return 0.012;}
    if(fillIdx = 4){return 0.018;}
    return 0.;
}

