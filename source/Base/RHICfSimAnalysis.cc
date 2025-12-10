#include "RHICfSimAnalysis.hh"

#include "DrawingUtil.hh"


RHICfSimAnalysis::RHICfSimAnalysis()
{
}

RHICfSimAnalysis::~RHICfSimAnalysis()
{
}

int RHICfSimAnalysis::Init()
{
    if(!RHICfOptContainer::Init()){return 0;}

    mCalRunType = RHICfOptContainer::GetRunType();
    for(int model=0; model<mModelNum; model++){
        if(model==rEPOSLHCR_F || model==QGSJETII04){continue;}
        GetOptContainer()->SetSimModel(model);

        for(int run=0; run<kRunNum; run++){
            if(mCalRunType != kALLRun && mCalRunType != run){continue;}

            TString runName = GetOptContainer()->GetRunTypeName(run);
            GetOptContainer()->SetRunType(runName);

            mSimDstReader[model][run] = new RHICfSimDstReader();
            mSimDstReader[model][run] -> Init();
        }
    }

    mBinning = new RHICfBinning();
    mBinning -> Init();

    return 1;
}

int RHICfSimAnalysis::Calculate()
{
    const int BoundaryNum = 5;
    const int binNum = 4;
    const int dleNum = 4;

    // ========= Binning ==========
    double ptBoundary[BoundaryNum];
    double xfBoundary[BoundaryNum];
    double xfBoundaryForPtBin[2];
    double ptBoundaryForXfBin[2];

    xfBoundaryForPtBin[0] = mBinning -> GetGlobalXfBinBoundary(2);
    xfBoundaryForPtBin[1] = mBinning -> GetGlobalXfBinBoundary(4);
    ptBoundaryForXfBin[0] = mBinning -> GetGlobalPtBinBoundary(2);
    ptBoundaryForXfBin[1] = mBinning -> GetGlobalPtBinBoundary(4);

    for(int bin=0; bin<BoundaryNum; bin++){
        ptBoundary[bin] = mBinning -> GetGlobalPtBinBoundary(bin);
        xfBoundary[bin] = mBinning -> GetGlobalXfBinBoundary(bin);
    }

    TH1D* hDataMass[kRunNum][dleNum][2];
    TH1D* hDataMassPt[kRunNum][dleNum][binNum];
    TH1D* hDataMassXf[kRunNum][dleNum][binNum];
    TH1D* hDataPtXf[kRunNum][dleNum][2];
    TH1D* hDataPt[kRunNum][dleNum][2];
    TH1D* hDataXf[kRunNum][dleNum][2];
    TH1D* hDataE[kRunNum][dleNum][2];
    TH1D* hDataBBC[kRunNum][2][2];
    TH1D* hDataBTOF[kRunNum];

    TH1D* hSimMass[kRunNum][mModelNum][dleNum][2];
    TH1D* hSimMassPt[kRunNum][mModelNum][dleNum][binNum];
    TH1D* hSimMassXf[kRunNum][mModelNum][dleNum][binNum];
    TH1D* hSimPt[kRunNum][mModelNum][dleNum][2];
    TH1D* hSimXf[kRunNum][mModelNum][dleNum][2];
    TH1D* hSimPtXf[kRunNum][mModelNum][dleNum][2];
    TH1D* hSimE[kRunNum][mModelNum][dleNum][2];
    TH1D* hSimBBC[kRunNum][mModelNum][dleNum][2][2];
    TH1D* hSimBTOF[kRunNum][mModelNum][dleNum];

    for(int run=0; run<kRunNum; run++){
        for(int e=0; e<dleNum; e++){
            for(int b=0; b<binNum; b++){
                hDataMassPt[run][e][b] = new TH1D(Form("hDataMassPt_%i_%i_%i", run, e, b), "", 100, 0., 300.);
                hDataMassXf[run][e][b] = new TH1D(Form("hDataMassXf_%i_%i_%i", run, e, b), "", 100, 0., 300.);
                for(int m=0; m<mModelNum; m++){
                    hSimMassPt[run][m][e][b] = new TH1D(Form("hSimMassPt_%i_%i_%i_%i", run, m, e, b), "", 100, 0., 300.);
                    hSimMassXf[run][m][e][b] = new TH1D(Form("hSimMassXf_%i_%i_%i_%i", run, m, e, b), "", 100, 0., 300.);
                }
            }
        
            for(int type=0; type<2; type++){
                hDataMass[run][e][type] = new TH1D(Form("hDataMass_%i_%i_%i", run, e, type), "", 100, 0., 300.);
                hDataPt[run][e][type] = new TH1D(Form("hDataPt_%i_%i_%i", run, e, type), "", 60, 0., 1.);
                hDataXf[run][e][type] = new TH1D(Form("hDataXf_%i_%i_%i", run, e, type), "", 60, 0., 1.);
                hDataPtXf[run][e][type] = new TH1D(Form("hDataPtXf_%i_%i_%i", run, e, type), "", 60, 0., 1.);
                hDataE[run][e][type] = new TH1D(Form("hDataE_%i_%i_%i", run, e, type), "", 60, 0., 250.);

                for(int m=0; m<mModelNum; m++){
                    hSimMass[run][m][e][type] = new TH1D(Form("hSimMass_%i_%i_%i_%i", run, m, e, type), "", 100, 0., 300.);
                    hSimPt[run][m][e][type] = new TH1D(Form("hSimPt_%i_%i_%i_%i", run, m, e, type), "", 60, 0., 1.);
                    hSimXf[run][m][e][type] = new TH1D(Form("hSimXf_%i_%i_%i_%i", run, m, e, type), "", 60, 0., 1.);
                    hSimPtXf[run][m][e][type] = new TH1D(Form("hSimPtXf_%i_%i_%i_%i", run, m, e, type), "", 60, 0., 1.);
                    hSimE[run][m][e][type] = new TH1D(Form("hSimE_%i_%i_%i_%i", run, m, e, type), "", 60, 0., 270.);
                }
            }
        }
        hDataBTOF[run] = new TH1D(Form("hDataBTOF_%i",run), "", 60, 0, 60);
        for(int ew=0; ew<2; ew++){
            for(int sl=0; sl<2; sl++){
                hDataBBC[run][ew][sl] = new TH1D(Form("hDataBBC_%i_%i_%i", run, ew, sl), "", 70, 0., 600.);

                for(int m=0; m<mModelNum; m++){
                    for(int e=0; e<dleNum; e++){
                        if(ew==0 && sl==0){hSimBTOF[run][m][e] = new TH1D(Form("hSimBTOF_%i_%i_%i", run, m, e), "", 60, 0, 60);}
                        hSimBBC[run][m][e][ew][sl] = new TH1D(Form("hSimBBC_%i_%i_%i_%i_%i", run, m, e, ew, sl), "", 70, 0., 600.);
                    }
                }
            }
        }
    }

    TH2D* hPi0Pos[kRunNum];
    TH2D* hSimPi0Pos[kRunNum][mModelNum];

    for(int run=0; run<kRunNum; run++){
        hPi0Pos[run] = new TH2D(Form("hPi0Pos_%i", run), "", 100, -200., 200., 100, -200., 200.);
        for(int m=0; m<mModelNum; m++){
            hSimPi0Pos[run][m] = new TH2D(Form("hSimPi0Pos_%i_%i", run, m), "", 100, -50., 50., 100, 0., 220.);
        }
    }

    double nDLEData[kRunNum][3][dleNum];
    double nDLEDataPt[kRunNum][binNum][2][dleNum];
    double nDLEDataXf[kRunNum][binNum][2][dleNum];
    double nDLESimPt[kRunNum][mModelNum][binNum][2][dleNum][dleNum];
    double nDLESimXf[kRunNum][mModelNum][binNum][2][dleNum][dleNum];

    memset(nDLEData, 0, sizeof(nDLEData));
    memset(nDLEDataPt, 0, sizeof(nDLEDataPt));
    memset(nDLEDataXf, 0, sizeof(nDLEDataXf));
    memset(nDLESimPt, 0, sizeof(nDLESimPt));
    memset(nDLESimXf, 0, sizeof(nDLESimXf));

    double nDLEALL[kRunNum][mModelNum][dleNum][dleNum]; // [model][DLE][process] for all
    double nDLE[kRunNum][mModelNum][dleNum][dleNum]; // [model][DLE][process] for mass window events (120 - 160)
    double nDLE_bkg[kRunNum][mModelNum][dleNum][dleNum]; // [model][DLE][process] for bkg window events (~100)
    double nDLE_tSig[kRunNum][mModelNum][dleNum][dleNum]; // [model][DLE][process] for truth signal events in mass window
    double nDLE_tBkg[kRunNum][mModelNum][dleNum][dleNum]; // [model][DLE][process] for truth bkg events in bkg window

    memset(nDLEALL, 0, sizeof(nDLEALL));
    memset(nDLE, 0, sizeof(nDLE));
    memset(nDLE_bkg, 0, sizeof(nDLE_bkg));
    memset(nDLE_tSig, 0, sizeof(nDLE_tSig));
    memset(nDLE_tBkg, 0, sizeof(nDLE_tBkg));





    Int_t RunNum;
    Int_t FillNum;
    Int_t RHICfRunType;

    Int_t DLEIdx;
    Int_t BTofMult;
    Int_t BBCSumADC[kBeamSideNum][2];
    
    Int_t RHICfParticleNum;
    Int_t Pi0Type[4];
    Int_t TowerIdx[4];
    Double_t Pos[4][2];
    Double_t Momentum[4][4];
    Double_t Mass[4];

    Double_t RHICfPlatedE[2][16];
    Int_t ZDCADC[kBeamSideNum][kZDCPmtNum];
    Int_t ZDCTDC[kBeamSideNum][kZDCPmtNum];
    Int_t ZDCSumADC[kBeamSideNum][2]; // [Attenuated, Unattenuated]
    Int_t ZDCSMD[kBeamSideNum][2][kSMDNum];  //[east, west][vertical, horizontal][smd]

    TFile* file[kRunNum];
    TTree* tree[kRunNum];

    int eventNum[kRunNum];

    TString picoEventPath = "/gpfs01/star/pwg/slee5/RHICfPicoEvent";
    TString particleRunName = RHICfOptContainer::GetParticleRunName();

    for(int run=0; run<kRunNum; run++){
        if(mCalRunType != kALLRun && mCalRunType != run){continue;}
        TString RunName = RHICfOptContainer::GetRunTypeName(run);

        file[run] = new TFile(Form("%s/RHICfPicoEvent_%s_%s.root", picoEventPath.Data(), particleRunName.Data(), RunName.Data()), "read");
        tree[run] = (TTree*)file[run] -> Get("event");
        tree[run] -> SetBranchAddress("RunNum", &RunNum);
        tree[run] -> SetBranchAddress("FillNum", &FillNum);
        tree[run] -> SetBranchAddress("RHICfRunType", &RHICfRunType);
        tree[run] -> SetBranchAddress("DLEIdx", &DLEIdx);
        tree[run] -> SetBranchAddress("BTofMult", &BTofMult);
        tree[run] -> SetBranchAddress("BBCSumADC", &BBCSumADC);
        // tree[run] -> SetBranchAddress("particleNum", &RHICfParticleNum);
        tree[run] -> SetBranchAddress("Pi0Type", &Pi0Type);
        tree[run] -> SetBranchAddress("TowerIdx", &TowerIdx);
        tree[run] -> SetBranchAddress("Position", &Pos);
        tree[run] -> SetBranchAddress("Momentum", &Momentum);
        tree[run] -> SetBranchAddress("Mass", &Mass);
        // tree[run] -> SetBranchAddress("RHICfPlatedE", &RHICfPlatedE);
        tree[run] -> SetBranchAddress("ZDCADC", &ZDCADC);
        tree[run] -> SetBranchAddress("ZDCTDC", &ZDCTDC);
        tree[run] -> SetBranchAddress("ZDCSumADC", &ZDCSumADC);
        tree[run] -> SetBranchAddress("ZDCSMD", &ZDCSMD);

        eventNum[run] = tree[run] -> GetEntries();
        // eventNum[run] = 0;
        for(int event=0; event<eventNum[run]; event++){
            if(event%10000 == 0){cout << "RHICfSimAnalysis::Calculate() Real data -- "<< RHICfOptContainer::GetRunTypeName(run) << " Event: " << event << " / " << eventNum[run] << endl;}
            tree[run] -> GetEntry(event);

            if(Mass[0] < 5.){continue;}

            double px = Momentum[0][0];
            double py = Momentum[0][1];
            double pz = Momentum[0][2];
            double e = Momentum[0][3];
            double pt = sqrt(px*px + py*py);
            double xf = pz/255.;
            double x = Pos[0][0]*1000.;
            double y = Pos[0][1]*1000.;

            if(Pi0Type[0] == 1){hPi0Pos[run] -> Fill(x, y);}
            

            hDataMass[run][DLEIdx][Pi0Type[0]] -> Fill(Mass[0]);
            hDataMass[run][kALLDLE][Pi0Type[0]] -> Fill(Mass[0]);
            if(120 < Mass[0] && Mass[0] < 160.){
                hDataPt[run][DLEIdx][Pi0Type[0]] -> Fill(pt);
                hDataPt[run][kALLDLE][Pi0Type[0]] -> Fill(pt);
                hDataXf[run][DLEIdx][Pi0Type[0]] -> Fill(xf);
                hDataXf[run][kALLDLE][Pi0Type[0]] -> Fill(xf);
                hDataE[run][DLEIdx][Pi0Type[0]] -> Fill(e);
                hDataE[run][kALLDLE][Pi0Type[0]] -> Fill(e);
                hDataPtXf[run][DLEIdx][0] -> Fill(pt);
                hDataPtXf[run][kALLDLE][0] -> Fill(pt);
                hDataPtXf[run][DLEIdx][1] -> Fill(xf);
                hDataPtXf[run][kALLDLE][1] -> Fill(xf);

                nDLEData[run][1][DLEIdx] +=1.;
                nDLEData[run][1][kALLDLE] +=1.;
            }
            if(Mass[0] < 100){
                nDLEData[run][2][DLEIdx] +=1.;
                nDLEData[run][2][kALLDLE] +=1.;
            }
            nDLEData[run][0][DLEIdx] +=1.;
            nDLEData[run][0][kALLDLE] +=1.;

            for(int b=0; b<binNum; b++){
                if(xfBoundaryForPtBin[0] <= xf && xf <= xfBoundaryForPtBin[1]){
                    if(ptBoundary[b] <= pt && pt <= ptBoundary[b+1]){
                        hDataMassPt[run][DLEIdx][b] -> Fill(Mass[0]);
                        hDataMassPt[run][kALLDLE][b] -> Fill(Mass[0]);

                        if(120 < Mass[0] && Mass[0] < 160.){
                            nDLEDataPt[run][b][0][DLEIdx] += 1.;
                            nDLEDataPt[run][b][0][kALLDLE] += 1.;
                        }
                        if(Mass[0] < 100.){
                            nDLEDataPt[run][b][1][DLEIdx] += 1.;
                            nDLEDataPt[run][b][1][kALLDLE] += 1.;
                        }
                    }
                }
                if(ptBoundaryForXfBin[0] <= pt && pt <= ptBoundaryForXfBin[1]){
                    if(xfBoundary[b] <= xf && xf <= xfBoundary[b+1]){
                        hDataMassXf[run][DLEIdx][b] -> Fill(Mass[0]);
                        hDataMassXf[run][kALLDLE][b] -> Fill(Mass[0]);

                        if(120 < Mass[0] && Mass[0] < 160.){
                            nDLEDataXf[run][b][0][DLEIdx] += 1.;
                            nDLEDataXf[run][b][0][kALLDLE] += 1.;
                        }
                        if(Mass[0] < 100.){
                            nDLEDataXf[run][b][1][DLEIdx] += 1.;
                            nDLEDataXf[run][b][1][kALLDLE] += 1.;
                        }
                    }
                }
            }

            hDataBTOF[run] -> Fill(BTofMult);
            for(int ew=0; ew<2; ew++){
                for(int sl=0; sl<2; sl++){
                    hDataBBC[run][ew][sl] -> Fill(BBCSumADC[ew][sl]);
                }
            }
        }
    }

    for(int m=0; m<mModelNum; m++){
        if(m==rEPOSLHCR_F || m==QGSJETII04){continue;}

        int simEventNum[kRunNum];
        memset(simEventNum, 0, sizeof(simEventNum));
        for(int run=0; run<kRunNum; run++){
            if(mCalRunType != kALLRun && mCalRunType != run){continue;}
            simEventNum[run] = mSimDstReader[m][run]->GetEventNum();
        }
        for(int run=0; run<kRunNum; run++){
            if(mCalRunType != kALLRun && mCalRunType != run){continue;}
            for(int event=0; event<simEventNum[run]; event++){
                if(event%10000 == 0){cout << "RHICfSimAnalysis::Calculate() -- " << RHICfOptContainer::GetSimModelName(m) << " " << RHICfOptContainer::GetRunTypeName(run) << " Event: " << event << " / " << simEventNum[run] << endl;}
                mMiniSimDst = mSimDstReader[m][run] -> GetMiniSimDst(event);

                // Event information
                int runType = mMiniSimDst -> GetRHICfRunType();
                int eventType = mMiniSimDst -> GetEventType(); // [0 = SD, 1 = DD, 2 = ND]
                int DLEIdx = mMiniSimDst -> GetDLEIdx(); // [0 = SDLE, 1 = DDLE, 2 = NDLE]

                // RHICf particle
                int particleNum = mMiniSimDst -> GetParticleNum();
                if(particleNum != 1){continue;}

                int pid = mMiniSimDst -> GetPID(); // 0 == pi0, 1 == neutron
                int type = mMiniSimDst -> GetPi0Type();
                double px = mMiniSimDst -> GetPx();
                double py = mMiniSimDst -> GetPy();
                double pz = mMiniSimDst -> GetPz();
                double pt = sqrt(px*px + py*py);
                double xf = pz/255.;
                double e = mMiniSimDst -> GetEnergy();
                double mass = mMiniSimDst -> GetMass();

                double x = mMiniSimDst -> GetPosX()*1000.;
                double y = mMiniSimDst -> GetPosY()*1000.;

                if(type == 1){hSimPi0Pos[run][m] -> Fill(x, y);}
                
                if(pid == 1){continue;}

                int truthFlag = GetTruthParFlag(type); // [1 == signal, 0 == bkg]

                nDLEALL[run][m][DLEIdx][eventType] += 1.;
                nDLEALL[run][m][DLEIdx][dleNum-1] += 1.;
                nDLEALL[run][m][dleNum-1][eventType] += 1.;
                nDLEALL[run][m][dleNum-1][dleNum-1] += 1.;
                
                if(120. < mass && mass < 160.){
                    nDLE[run][m][DLEIdx][eventType] += 1.;
                    nDLE[run][m][DLEIdx][dleNum-1] += 1.;
                    nDLE[run][m][dleNum-1][eventType] += 1.;
                    nDLE[run][m][dleNum-1][dleNum-1] += 1.;
                    if(truthFlag == 1){
                        nDLE_tSig[run][m][DLEIdx][eventType] += 1.;
                        nDLE_tSig[run][m][DLEIdx][dleNum-1] += 1.;
                        nDLE_tSig[run][m][dleNum-1][eventType] += 1.;
                        nDLE_tSig[run][m][dleNum-1][dleNum-1] += 1.;
                    }
                }
                if(mass < 100.){
                    nDLE_bkg[run][m][DLEIdx][eventType] += 1.;
                    nDLE_bkg[run][m][DLEIdx][dleNum-1] += 1.;
                    nDLE_bkg[run][m][dleNum-1][eventType] += 1.;
                    nDLE_bkg[run][m][dleNum-1][dleNum-1] += 1.;

                    if(truthFlag == 0){
                        nDLE_tBkg[run][m][DLEIdx][eventType] += 1.;
                        nDLE_tBkg[run][m][DLEIdx][dleNum-1] += 1.;
                        nDLE_tBkg[run][m][dleNum-1][eventType] += 1.;
                        nDLE_tBkg[run][m][dleNum-1][dleNum-1] += 1.;
                    }
                }

                // =============================
                // =============================
                // =============================

                hSimMass[run][m][DLEIdx][type] -> Fill(mass);
                hSimMass[run][m][kALLDLE][type] -> Fill(mass);
                if(120. < mass  && mass < 160.){
                    hSimPt[run][m][DLEIdx][type] -> Fill(pt);
                    hSimPt[run][m][kALLDLE][type] -> Fill(pt);
                    hSimXf[run][m][DLEIdx][type] -> Fill(xf);
                    hSimXf[run][m][kALLDLE][type] -> Fill(xf);
                    hSimE[run][m][DLEIdx][type] -> Fill(e);
                    hSimE[run][m][kALLDLE][type] -> Fill(e);
                    hSimPtXf[run][m][DLEIdx][0] -> Fill(pt);
                    hSimPtXf[run][m][kALLDLE][0] -> Fill(pt);
                    hSimPtXf[run][m][DLEIdx][1] -> Fill(xf);
                    hSimPtXf[run][m][kALLDLE][1] -> Fill(xf);
                }

                for(int b=0; b<binNum; b++){
                    if(xfBoundaryForPtBin[0] <= xf && xf <= xfBoundaryForPtBin[1]){
                        if(ptBoundary[b] <= pt && pt <= ptBoundary[b+1]){
                            hSimMassPt[run][m][DLEIdx][b] -> Fill(mass);
                            hSimMassPt[run][m][kALLDLE][b] -> Fill(mass);
                            
                            if(120. < mass && mass < 160.){
                                nDLESimPt[run][m][b][0][DLEIdx][eventType] += 1.;
                                nDLESimPt[run][m][b][0][kALLDLE][eventType] += 1.;
                            }
                            if(mass < 100.){
                                nDLESimPt[run][m][b][1][DLEIdx][eventType] += 1.;
                                nDLESimPt[run][m][b][1][kALLDLE][eventType] += 1.;
                            }
                        }
                    }
                    if(ptBoundaryForXfBin[0] <= pt && pt <= ptBoundaryForXfBin[1]){
                        if(xfBoundary[b] <= xf && xf <= xfBoundary[b+1]){
                            hSimMassXf[run][m][DLEIdx][b] -> Fill(mass);
                            hSimMassXf[run][m][kALLDLE][b] -> Fill(mass);

                            if(120. < mass && mass < 160.){
                                nDLESimXf[run][m][b][0][DLEIdx][eventType] += 1.;
                                nDLESimXf[run][m][b][0][kALLDLE][eventType] += 1.;
                            }
                            if(mass < 100.){
                                nDLESimXf[run][m][b][1][DLEIdx][eventType] += 1.;
                                nDLESimXf[run][m][b][1][kALLDLE][eventType] += 1.;
                            }
                        }
                    }
                }

                hSimBTOF[run][m][DLEIdx] -> Fill(mMiniSimDst->GetBTofMult());
                hSimBTOF[run][m][kALLDLE] -> Fill(mMiniSimDst->GetBTofMult());
                for(int ew=0; ew<2; ew++){
                    for(int sl=0; sl<2; sl++){
                        hSimBBC[run][m][DLEIdx][ew][sl] -> Fill(mMiniSimDst->GetBBCSumADC(ew, sl));
                        hSimBBC[run][m][kALLDLE][ew][sl] -> Fill(mMiniSimDst->GetBBCSumADC(ew, sl));
                    }
                }
            }
        }
    }
    // ==============================================================================================
    // ==============================================================================================
    // ==============================================================================================

    DrawingUtil* drawUtil = new DrawingUtil("../simFigure");

    for(int run=0; run<kRunNum; run++){
        drawUtil -> Clear("all");
        
        drawUtil -> SetText(true, "p+p#rightarrow#pi^{0}+X @ #sqrt{s} = 510 GeV");
        drawUtil -> SetText(true, Form("%s run", RHICfOptContainer::GetRunTypeName(run).Data()));

        drawUtil -> SetLegend(true, hDataMass[run][kALLDLE][0], "Data");
        for(int m=0; m<mModelNum; m++){
            if(hSimMass[run][m][kALLDLE][0]->GetEntries() < 1.){continue;}
            drawUtil -> SetLegend(true, hSimMass[run][m][kALLDLE][0], Form("%s", RHICfOptContainer::GetSimModelName(m).Data()));
        }

        drawUtil -> SetCanvas(2);
        for(int type=0; type<2; type++){
            drawUtil -> SetText(0, Form("Type%i #pi^{0} mass", type+1), type);
            drawUtil -> SetTH1D(hDataMass[run][kALLDLE][type], ";M_{#gamma#gamma} [GeV/c^{2}];", type, kBlack);
            for(int m=0; m<mModelNum; m++){
                drawUtil -> SetTH1D(hSimMass[run][m][kALLDLE][type], "", type);
            }
        }
        drawUtil -> DrawHistWithRatio(0, "n", "MC/Data");
        drawUtil -> SaveFigure(Form("Mass_%s", RHICfOptContainer::GetRunTypeName(run).Data()));

        // ==============================================================================================
        // ==============================================================================================
        // ==============================================================================================

        drawUtil -> Clear();
        drawUtil -> SetCanvas(binNum*2);
        int binNumCIdx = 0;
        for(int b=0; b<binNum; b++){
            drawUtil -> SetText(0, Form("%.2f<p_{T}<%.2f", ptBoundary[b], ptBoundary[b+1]), binNumCIdx);
            drawUtil -> SetText(0, Form("%.2f<x_{F}<%.2f", xfBoundaryForPtBin[0], xfBoundaryForPtBin[1]), binNumCIdx);
            drawUtil -> SetTH1D(hDataMassPt[run][kALLDLE][b], ";M_{#gamma#gamma} [GeV/c^{2}];", binNumCIdx, kBlack);
            for(int m=0; m<mModelNum; m++){
                drawUtil -> SetTH1D(hSimMassPt[run][m][kALLDLE][b], "", binNumCIdx);
            }
            binNumCIdx++;
        }
        for(int b=0; b<binNum; b++){
            drawUtil -> SetText(0, Form("%.2f<p_{T}<%.2f", ptBoundaryForXfBin[0], ptBoundaryForXfBin[1]), binNumCIdx);
            drawUtil -> SetText(0, Form("%.2f<x_{F}<%.2f", xfBoundary[b], xfBoundary[b+1]), binNumCIdx);
            drawUtil -> SetTH1D(hDataMassXf[run][kALLDLE][b], ";M_{#gamma#gamma} [GeV/c^{2}];", binNumCIdx, kBlack);
            for(int m=0; m<mModelNum; m++){
                drawUtil -> SetTH1D(hSimMassXf[run][m][kALLDLE][b], "", binNumCIdx);
            }
            binNumCIdx++;
        }
        drawUtil -> DrawHistWithRatio(0, "n", "MC/Data");
        drawUtil -> SaveFigure(Form("MassPtXf_%s", RHICfOptContainer::GetRunTypeName(run).Data()));

        // ==============================================================================================
        // ==============================================================================================
        // ==============================================================================================

        drawUtil -> Clear();
        drawUtil -> SetCanvas(8);
        for(int type=0; type<2; type++){
            drawUtil -> SetText(0, Form("Type%i #pi^{0}", type+1), type);
            drawUtil -> SetTH1D(hDataPt[run][kALLDLE][type], ";p_{T} [GeV/c];", type, kBlack);
            
            drawUtil -> SetText(0, Form("Type%i #pi^{0}", type+1), type+3);
            drawUtil -> SetTH1D(hDataXf[run][kALLDLE][type], ";x_{F};", type+3, kBlack);

            drawUtil -> SetText(0, Form("Type%i #pi^{0}", type+1), type+6);
            drawUtil -> SetTH1D(hDataE[run][kALLDLE][type], ";E [GeV];", type+6, kBlack);

            for(int m=0; m<mModelNum; m++){
                drawUtil -> SetTH1D(hSimPt[run][m][kALLDLE][type], "", type);
                drawUtil -> SetTH1D(hSimXf[run][m][kALLDLE][type], "", type+3);
                drawUtil -> SetTH1D(hSimE[run][m][kALLDLE][type], "", type+6);
            }
        }
        drawUtil -> SetTH1D(hDataPtXf[run][kALLDLE][0], ";p_{T} [GeV/c];", 2, kBlack);
        for(int m=0; m<mModelNum; m++){
            drawUtil -> SetTH1D(hSimPtXf[run][m][kALLDLE][0], "", 2);
        }
        drawUtil -> SetTH1D(hDataPtXf[run][kALLDLE][1], ";x_{F};", 5, kBlack);
        for(int m=0; m<mModelNum; m++){
            drawUtil -> SetTH1D(hSimPtXf[run][m][kALLDLE][1], "", 5);
        }
        drawUtil -> DrawHistWithRatio(0, "n", "MC/Data");
        drawUtil -> SaveFigure(Form("Pi0Kinematic_%s", RHICfOptContainer::GetRunTypeName(run).Data()));

        // ==============================================================================================
        // ==============================================================================================
        // ==============================================================================================
        drawUtil -> Clear();
        drawUtil -> SetCanvas(5);
        drawUtil -> SetText(0,"B-TOF Mult.", 0);
        drawUtil -> SetTH1D(hDataBTOF[run], ";(Mult);", 0, kBlack);
        for(int m=0; m<mModelNum; m++){
            drawUtil -> SetTH1D(hSimBTOF[run][m][kALLDLE], "", 0);
        }
        int tmpIdx = 1;
        for(int ew=0; ew<2; ew++){
            TString ewName = (ew==0)? "East" : "West";
            for(int sl=0; sl<2; sl++){
                TString slName = (sl==0)? "Small" : "Large";

                drawUtil -> SetText(0,Form("%s BBC %s", ewName.Data(), slName.Data()), tmpIdx);
                drawUtil -> SetTH1D(hDataBBC[run][ew][sl], ";(ADC);", tmpIdx, kBlack);
                for(int m=0; m<mModelNum; m++){
                    drawUtil -> SetTH1D(hSimBBC[run][m][kALLDLE][ew][sl], "", tmpIdx);
                }
                tmpIdx++;
            }
        }
        drawUtil -> DrawHistWithRatio(0, "n ylog", "MC/Data");
        drawUtil -> SaveFigure(Form("STARDet_%s", RHICfOptContainer::GetRunTypeName(run).Data()));


        // ==============================================================================================
        // ==============================================================================================
        // ==============================================================================================
        drawUtil -> Clear("all");
        drawUtil -> SetCanvas(12);
        drawUtil -> SetXLabelName(0, "SD");
        drawUtil -> SetXLabelName(1, "DD");
        drawUtil -> SetXLabelName(2, "ND");
        drawUtil -> SetText(true, "p+p#rightarrow#pi^{0}+X @ #sqrt{s} = 510 GeV", -1, -1, -1, 0.035);
        drawUtil -> SetText(true, Form("%s run", RHICfOptContainer::GetRunTypeName(run).Data()), -1, -1, -1, 0.045);

        int tmpCIdx = 0;
        for(int dle=0; dle<dleNum-1; dle++){
            for(int m=0; m<mModelNum; m++){
                drawUtil -> SetTH1DRatio(dleNum-1, nDLE[run][m][dle], nDLE[run][m][dle][0]+nDLE[run][m][dle][1]+nDLE[run][m][dle][2], ";Process;Truth Process Ratio", tmpCIdx);
                drawUtil -> SetTH1DRatio(dleNum-1, nDLE_bkg[run][m][dle], nDLE_bkg[run][m][dle][0]+nDLE_bkg[run][m][dle][1]+nDLE_bkg[run][m][dle][2], "", tmpCIdx);
                drawUtil -> SetTH1DRatio(dleNum-1, nDLE_tSig[run][m][dle], nDLE_tSig[run][m][dle][0]+nDLE_tSig[run][m][dle][1]+nDLE_tSig[run][m][dle][2], "", tmpCIdx);
                drawUtil -> SetTH1DRatio(dleNum-1, nDLE_tBkg[run][m][dle], nDLE_tBkg[run][m][dle][0]+nDLE_tBkg[run][m][dle][1]+nDLE_tBkg[run][m][dle][2], "", tmpCIdx);

                drawUtil -> SetText(false, Form("%s", RHICfOptContainer::GetSimModelName(m).Data()), tmpCIdx);
                drawUtil -> SetText(false, Form("%s", RHICfOptContainer::GetDLEName(dle).Data()), tmpCIdx);

                if(tmpCIdx==0){
                    cout << drawUtil->GetTH1D(0, 0) -> GetName() << endl;
                    drawUtil -> SetLegend(true, drawUtil->GetTH1D(0, 0), "Signal region");
                    drawUtil -> SetLegend(true, drawUtil->GetTH1D(0, 1), "Bkg. region");
                    drawUtil -> SetLegend(true, drawUtil->GetTH1D(0, 2), "Truth Signal");
                    drawUtil -> SetLegend(true, drawUtil->GetTH1D(0, 3), "Truth Bkg.");
                }
                tmpCIdx++;
            }
        }
        drawUtil -> DrawHist();
        drawUtil -> SaveFigure(Form("dleRatio_massWindow_%s", RHICfOptContainer::GetRunTypeName(run).Data()));

        // ==============================================================================================
        // ==============================================================================================
        // ==============================================================================================
        drawUtil -> Clear("all");
        drawUtil -> SetCanvas(5);
        drawUtil -> SetText(true, "p+p#rightarrow#pi^{0}+X @ #sqrt{s} = 510 GeV", -1, -1, -1, 0.035);
        drawUtil -> SetText(true, Form("%s run", RHICfOptContainer::GetRunTypeName(run).Data()), -1, -1, -1, 0.045);
        drawUtil -> SetLegend(true, hDataMass[run][kALLDLE][0], "Data");
        for(int m=0; m<mModelNum; m++){
            if(hSimMass[run][m][kALLDLE][0]->GetEntries() < 1.){continue;}
            drawUtil -> SetLegend(true, hSimMass[run][m][kALLDLE][0], Form("%s", RHICfOptContainer::GetSimModelName(m).Data()));
        }

        drawUtil -> SetXLabelName(0, "SDLE");
        drawUtil -> SetXLabelName(1, "DDLE");
        drawUtil -> SetXLabelName(2, "NDLE");

        drawUtil -> SetText(false, "All region", 0);
        drawUtil -> SetText(false, "120 < M_{#gamma#gamma} < 160 GeV/c^{2}", 1);
        drawUtil -> SetText(false, "120 < M_{#gamma#gamma} < 160 GeV/c^{2}", 2);
        drawUtil -> SetText(false, "M_{#gamma#gamma} < 100 GeV/c^{2}", 3);
        drawUtil -> SetText(false, "M_{#gamma#gamma} < 100 GeV/c^{2}", 4);

        drawUtil -> SetTH1DRatio(dleNum-1, nDLEData[run][0], nDLEData[run][0][3], ";DLE;", 0, kBlack);
        drawUtil -> SetTH1DRatio(dleNum-1, nDLEData[run][1], nDLEData[run][1][3], ";DLE;", 1, kBlack);
        drawUtil -> SetTH1DRatio(dleNum-1, nDLEData[run][1], nDLEData[run][1][3], ";DLE;", 2, kBlack);
        drawUtil -> SetTH1DRatio(dleNum-1, nDLEData[run][2], nDLEData[run][2][3], ";DLE;", 3, kBlack);
        drawUtil -> SetTH1DRatio(dleNum-1, nDLEData[run][2], nDLEData[run][2][3], ";DLE;", 4, kBlack);
        for(int i=0; i<5; i++){
            drawUtil->GetTH1D(i, 0)->SetLineWidth(2);
        }

        for(int m=0; m<mModelNum; m++){
            double dleall[3];
            double dleSig[3];
            double dleTSig[3];
            double dleBkg[3];
            double dleTBkg[3];
            for(int i=0; i<3; i++){
                dleall[i] = nDLEALL[run][m][i][3];
                dleSig[i] = nDLE[run][m][i][3];
                dleTSig[i] = nDLE_tSig[run][m][i][3];
                dleBkg[i] = nDLE_bkg[run][m][i][3];
                dleTBkg[i] = nDLE_tBkg[run][m][i][3];
            }
            drawUtil -> SetTH1DRatio(dleNum-1, dleall, nDLEALL[run][m][3][3], ";DLE;", 0);
            drawUtil -> SetTH1DRatio(dleNum-1, dleSig, nDLE[run][m][3][3], ";DLE;", 1);
            drawUtil -> SetTH1DRatio(dleNum-1, dleTSig, nDLE_tSig[run][m][3][3], ";DLE;", 2);
            drawUtil -> SetTH1DRatio(dleNum-1, dleBkg, nDLE_bkg[run][m][3][3], ";DLE;", 3);
            drawUtil -> SetTH1DRatio(dleNum-1, dleTBkg, nDLE_tBkg[run][m][3][3], ";DLE;", 4);
        }

        drawUtil -> DrawHist();
        drawUtil -> SaveFigure(Form("dleRatio_massWindow2_%s", RHICfOptContainer::GetRunTypeName(run).Data()));

    

        // ==============================================================================================
        // ==============================================================================================
        // ==============================================================================================
        drawUtil -> Clear("leg");   
        drawUtil -> SetCanvas(6);
        drawUtil -> SetText(false, "120 < M_{#gamma#gamma} < 160 GeV/c^{2}");
        double x[binNum];
        double y[binNum];
        double total[binNum];
        for(int dle=0; dle<3; dle++){
            drawUtil -> SetText(false, Form("%s", RHICfOptContainer::GetDLEName(dle).Data()), dle);
            for(int m=0; m<mModelNum; m++){
                for(int bin=0; bin<binNum; bin++){
                    x[bin] = fabs(ptBoundary[bin] - ptBoundary[bin+1])/2. + ptBoundary[bin];
                    total[bin] = nDLESimPt[run][m][bin][0][dle][0]+nDLESimPt[run][m][bin][0][dle][1]+nDLESimPt[run][m][bin][0][dle][2];
                    y[bin] = nDLESimPt[run][m][bin][0][dle][dle];
                }
                drawUtil -> SetGraphRatio(binNum, x, y, total, ";p_{T} [GeV/c];Intended Process Ratio", dle);
            }
        }
        for(int dle=0; dle<3; dle++){
            drawUtil -> SetText(false, Form("%s", RHICfOptContainer::GetDLEName(dle).Data()), dle+3);
            for(int m=0; m<mModelNum; m++){
                for(int bin=0; bin<binNum; bin++){
                    x[bin] = fabs(xfBoundary[bin] - xfBoundary[bin+1])/2. + xfBoundary[bin];
                    total[bin] = nDLESimXf[run][m][bin][0][dle][0]+nDLESimXf[run][m][bin][0][dle][1]+nDLESimXf[run][m][bin][0][dle][2];
                    y[bin] = nDLESimXf[run][m][bin][0][dle][dle];
                }
                drawUtil -> SetGraphRatio(binNum, x, y, total, ";x_{F};Intended Process Ratio", dle+3);
            }
        }
        for(int m=0; m<mModelNum; m++){
            drawUtil -> SetLegend(true, drawUtil->GetGraph(0, m), Form("%s", RHICfOptContainer::GetSimModelName(m).Data()));
        }

        drawUtil -> DrawGraph("ratio");
        drawUtil -> SaveFigure(Form("DLE_Summary_massWindow_%s", RHICfOptContainer::GetRunTypeName(run).Data()));



        drawUtil -> Clear();   
        drawUtil -> SetCanvas(6);
        drawUtil -> SetText(false, "M_{#gamma#gamma} < 100 GeV/c^{2}");
        for(int dle=0; dle<3; dle++){
            drawUtil -> SetText(false, Form("%s", RHICfOptContainer::GetDLEName(dle).Data()), dle);
            for(int m=0; m<mModelNum; m++){
                for(int bin=0; bin<binNum; bin++){
                    x[bin] = fabs(ptBoundary[bin] - ptBoundary[bin+1])/2. + ptBoundary[bin];
                    total[bin] = nDLESimPt[run][m][bin][1][dle][0]+nDLESimPt[run][m][bin][1][dle][1]+nDLESimPt[run][m][bin][1][dle][2];
                    y[bin] = nDLESimPt[run][m][bin][1][dle][dle];
                }
                drawUtil -> SetGraphRatio(binNum, x, y, total, ";p_{T} [GeV/c];Intended Process Ratio", dle);
            }
        }
        for(int dle=0; dle<3; dle++){
            drawUtil -> SetText(false, Form("%s", RHICfOptContainer::GetDLEName(dle).Data()), dle+3);
            for(int m=0; m<mModelNum; m++){
                for(int bin=0; bin<binNum; bin++){
                    x[bin] = fabs(xfBoundary[bin] - xfBoundary[bin+1])/2. + xfBoundary[bin];
                    total[bin] = nDLESimXf[run][m][bin][1][dle][0]+nDLESimXf[run][m][bin][1][dle][1]+nDLESimXf[run][m][bin][1][dle][2];
                    y[bin] = nDLESimXf[run][m][bin][1][dle][dle];
                }
                drawUtil -> SetGraphRatio(binNum, x, y, total, ";x_{F};Intended Process Ratio", dle+3);
            }
        }
        drawUtil -> DrawGraph("ratio");
        drawUtil -> SaveFigure(Form("DLE_Summary_massBkg_%i", RHICfOptContainer::GetRunTypeName(run).Data()));

    }

    TCanvas* cPos = new TCanvas("cPos", "", 600.*3., 600.*2.);
    for(int run=0; run<3; run++){
        cPos -> Clear();
        cPos -> Divide(3,2);
        cPos -> cd(1);
        hPi0Pos[run] -> SetTitle("data;x;y");
        hPi0Pos[run] -> Draw("colz");

        for(int i=0; i<4; i++){
            cPos->cd(i+2);
            hSimPi0Pos[run][i] -> SetTitle(Form("%s;x;y", RHICfOptContainer::GetSimModelName(i).Data()));
            hSimPi0Pos[run][i] -> Draw("colz");
        }
        cPos -> Draw();
        cPos -> SaveAs(Form("../simFigure/Pi0Pos_%s.pdf", RHICfOptContainer::GetRunTypeName(run).Data()));
    }

    return 1;
}

int RHICfSimAnalysis::Finish()
{
    return 1;
}

int RHICfSimAnalysis::GetTruthParFlag(int type)
{
    vector<pair<int, int> > hitParticle[rTowerNum];
    vector<int> hitParNum[rTowerNum];

    // Loop for RHICf incident truth tracks
    for(int it=0; it<rTowerNum; it++){
        int truthNum = mMiniSimDst -> GetRHICfTruthNum(it);
        for(int truth=0; truth<truthNum; truth++){
            int truthIdx = mMiniSimDst -> GetRHICfTruthId(it, truth);

            int pid = mMiniSimDst -> GetSimTrackPid(truthIdx);
            double energy = mMiniSimDst -> GetSimTrackEnergy(truthIdx);
            int parentTrkId = mMiniSimDst -> GetSimTrackParentId(truthIdx);

            int parentPid = -1;
            // Find a track's origin of primary particles
            if(parentTrkId !=-1){
                parentPid = mMiniSimDst -> GetSimTrackPid(parentTrkId);

                while(true){
                    if(parentPid == 111){break;}
                    int grandParentTrkId = mMiniSimDst -> GetSimTrackParentId(parentTrkId);
                    if(grandParentTrkId == -1){break;}
                    int grandParentPid = mMiniSimDst -> GetSimTrackPid(grandParentTrkId);
                    parentTrkId = grandParentTrkId;
                    parentPid = grandParentPid;
                }
            }
            else{
                parentTrkId = truthIdx;
                parentPid = pid;
            }

            // Save the RHICf incident primary particles, and removed duplicated tracks
            bool isSaved = false;
            for(int i=0; i<hitParticle[it].size(); i++){
                if(hitParticle[it][i].first == parentTrkId){
                    isSaved = true; 
                    hitParNum[it][i] = hitParNum[it][i] + 1;
                    break;
                }
            }
            if(!isSaved){
                hitParticle[it].push_back(make_pair(parentTrkId, parentPid));
                hitParNum[it].push_back(1);
            }
        }
    }

    if(type == 0){
        for(int i=0; i<hitParticle[0].size(); i++){
            int tsID = hitParticle[0][i].first;
            int tsPID = hitParticle[0][i].second;

            if(tsPID == 111){
                for(int j=0; j<hitParticle[1].size(); j++){
                    int tlID = hitParticle[1][j].first;
                    int tlPID = hitParticle[1][j].second;

                    if(tlPID == 111){
                        if(tsID == tlID){return 1;}
                    }
                }
            }
        }
    }
    if(type == 1){
        for(int i=0; i<hitParticle[1].size(); i++){
            int pid = hitParticle[1][i].second;
            int num =  hitParNum[1][i];
            if(pid == 111 && num > 1){
                return 1;
            }
        }
    }

    return 0;
}
