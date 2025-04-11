#include "RHICfFigureDrawing.hh"

RHICfFigureDrawing::RHICfFigureDrawing() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
}

RHICfFigureDrawing::~RHICfFigureDrawing()
{
}

void RHICfFigureDrawing::Init()
{   
    cout << "RHICfFigureDrawing::Init() -- Done." << endl;
}

void RHICfFigureDrawing::SetBinning(RHICfBinning* binning){mBinning = binning;}
void RHICfFigureDrawing::SetMassFitting(RHICfMassFitting* fitting){mFitting = fitting;}
void RHICfFigureDrawing::SetDilution(RHICfDilutionFactor* dilution){mDilution = dilution;}
void RHICfFigureDrawing::SetAsymmetry(RHICfAsymmetry* asymmetry){mAsymmetry = asymmetry;}

void RHICfFigureDrawing::DrawMassHist()
{
    TString figurePath = mOptContainer -> GetFigurePath() + "/MassFitting";
    TString ParticleTypeName = mOptContainer -> GetParticleRunName();
    TCanvas* cMass = new TCanvas("cMass", "", 2400, 1200);

    TLatex* latex = new TLatex();

    int runIdx = mOptContainer->GetRunType();
    for(int run=0; run<kRunNum; run++){
        if(runIdx != kALLRun && runIdx != run){continue;}

        TString runName = mOptContainer -> GetRunTypeName(run);
        cMass -> Clear();
        cMass -> Divide(4,2);
        int tmpcIdx = 0;
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer->GetDLEName(dle);
                const int cIdx = tmpcIdx+1;
                cMass -> cd(cIdx);

                TH1D* mass = mFitting -> GetMassHist(run, type, dle);
                mass -> GetYaxis()->SetRangeUser(0, 1.3*mass->GetBinContent(mass->GetMaximumBin()));
                mass -> Draw();

                latex -> DrawLatexNDC(0.15, 0.85, GetSystemString());
                if(mOptContainer->GetParticleRunIdx() == kPi0Run){latex -> DrawLatexNDC(0.15, 0.79, Form("#pi^{0} type %i", type+1));}
                latex -> DrawLatexNDC(0.15, 0.73, Form("%s", dleName.Data()));

                latex -> DrawLatexNDC(0.53, 0.65, Form("<M> = %.2f", mFitting->GetMassMean(run, type, dle)));
                latex -> DrawLatexNDC(0.53, 0.59, Form("#sigma = %.2f", mFitting->GetMassSigma(run, type, dle)));
                latex -> DrawLatexNDC(0.53, 0.53, Form("All N = %i", int(mFitting->GetMassAllCounts(run, type, dle))));
                latex -> DrawLatexNDC(0.53, 0.47, Form("Sig. N = %i", int(mFitting->GetMassSignalCounts(run, type, dle))));
                latex -> DrawLatexNDC(0.53, 0.42, Form("Bkg. N = %i", int(mFitting->GetMassBkgCounts(run, type, dle))));

                tmpcIdx++;
            }
        }
        cMass -> Draw();
        cMass -> SaveAs(Form("%s/Mass_%s_%s.pdf", figurePath.Data(), ParticleTypeName.Data(), runName.Data()));
    }

    TCanvas* cMass_kinematics;
    for(int run=0; run<kRunNum; run++){
        if(runIdx != kALLRun && runIdx != run){continue;}
        TString runName = mOptContainer -> GetRunTypeName(run);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer->GetDLEName(dle);

                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                if(ptNum <= 0){continue;}
                
                cMass_kinematics = new TCanvas(Form("cMass_kinematics_%i_%i_%i", run, type, dle), "", double(xfNum*700), double(ptNum*700));
                cMass_kinematics -> Clear();
                cMass_kinematics -> Divide(xfNum, ptNum);

                int tmpcIdx = 0;
                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        const int cIdx = tmpcIdx+1;
                        cMass_kinematics -> cd(cIdx);

                        TH1D* mass = mFitting -> GetMassHist(run, type, dle, pt, xf);
                        mass -> GetYaxis()->SetRangeUser(0, 1.3*mass->GetBinContent(mass->GetMaximumBin()));
                        mass -> Draw();

                        latex -> DrawLatexNDC(0.15, 0.85, GetSystemString());
                        if(mOptContainer->GetParticleRunIdx() == kPi0Run){latex -> DrawLatexNDC(0.15, 0.79, Form("#pi^{0} type %i", type+1));}
                        latex -> DrawLatexNDC(0.15, 0.73, Form("%s", dleName.Data()));

                        latex -> DrawLatexNDC(0.5, 0.79, Form("%.2f < p_{T} < %.2f", mBinning->GetPtBinBoundary(run, type, dle, pt), mBinning->GetPtBinBoundary(run, type, dle, pt+1) ));
                        latex -> DrawLatexNDC(0.5, 0.73, Form("%.2f < x_{F} < %.2f", mBinning->GetXfBinBoundary(run, type, dle, xf), mBinning->GetXfBinBoundary(run, type, dle, xf+1) ));

                        latex -> DrawLatexNDC(0.53, 0.65, Form("<M> = %.2f", mFitting->GetMassMean(run, type, dle, pt, xf)));
                        latex -> DrawLatexNDC(0.53, 0.59, Form("#sigma = %.2f", mFitting->GetMassSigma(run, type, dle, pt, xf)));
                        latex -> DrawLatexNDC(0.53, 0.53, Form("All N = %i", int(mFitting->GetMassAllCounts(run, type, dle, pt, xf))));
                        latex -> DrawLatexNDC(0.53, 0.47, Form("Sig. N = %i", int(mFitting->GetMassSignalCounts(run, type, dle, pt, xf))));
                        latex -> DrawLatexNDC(0.53, 0.42, Form("Bkg. N = %i", int(mFitting->GetMassBkgCounts(run, type, dle, pt, xf))));

                        tmpcIdx++;
                    }
                }
                cMass_kinematics -> Draw();
                cMass_kinematics -> SaveAs(Form("%s/Mass_%s_%s_type%i_%s_kinematics.pdf", figurePath.Data(), ParticleTypeName.Data(), runName.Data(), type+1, dleName.Data()));
            }
        }
    }
}

void RHICfFigureDrawing::DrawBinningHist()
{
    TString figurePath = mOptContainer -> GetFigurePath() + "/Binning";
    TString ParticleTypeName = mOptContainer -> GetParticleRunName();
    TCanvas* cBinning = new TCanvas("cBinning", "", 2200., 600.);
    TLatex* latex = new TLatex();

    TGraph* graph = new TGraph();
    graph -> SetMarkerColor(kBlack);
    graph -> SetMarkerSize(1.3);
    graph -> SetMarkerStyle(20);

    int runIdx = mOptContainer->GetRunType();
    for(int run=0; run<kRunNum; run++){
        if(runIdx != kALLRun && runIdx != run){continue;}
        TString runName = mOptContainer -> GetRunTypeName(run);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer->GetDLEName(dle);

                cBinning -> Clear();
                cBinning -> Divide(3, 1);
                
                cBinning -> cd(1);
                TH1D* ptHist = mBinning -> GetPtHist(run, type, dle);
                ptHist -> GetYaxis()->SetRangeUser(0, 1.3*ptHist->GetBinContent(ptHist->GetMaximumBin()));
                ptHist -> Draw("hist");

                latex -> DrawLatexNDC(0.15, 0.85, GetSystemString());
                if(mOptContainer->GetParticleRunIdx() == kPi0Run){latex -> DrawLatexNDC(0.15, 0.79, Form("#pi^{0} type %i", type+1));}
                latex -> DrawLatexNDC(0.15, 0.73, Form("%s", dleName.Data()));

                cBinning -> cd(2);
                TH1D* xfHist = mBinning -> GetXfHist(run, type, dle);
                xfHist -> GetYaxis()->SetRangeUser(0, 1.3*xfHist->GetBinContent(xfHist->GetMaximumBin()));
                xfHist -> Draw("hist");
                
                latex -> DrawLatexNDC(0.15, 0.85, GetSystemString());
                if(mOptContainer->GetParticleRunIdx() == kPi0Run){latex -> DrawLatexNDC(0.15, 0.79, Form("#pi^{0} type %i", type+1));}
                latex -> DrawLatexNDC(0.15, 0.73, Form("%s", dleName.Data()));

                cBinning -> cd(3);
                TH2D* ptxfHist = mBinning -> GetKinematicsHist(run, type, dle);
                ptxfHist -> Draw("colz");

                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                for(int pt=0; pt<ptNum+1; pt++){
                    double ptBoundary = mBinning -> GetPtBinBoundary(run, type, dle, pt);
                    TLine* line = new TLine(0., ptBoundary, 1., ptBoundary);
                    line -> SetLineColor(kBlack);
                    line -> SetLineWidth(2);
                    line -> SetLineStyle(2);
                    line -> Draw("same");
                }

                for(int xf=0; xf<xfNum; xf++){
                    double xfBoundary = mBinning -> GetXfBinBoundary(run, type, dle, xf);
                    TLine* line = new TLine(xfBoundary, 0., xfBoundary, 1.);
                    line -> SetLineColor(kBlack);
                    line -> SetLineWidth(2);
                    line -> SetLineStyle(2);
                    line -> Draw("same");
                }

                graph -> Set(0);
                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        double ptMean = mBinning -> GetPtBinMean(run, type, dle, pt, xf);
                        double xfMean = mBinning -> GetXfBinMean(run, type, dle, pt, xf);
                        graph -> SetPoint(graph->GetN(), xfMean, ptMean);
                    }
                }
                graph -> Draw("same, p");

                cBinning -> Draw();
                cBinning -> SaveAs(Form("%s/Binning_%s_%s_type%i_%s.pdf", figurePath.Data(), ParticleTypeName.Data(), runName.Data(), type+1, dleName.Data()));
            }
        }
    }
}

void RHICfFigureDrawing::DrawDilutionHist()
{
    TString figurePath = mOptContainer -> GetFigurePath() + "/Dilution";
    TString ParticleTypeName = mOptContainer -> GetParticleRunName();

    TCanvas* cDilution;
    TLatex* latex = new TLatex();

    int runIdx = mOptContainer->GetRunType();
    for(int run=0; run<kRunNum; run++){
        if(runIdx != kALLRun && runIdx != run){continue;}
        TString runName = mOptContainer -> GetRunTypeName(run);
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer->GetDLEName(dle);

                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                if(ptNum <= 0){continue;}
                
                cDilution = new TCanvas(Form("cDilution_%i_%i_%i", run, type, dle), "", double(xfNum*700), double(ptNum*700));
                cDilution -> Clear();
                cDilution -> Divide(xfNum, ptNum);

                int tmpcIdx = 0;
                for(int pt=0; pt<ptNum; pt++){
                    for(int xf=0; xf<xfNum; xf++){
                        const int cIdx = tmpcIdx+1;
                        cDilution -> cd(cIdx);

                        TH1D* dilution = mDilution -> GetDilutionHist(run, type, dle, pt, xf);
                        dilution -> GetYaxis()->SetRangeUser(0, 1.5*dilution->GetBinContent(dilution->GetMaximumBin()));
                        dilution -> Draw("hist");

                        latex -> DrawLatexNDC(0.15, 0.85, GetSystemString());
                        if(mOptContainer->GetParticleRunIdx() == kPi0Run){latex -> DrawLatexNDC(0.15, 0.79, Form("#pi^{0} type %i", type+1));}
                        latex -> DrawLatexNDC(0.15, 0.73, Form("%s", dleName.Data()));

                        latex -> DrawLatexNDC(0.5, 0.79, Form("%.2f < p_{T} < %.2f", mBinning->GetPtBinBoundary(run, type, dle, pt), mBinning->GetPtBinBoundary(run, type, dle, pt+1) ));
                        latex -> DrawLatexNDC(0.5, 0.73, Form("%.2f < x_{F} < %.2f", mBinning->GetXfBinBoundary(run, type, dle, xf), mBinning->GetXfBinBoundary(run, type, dle, xf+1) ));
                        latex -> DrawLatexNDC(0.5, 0.67, Form("D_{#phi} = %.3f", mDilution->GetDilutionFactor(run, type, dle, pt, xf)));

                        tmpcIdx++;
                    }
                }
                cDilution -> Draw();
                cDilution -> SaveAs(Form("%s/Dilution_%s_%s_type%i_%s.pdf", figurePath.Data(), ParticleTypeName.Data(), runName.Data(), type+1, dleName.Data()));
            }
        }
    }
}

void RHICfFigureDrawing::DrawAsymmetryGraph()
{
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
    if(particleTypeIdx == kPi0Run){
        legendPosY[0] = 0.65;
        legendPosY[1] = 0.89;
        drawANRange[0] = -0.3;
        drawANRange[1] = 0.3;
        latexPosY[0] = 0.13;
        latexPosY[1] = 0.06;
        anTypeNum = 3;
    }
    if(particleTypeIdx == kNeutronRun){
        legendPosY[0] = 0.12;
        legendPosY[1] = 0.4;
        drawANRange[0] = -0.35;
        drawANRange[1] = 0.1;
        latexPosY[0] = 0.85;
        latexPosY[1] = -0.06;
        anTypeNum = 1;
    }

    TCanvas* cAsymmetry = new TCanvas("cAsymmetry", "", 600.*2., 600.*3.);
    TLatex* latex = new TLatex();
    TLegend* legend = new TLegend(0.12, legendPosY[0], 0.4, legendPosY[1]);
    legend -> SetBorderSize(0);

    TLine* line = new TLine(0., 0., 1., 0.);
    line -> SetLineStyle(2);
    line -> SetLineColor(kBlack);

    const int ptNum = mBinning -> GetGlobalPtBinNum();
    const int xfNum = mBinning -> GetGlobalXfBinNum();

    TGraph* base = new TGraph();
    base -> SetPoint(0, 0., -1.);
    base -> SetPoint(1, 1.5, 1.);
    base -> SetMarkerSize(0);
    base -> SetMarkerStyle(20);
    base -> SetLineWidth(0);
    base -> GetYaxis() -> SetRangeUser(drawANRange[0], drawANRange[1]); 
    base -> GetXaxis() -> SetRangeUser(0., 1.);

    for(int dle=0; dle<kDLENum; dle++){
        TString dleName = mOptContainer->GetDLEName(dle);
        if(dle == kDLENum-1){dleName = "Inclusive";}

        for(int ptxf=0; ptxf<2; ptxf++){
            bool isPt = (ptxf==0)? true : false;
            int binNum = (isPt)? xfNum : ptNum;

            for(int anType=0; anType<anTypeNum; anType++){
                TString anTypeName = "";
                if(anType == 0){anTypeName = "Inclusive";}
                if(anType == 1){anTypeName = "Bkg";}
                if(anType == 2){anTypeName = "Exclusive";}
            
                cAsymmetry -> Clear();
                cAsymmetry -> Divide(2, 3);

                int tmpcIdx = 0;
                for(int bin=0; bin<binNum; bin++){
                    const int cIdx = tmpcIdx+1;
                    cAsymmetry -> cd(cIdx);
                    TString xTitle = (isPt)? "p_{T} [GeV/c]" : "x_{F}";
                    base -> SetTitle(Form("; %s; A_{N}", xTitle.Data()));
                    base -> Draw("ap");

                    legend -> Clear("ICESM");

                    for(int run=0; run<kRunNum; run++){
                        TString runName = mOptContainer -> GetRunTypeName(run);
                        for(int type=0; type<kTypeNum; type++){
                        
                            TGraphErrors* graph = mAsymmetry -> GetANGraph(isPt, anType, run, type, dle, bin);
                            graph -> SetMarkerColor(runColor[run]);
                            graph -> SetMarkerStyle(typeStyle[type]);
                            graph -> SetMarkerSize(1.2);
                            graph -> SetLineWidth(1.5);
                            graph -> SetLineColor(runColor[run]);
                            graph -> Draw("p, same");

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
                        latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning->GetGlobalXfBinBoundary(bin), binName.Data(), mBinning->GetGlobalXfBinBoundary(bin+1) ));
                    }
                    else{
                        binName = "p_{T}";
                        latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning->GetGlobalPtBinBoundary(bin), binName.Data(), mBinning->GetGlobalPtBinBoundary(bin+1) ));
                    }

                    tmpcIdx++;
                }
                cAsymmetry -> Draw();
                cAsymmetry -> SaveAs(Form("%s/AN%s_%s_%s_kinematics_%i.pdf", figurePath.Data(), anTypeName.Data(), ParticleTypeName.Data(), dleName.Data(), ptxf));
            }
        }
    }

    delete cAsymmetry;
    cAsymmetry = new TCanvas("cAsymmetry_all", "", 600.*2., 600.*2.);

    for(int ptxf=0; ptxf<2; ptxf++){
        bool isPt = (ptxf==0)? true : false;
        int binNum = (isPt)? xfNum : ptNum;

        for(int anType=0; anType<anTypeNum; anType++){
            TString anTypeName = "";
            if(anType == 0){anTypeName = "Inclusive";}
            if(anType == 1){anTypeName = "Bkg";}
            if(anType == 2){anTypeName = "Exclusive";}

            cAsymmetry -> Clear();
            cAsymmetry -> Divide(2, 2);
            
            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer->GetDLEName(dle);
                if(dle == kDLENum-1){dleName = "Inclusive";}

                const int cIdx = dle+1;
                cAsymmetry -> cd(cIdx);
                TString xTitle = (isPt)? "p_{T} [GeV/c]" : "x_{F}";
                base -> SetTitle(Form("; %s; A_{N}", xTitle.Data()));
                base -> Draw("ap");

                legend -> Clear("ICESM");

                for(int bin=0; bin<binNum; bin++){
                    TGraphErrors* graph = mAsymmetry -> GetANGraph(isPt, anType, dle, bin);
                    graph -> SetMarkerColor(binColor[bin]);
                    graph -> SetMarkerStyle(typeStyle[0]);
                    graph -> SetMarkerSize(1.2);
                    graph -> SetLineWidth(1.5);
                    graph -> SetLineColor(binColor[bin]);
                    graph -> Draw("p, same");

                    TString binName = "";
                    if(isPt){
                        binName = "x_{F}";
                        legend -> AddEntry(graph, Form("%.2f < %s < %.2f", mBinning->GetGlobalXfBinBoundary(bin), binName.Data(), mBinning->GetGlobalXfBinBoundary(bin+1)), "p");
                    }
                    else{
                        binName = "p_{T}";
                        legend -> AddEntry(graph, Form("%.2f < %s < %.2f", mBinning->GetGlobalPtBinBoundary(bin), binName.Data(), mBinning->GetGlobalPtBinBoundary(bin+1)), "p");
                    }
                }

                legend -> Draw("same");
                line -> Draw("same");
                latex -> DrawLatexNDC(0.13, latexPosY[0], GetSystemString());
                latex -> DrawLatexNDC(0.13, latexPosY[0]+latexPosY[1], Form("%s", dleName.Data()));
                latex -> DrawLatexNDC(0.33, latexPosY[0]+latexPosY[1], "#eta > 6");
            }
            cAsymmetry -> Draw();
            cAsymmetry -> SaveAs(Form("%s/AN%s_%s_%i.pdf", figurePath.Data(), anTypeName.Data(), ParticleTypeName.Data(), ptxf));
        }
    }

    delete cAsymmetry;
    cAsymmetry = new TCanvas("cAsymmetry_all", "", 600.*2., 600.*3.);

    for(int ptxf=0; ptxf<2; ptxf++){
        bool isPt = (ptxf==0)? true : false;
        int binNum = (isPt)? xfNum : ptNum;

        for(int anType=0; anType<anTypeNum; anType++){
            TString anTypeName = "";
            if(anType == 0){anTypeName = "Inclusive";}
            if(anType == 1){anTypeName = "Bkg";}
            if(anType == 2){anTypeName = "Exclusive";}

            cAsymmetry -> Clear();
            cAsymmetry -> Divide(2, 3);
            
            int tmpcIdx = 0;
            for(int bin=0; bin<binNum; bin++){
                const int cIdx = tmpcIdx+1;
                cAsymmetry -> cd(cIdx);
                TString xTitle = (isPt)? "p_{T} [GeV/c]" : "x_{F}";
                base -> SetTitle(Form("; %s; A_{N}", xTitle.Data()));
                base -> Draw("ap");

                legend -> Clear("ICESM");

                for(int dle=0; dle<kDLENum; dle++){
                    TString dleName = mOptContainer->GetDLEName(dle);
                    if(dle == kDLENum-1){dleName = "Inclusive";}

                    TGraphErrors* graph = mAsymmetry -> GetANGraph(isPt, anType, dle, bin);
                    graph -> SetMarkerColor(markerColorDLE[dle]);
                    graph -> SetMarkerStyle(markerStyleDLE[dle]);
                    graph -> SetMarkerSize(1.2);
                    graph -> SetLineWidth(1.5);
                    graph -> SetLineColor(markerColorDLE[dle]);
                    graph -> Draw("p, same");
                    legend -> AddEntry(graph, Form("%s", dleName.Data()), "p");
                }

                legend -> Draw("same");
                line -> Draw("same");
                latex -> DrawLatexNDC(0.13, latexPosY[0], GetSystemString());
                latex -> DrawLatexNDC(0.13, latexPosY[0]+latexPosY[1], "#eta > 6");
                TString binName = "";
                if(isPt){
                    binName = "x_{F}";
                    latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning->GetGlobalXfBinBoundary(bin), binName.Data(), mBinning->GetGlobalXfBinBoundary(bin+1) ));
                }
                else{
                    binName = "p_{T}";
                    latex -> DrawLatexNDC(0.5, latexPosY[0]+latexPosY[1], Form("%.2f < %s < %.2f", mBinning->GetGlobalPtBinBoundary(bin), binName.Data(), mBinning->GetGlobalPtBinBoundary(bin+1) ));
                }
                tmpcIdx++;
            }
            cAsymmetry -> Draw();
            cAsymmetry -> SaveAs(Form("%s/AN%s_%s_%i_DLECompare.pdf", figurePath.Data(), anTypeName.Data(), ParticleTypeName.Data(), ptxf));
        }
    }



    delete cAsymmetry;
    cAsymmetry = new TCanvas("cAsymmetry_Summary", "", 600.*2., 600.);

    for(int anType=0; anType<anTypeNum; anType++){
        TString anTypeName = "";
        if(anType == 0){anTypeName = "Inclusive";}
        if(anType == 1){anTypeName = "Bkg";}
        if(anType == 2){anTypeName = "Exclusive";}

        cAsymmetry -> Clear();
        cAsymmetry -> Divide(2, 1);

        for(int ptxf=0; ptxf<2; ptxf++){
            bool isPt = (ptxf==0)? true : false;
            int binNum = (isPt)? xfNum : ptNum;

            const int cIdx = ptxf+1;
            cAsymmetry -> cd(cIdx);
            TString xTitle = (isPt)? "p_{T} [GeV/c]" : "x_{F}";
            base -> SetTitle(Form("; %s; A_{N}", xTitle.Data()));
            base -> GetYaxis() -> SetRangeUser(drawANRange[0], drawANRange[1]+0.1); 
            base -> Draw("ap");

            legend -> Clear("ICESM");

            for(int dle=0; dle<kDLENum; dle++){
                TString dleName = mOptContainer->GetDLEName(dle);
                if(dle == kDLENum-1){dleName = "Inclusive";}

                TGraphErrors* graph = mAsymmetry -> GetANSummaryGraph(isPt, anType, dle);
                graph -> SetMarkerStyle(markerStyleDLE[dle]);
                graph -> SetMarkerColor(markerColorDLE[dle]);
                graph -> SetLineColor(markerColorDLE[dle]);
                graph -> SetMarkerSize(1.2);
                graph -> SetLineWidth(1.5);
                graph -> Draw("same, p");

                legend -> AddEntry(graph, Form("%s", dleName.Data()), "p");
            }

            legend -> Draw("same");
            line -> Draw("same");
            latex -> DrawLatexNDC(0.13, latexPosY[0], GetSystemString());
            latex -> DrawLatexNDC(0.13, latexPosY[0]+latexPosY[1], "#eta > 6");
        }

        cAsymmetry -> Draw();
        cAsymmetry -> SaveAs(Form("%s/ANSummary_%s_%s.pdf", figurePath.Data(), anTypeName.Data(), ParticleTypeName.Data()));
    }
}

TString RHICfFigureDrawing::GetSystemString()
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
