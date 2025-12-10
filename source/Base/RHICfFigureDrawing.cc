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
// void RHICfFigureDrawing::SetAsymmetry(RHICfAsymmetry* asymmetry){mAsymmetry = asymmetry;}

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
                latex -> DrawLatexNDC(0.53, 0.53, Form("#chi^{2}/NDF = %.2f", mFitting->GetMassFitChi2(run, type, dle)));
                latex -> DrawLatexNDC(0.53, 0.47, Form("B/S ratio = %.3f", mFitting->GetBSRatio(run, type, dle)));
                latex -> DrawLatexNDC(0.53, 0.41, Form("          #pm %.3f", mFitting->GetBSRatioErr(run, type, dle)));

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
                        latex -> DrawLatexNDC(0.53, 0.53, Form("#chi^{2}/NDF = %.2f", mFitting->GetMassFitChi2(run, type, dle, pt, xf)));
                        latex -> DrawLatexNDC(0.53, 0.47, Form("B/S ratio = %.3f", mFitting->GetBSRatio(run, type, dle, pt, xf)));
                        latex -> DrawLatexNDC(0.53, 0.41, Form("          #pm %.3f", mFitting->GetBSRatioErr(run, type, dle, pt, xf)));
                        
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
                ptxfHist -> SetTitle("; x_{F}; p_{T} [GeV/c]");
                ptxfHist -> Draw("colz");
                latex -> DrawLatexNDC(0.15, 0.85, GetSystemString());
                if(mOptContainer->GetParticleRunIdx() == kPi0Run){latex -> DrawLatexNDC(0.15, 0.79, Form("%s run, #pi^{0} type %i", runName.Data(), type+1));}
                latex -> DrawLatexNDC(0.15, 0.73, Form("%s", dleName.Data()));

                int ptNum = mBinning -> GetPtBinNum(run, type, dle);
                int xfNum = mBinning -> GetXfBinNum(run, type, dle);
                // for(int pt=0; pt<ptNum+1; pt++){
                //     double ptBoundary = mBinning -> GetPtBinBoundary(run, type, dle, pt);
                //     TLine* line = new TLine(0., ptBoundary, 1., ptBoundary);
                //     line -> SetLineColor(kBlack);
                //     line -> SetLineWidth(2);
                //     line -> SetLineStyle(2);
                //     line -> Draw("same");
                // }

                // for(int xf=0; xf<xfNum; xf++){
                //     double xfBoundary = mBinning -> GetXfBinBoundary(run, type, dle, xf);
                //     TLine* line = new TLine(xfBoundary, 0., xfBoundary, 1.);
                //     line -> SetLineColor(kBlack);
                //     line -> SetLineWidth(2);
                //     line -> SetLineStyle(2);
                //     line -> Draw("same");
                // }

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
                        latex -> DrawLatexNDC(0.65, 0.67, Form("D_{#phi} = %.3f", mDilution->GetDilutionFactor(run, type, dle, pt, xf)));
                        latex -> DrawLatexNDC(0.65, 0.61, Form("   #pm %.3f", mDilution->GetDilutionFactorErr(run, type, dle, pt, xf)));
                        tmpcIdx++;
                    }
                }
                cDilution -> Draw();
                cDilution -> SaveAs(Form("%s/Dilution_%s_%s_type%i_%s.pdf", figurePath.Data(), ParticleTypeName.Data(), runName.Data(), type+1, dleName.Data()));
            }
        }
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
