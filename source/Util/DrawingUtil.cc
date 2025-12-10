#include "DrawingUtil.hh"

DrawingUtil::DrawingUtil(TString figurePath)
{
    mFigurePath = figurePath;
    if(mFigurePath[mFigurePath.Sizeof()-2] != '/'){mFigurePath = mFigurePath + "/";}

    mColorPaletteIdx = 0;
    mMarkerStyleIdx = 0;
}

DrawingUtil::~DrawingUtil()
{
}

void DrawingUtil::Init()
{

}

void DrawingUtil::Clear(TString opt)
{
    if(mCanvas != nullptr){delete mCanvas;}
    mDrawOptArr.clear();
    mLabelName.clear();
    mHist.clear();
    mGraph.clear();

    bool isClearLegend = false;
    bool isClearText = false;
    opt.ToUpper();
    if(opt.Index("ALL") != -1){
        isClearLegend = true;
        isClearText = true;
    }
    if(opt.Index("TEXT") != -1 || opt.Index("LATEX") != -1 || opt.Index("TXT") != -1){isClearText = true;}
    if(opt.Index("LEG") != -1){isClearLegend = true;}

    int latexNum = mLatex.size();
    for(int i=0; i<latexNum; i++){
        if(mLatex[i].persistency && !isClearText){continue;}
        mLatex[i].Clear();
    }
    int legendNum = mLegend.size();
    for(int i=0; i<legendNum; i++){
        if(mLegend[i].persistency && !isClearLegend){continue;}
        mLegend[i].Clear();
    }
}

void DrawingUtil::DrawHist(TString opt)
{
    InitDrawOption(opt);

    int cNum = 0;
    InitHist(cNum);

    int histNum = mHist.size();
    for(int c=0; c<cNum; c++){
        mCanvas -> cd(c+1);

        if(GetDrawFlag("xLog")){gPad -> SetLogx();}
        if(GetDrawFlag("yLog")){gPad -> SetLogy();}

        gPad -> SetTopMargin(0.03);
        gPad -> SetLeftMargin(0.135);
        gPad -> SetBottomMargin(0.11);

        double maxBin;
        for(int i=0; i<histNum; i++){
            if(mHist[i].base.cIdx == c && mHist[i].base.objIdx == -1){
                mHist[i].hist -> GetYaxis()->SetTitleSize(0.045);
                mHist[i].hist -> GetXaxis()->SetTitleSize(0.045);
                mHist[i].hist -> GetYaxis()->SetTitleOffset(1.6);
                mHist[i].hist -> GetXaxis()->SetTitleOffset(1.);
                mHist[i].hist -> Draw();

                int nBinsX = mHist[i].hist ->GetNbinsX();
                maxBin = mHist[i].hist -> GetBinLowEdge(nBinsX) + mHist[i].hist->GetBinWidth(nBinsX);
            }
        }
        for(int i=0; i<histNum; i++){
            if(mHist[i].base.cIdx == c && mHist[i].base.objIdx != -1){
                mHist[i].hist -> Draw("same");
            }
        }
        DrawText(c);
        DrawLegend(c);

        if(GetDrawFlag("ratio")){
            TLine* ratioLine = new TLine(0., 1., maxBin, 1.);
            ratioLine -> SetLineColor(kBlack);
            ratioLine -> SetLineStyle(2);
            ratioLine -> SetLineWidth(1);
            ratioLine -> Draw("same");
        }
    }
}

void DrawingUtil::DrawHistWithRatio(int baseIdx, TString opt, TString yTitle)
{
    InitDrawOption(opt);
    int cNum = 0;
    InitHist(cNum);

    int histNum = mHist.size();
    for(int c=0; c<cNum; c++){
        mCanvas -> cd(c+1);

        TH1D* bkgRatioHist;
        TH1D* baseRatioHist;
        vector<TH1D*> ratioHistArr;

        TString canvasName = mCanvas -> GetName();
        TPad* mainPad = new TPad(Form("%s_mainPad_c%i", canvasName.Data(), c), "", 0.1, 0.3, 0.9, 0.9);  
        mainPad -> Draw();
        mainPad -> SetTopMargin(0.03);
        mainPad -> SetLeftMargin(0.135);
        mainPad -> SetBottomMargin(0.2);
        mainPad -> SetGrid(0,0);

        TPad* ratioPad = new TPad(Form("%s_ratioPad_c%i", canvasName.Data(), c), "", 0.1, 0.1, 0.9, 0.415+0.005);  
        ratioPad -> Draw();
        ratioPad -> SetTopMargin(0);
        ratioPad -> SetLeftMargin(0.135);
        ratioPad -> SetBottomMargin(0.3);
        ratioPad -> SetGrid(0,0);

        mainPad -> cd();
        if(GetDrawFlag("xLog")){mainPad -> SetLogx();}
        if(GetDrawFlag("yLog")){gPad -> SetLogy();}

        for(int i=0; i<histNum; i++){
            if(mHist[i].base.cIdx == c && mHist[i].base.objIdx == -1){
                mHist[i].hist -> GetYaxis()->SetTitleSize(0.045);
                mHist[i].hist -> GetXaxis()->SetTitleSize(0);
                mHist[i].hist -> GetYaxis()->SetTitleOffset(1.6);
                mHist[i].hist -> GetXaxis()->SetTitleOffset(0.);
                mHist[i].hist -> Draw();

                bkgRatioHist = (TH1D*)mHist[i].hist -> Clone(Form("%s_bkgRatioHist", mHist[i].hist->GetName()));
            }
        }
        for(int i=0; i<histNum; i++){
            if(mHist[i].base.cIdx == c && mHist[i].base.objIdx != -1){
                mHist[i].hist -> Draw("same");

                if(mHist[i].base.objIdx == baseIdx){
                    mHist[i].hist -> SetLineWidth(2);
                    baseRatioHist = (TH1D*)mHist[i].hist -> Clone(Form("%s_RatioBaseHist_c%i", mHist[i].hist, c));
                }
                else{
                    mHist[i].hist -> SetLineWidth(1);
                    TH1D* ratioHist = (TH1D*)mHist[i].hist -> Clone(Form("%s_RatioHist_c%i_%i", mHist[i].hist, c, i));
                    ratioHistArr.push_back(ratioHist);
                }
            }
        }

        DrawText(c);
        DrawLegend(c);

        ratioPad -> cd();
        bkgRatioHist -> GetYaxis()->SetTitleOffset(0.95);
        bkgRatioHist -> GetXaxis()->SetTitleOffset(1.2);
        bkgRatioHist -> GetYaxis() -> SetTitleSize(0.08);
        bkgRatioHist -> GetXaxis() -> SetTitleSize(0.08);
        bkgRatioHist -> GetYaxis() -> SetLabelSize(0.06);
        bkgRatioHist -> GetXaxis() -> SetLabelSize(0.07);
        bkgRatioHist -> GetXaxis() -> SetLabelOffset(0.02);
        bkgRatioHist -> GetYaxis()->SetRangeUser(0.001, 1.999);
               
        TString xTitleTmp = bkgRatioHist -> GetXaxis()->GetTitle();
        bkgRatioHist -> SetTitle(Form(";%s;%s", xTitleTmp.Data(), yTitle.Data()));
        bkgRatioHist -> Draw();

        int ratioHistNum = ratioHistArr.size();
        int binNum = baseRatioHist -> GetNbinsX();

        for(int i=0; i<ratioHistNum; i++){
            for(int b=0; b<binNum; b++){
                double entry = ratioHistArr[i] -> GetBinContent(b+1);
                double baseEntry = baseRatioHist -> GetBinContent(b+1);

                double ratio = entry / baseEntry;
                if(baseEntry < 1.e-4 || entry < 1.e-4){ratio = 0.;}
                if(ratio >= 1.999){ratio = 2.;}
                ratioHistArr[i] -> SetBinContent(b+1, ratio);
            }
            ratioHistArr[i] -> SetMarkerColor(ratioHistArr[i]->GetLineColor());
            ratioHistArr[i] -> Draw("same");
        }

        double nBinsX = bkgRatioHist -> GetNbinsX();
        double binLower = bkgRatioHist -> GetBinLowEdge(1);
        double binHigher = bkgRatioHist -> GetBinLowEdge(nBinsX) + bkgRatioHist->GetBinWidth(nBinsX);

        TLine* ratioLine = new TLine(binLower, 1., binHigher, 1.);
        ratioLine -> SetLineColor(kBlack);
        ratioLine -> SetLineStyle(2);
        ratioLine -> SetLineWidth(1);
        ratioLine -> Draw("same");
    }
}

void DrawingUtil::DrawGraph(TString opt)
{
    InitDrawOption(opt);
    int cNum = 0;
    InitGraph(cNum);

    int graphNum = mGraph.size();
    for(int c=0; c<cNum; c++){
        mCanvas -> cd(c+1);
        gPad -> SetTopMargin(0.03);
        gPad -> SetLeftMargin(0.135);
        gPad -> SetBottomMargin(0.11);

        double minX, maxX, minY, maxY;
        for(int i=0; i<graphNum; i++){
            if(mGraph[i].base.cIdx == c && mGraph[i].base.objIdx == -1){
                mGraph[i].graph -> GetYaxis()->SetTitleSize(0.045);
                mGraph[i].graph -> GetXaxis()->SetTitleSize(0.045);
                mGraph[i].graph -> GetYaxis()->SetTitleOffset(1.6);
                mGraph[i].graph -> GetXaxis()->SetTitleOffset(1.);
                mGraph[i].graph -> Draw("ap");

                mGraph[i].graph->GetPoint(0, minX, minY);
                mGraph[i].graph->GetPoint(mGraph[i].graph->GetN()-1, maxX, maxY);
            }
            for(int i=0; i<graphNum; i++){
                if(mGraph[i].base.cIdx == c && mGraph[i].base.objIdx != -1){
                    mGraph[i].graph -> Draw("same, p");
                }
            }
        }
        DrawText(c);
        DrawLegend(c);

        if(GetDrawFlag("Ratio")){
            TLine* ratioLine = new TLine(minX, 1., maxX, 1.);
            ratioLine -> SetLineColor(kBlack);
            ratioLine -> SetLineStyle(2);
            ratioLine -> SetLineWidth(1);
            ratioLine -> Draw("same");
        }
    }
}


void DrawingUtil::SaveFigure(TString name)
{
    TString figureFormat = ".pdf";

    TString saveName = mFigurePath;
    if(name.Index("/") != -1){saveName = name;}
    else{saveName = saveName + name;}

    TObjArray *tokens = saveName.Tokenize(".");
    int entries = tokens -> GetEntries();
    if(entries == 0){
        saveName = saveName + figureFormat;
    }
    else{
        TString tmpFormat = ((TObjString *)tokens->At(entries-1))->GetString();
        tmpFormat.ToUpper();
        if(tmpFormat.Index("PDF") == -1 && tmpFormat.Index("PNG") == -1 && tmpFormat.Index("PDF") == -1){
            saveName = saveName + figureFormat;
        }
    }

    mCanvas -> Draw();
    mCanvas -> SaveAs(Form("%s", saveName.Data()));
}

void DrawingUtil::SetCanvas(int divideNum, TString name, double size)
{   
    int xNum = 1;
    int yNum = 1;
    if(divideNum == 4){xNum = 2;}
    else if(divideNum < 4){xNum = divideNum;}
    else if(divideNum < 10){xNum = 3;}
    else if(divideNum > 10){xNum = 4;}
    double tmp = double(divideNum)/double(xNum);
    yNum = int(ceil(tmp));
    if(name == ""){name = "canvas1";}
    mCanvas = new TCanvas(Form("%s_%i_%i", name.Data(), xNum, yNum), "", double(xNum)*size, double(yNum)*size);
    mCanvas -> Divide(xNum, yNum);
}

void DrawingUtil::SetXLabelName(int i, TString name)
{
    mLabelName.push_back(make_pair(i, name));
}

void DrawingUtil::SetTH1DRatio(int num, double* y, double base, TString title, int cIdx, int colorIdx, int markerIdx)
{
    InitDrawOption("ratio");
    TH1D* hist = new TH1D(Form("th1dRatio_%s_c%i_%.2f", title.Data(), cIdx, base), "", num, 0., double(num));
    hist -> Sumw2();

    for(int b=0; b<num; b++){
        double ratio = y[b]/base;
        double error = sqrt(y[b])/base;
        if(y[b] < 1. || base < 1.){
            ratio = 0.;
            error = 0.;
        }
        hist -> SetBinContent(b+1, ratio);
        hist -> SetBinError(b+1, error);
    }

    SetTH1D(hist, title, cIdx, colorIdx, markerIdx);
}

void DrawingUtil::SetTH1D(TH1D* hist, TString title, int cIdx, int colorIdx, int markerIdx)
{
    int histNum = mHist.size();
    int objNum = 0;
    for(int i=0; i<histNum; i++){
        if(mHist[i].base.cIdx == cIdx){
            objNum++;
        }
    }

    InfoHist histInfo;

    TH1D* clone = (TH1D*)hist -> Clone(Form("%s_clone_th1d_c%i_obj%i", TString(hist->GetName()).Data(), cIdx, objNum));
    histInfo.hist = clone;

    int color = (colorIdx == -1)? GetColor(objNum) : colorIdx;
    hist -> SetLineColor(color);
    hist -> SetMarkerColor(color);
    hist -> SetLineWidth(2.5);
    histInfo.hist -> SetLineColor(color);
    histInfo.hist -> SetMarkerColor(color);

    if(markerIdx != -1){
        hist -> SetMarkerStyle(GetMarker(markerIdx));
        histInfo.hist -> SetMarkerStyle(GetMarker(markerIdx));
    }
    histInfo.base.Clear();
    histInfo.base.title = title;
    histInfo.base.cIdx = cIdx;
    histInfo.base.objIdx = objNum;
    histInfo.base.colorIdx = colorIdx;
    histInfo.base.markerIdx = markerIdx;

    mHist.push_back(histInfo);
}

TH1D* DrawingUtil::GetTH1D(int cIdx, int objIdx)
{
    int histNum = mHist.size();
     for(int i=0; i<histNum; i++){
        if(mHist[i].base.cIdx == cIdx){
            if(mHist[i].base.objIdx == objIdx){
                return mHist[i].hist;
            }
        }
    }
    return 0;
}

void DrawingUtil::SetGraph(int num, double* x, double* y, TString title, int cIdx, int colorIdx, int markerIdx)
{
    TGraphErrors* graph = SetGraph(title, cIdx, colorIdx, markerIdx);
    for(int i=0; i<num; i++){
        graph -> SetPoint(graph->GetN(), x[i], y[i]);
    }
}

void DrawingUtil::SetGraphRatio(int num, double* x, double* y, double* yBase, TString title, int cIdx, int colorIdx, int markerIdx)
{
    TGraphErrors* graph = SetGraph(title, cIdx, colorIdx, markerIdx);
    for(int i=0; i<num; i++){
        double ratio = y[i]/yBase[i];
        double ratioError = sqrt(y[i])/yBase[i];
        if(yBase[i] < 1e-5 || y[i] < 1e-5){ratio = 0.; ratioError= 0.;}
        graph -> SetPoint(graph->GetN(), x[i], ratio);
        graph -> SetPointError(graph->GetN()-1, 0., ratioError);
    }
}

void DrawingUtil::SetGraphError(int num, double* x, double* y, double* ex, double* ey, TString title, int cIdx, int colorIdx, int markerIdx)
{
    TGraphErrors* graph = SetGraph(title, cIdx, colorIdx, markerIdx);
    for(int i=0; i<num; i++){
        graph -> SetPoint(graph->GetN(), x[i], y[i]);
        graph -> SetPointError(graph->GetN()-1, ex[i], ey[i]);
    }
}

TGraphErrors* DrawingUtil::SetGraph(TString title, int cIdx, int colorIdx, int markerIdx)
{
    int graphNum = mGraph.size();
    int objNum = 0;
    for(int i=0; i<graphNum; i++){
        if(mGraph[i].base.cIdx == cIdx){
            objNum++;
        }
    }

    TGraphErrors* graph = new TGraphErrors();
    graph -> SetName(Form("tmpGraph_c%i_obj%i", cIdx, objNum));

    int color = (colorIdx == -1)? GetColor(objNum) : colorIdx;
    graph -> SetLineColor(color);
    graph -> SetMarkerColor(color);
    graph -> SetLineWidth(1.5);
    graph -> SetMarkerSize(2);
    graph -> SetLineColor(color);
    graph -> SetMarkerColor(color);

    int marker = (markerIdx == -1)? GetMarker(objNum) : markerIdx;
    graph -> SetMarkerStyle(marker);

    InfoGraph graphInfo;
    graphInfo.graph = graph;
    graphInfo.base.Clear();
    graphInfo.base.title = title;
    graphInfo.base.cIdx = cIdx;
    graphInfo.base.objIdx = objNum;
    graphInfo.base.colorIdx = colorIdx;
    graphInfo.base.markerIdx = markerIdx;

    mGraph.push_back(graphInfo);
    return graph;
}

TGraphErrors* DrawingUtil::GetGraph(int cIdx, int objIdx)
{
    int graphNum = mGraph.size();
     for(int i=0; i<graphNum; i++){
        if(mGraph[i].base.cIdx == cIdx){
            if(mGraph[i].base.objIdx == objIdx){
                return mGraph[i].graph;
            }
        }
    }
    return 0; 
}

void DrawingUtil::SetText(bool isPersistant, TString text, int cIdx, double x, double y, double size, int font)
{
    InfoText textInfo;
    textInfo.Clear();
    textInfo.persistency = isPersistant;
    textInfo.text = text;
    textInfo.cIdx = cIdx;
    textInfo.x = x;
    textInfo.y = y;
    textInfo.size = size;
    textInfo.font = font;
    mLatex.push_back(textInfo);
}

void DrawingUtil::SetLegend(bool isPersistant, TObject* obj, TString text, int cIdx, TString opt)
{   
    InfoText textInfo;
    textInfo.Clear();
    textInfo.persistency = isPersistant;
    textInfo.text = text;
    textInfo.cIdx = cIdx;

    if(opt == ""){
        TString name = obj->GetName();
        name.ToUpper();
        if(name.Index("GRAPH") != -1){
            opt = "p";
        }
        else{
            opt = "l";
        }
    }

    textInfo.opt = opt;
    textInfo.obj = obj;
    mLegend.push_back(textInfo);
}

int DrawingUtil::GetColor(int idx)
{
    int Palette1[10] = {2, 4, 8, 95, 51, 91, 66, 28, 14, 1};
    int PaletteRainbow[7] = {2, 94, 90, 8, 4, kBlue+2, 51};

    if(idx < 0){return -1;}
    if(mColorPaletteIdx == 0){
        if(idx < 10){
            return Palette1[idx];
        }
        if(idx >= 10){
            return Palette1[9];
        }
    }
    else if(mColorPaletteIdx >= 1){
        if(idx < 7){
            return PaletteRainbow[idx];
        }
        if(idx >= 7){
            return 1;
        }
    }
    return -1;
}
int DrawingUtil::GetMarker(int idx)
{
    int markerStyle1[10] = {20, 21, 22, 23, 29, 33, 34, 47, 49, 39};
    int markerStyle2[10] = {24, 25, 26, 32, 30, 27, 28, 46, 36, 38};

    if(idx < 0){return -1;}
    if(mColorPaletteIdx == 0){
        if(idx < 10){
            return markerStyle1[idx];
        }
        if(idx >= 10){
            return markerStyle2[idx-10];
        }
    }
    else if(mColorPaletteIdx >= 1){
        if(idx < 10){
            return markerStyle2[idx];
        }
        if(idx >= 10){
            return markerStyle1[idx-10];
        }
    }
    return -1;
}

void DrawingUtil::InitHist(int& cNum)
{
    cNum = -1;
    int histNum = mHist.size();
    for(int i=0; i<histNum; i++){
        if(cNum < mHist[i].base.cIdx+1){
            cNum = mHist[i].base.cIdx+1;
        }
    }
    if(cNum <= 0){cNum = 1;}

    for(int c=0; c<cNum; c++){
        InfoHist histInfo;

        int bins = 0;
        double minBin = 99999.;
        double maxBin = -99999.;
        double maxHeight = 0;
        TString histName;
        TString titleX;
        TString titleY;

        for(int i=0; i<histNum; i++){
            if(mHist[i].base.cIdx == c){
                if(GetDrawFlag("norm")){
                    NormalizedTH1D(mHist[i].hist);
                }

                double nBinsX = mHist[i].hist -> GetNbinsX();
                double binLower = mHist[i].hist -> GetBinLowEdge(1);
                double binHigher = mHist[i].hist -> GetBinLowEdge(nBinsX) + mHist[i].hist->GetBinWidth(nBinsX);
                double maximum = mHist[i].hist -> GetMaximum();

                if(bins < nBinsX){bins = nBinsX;}
                if(minBin > binLower){minBin = binLower;}
                if(maxBin < binHigher){maxBin = binHigher;}
                if(maxHeight < maximum){
                    histName = mHist[i].hist->GetName();
                    maxHeight = maximum;
                }
                TString tmpName, tmpXTitle, tmpYTitle;
                bool validTitle = GetTitle(mHist[i].base.title, tmpName, tmpXTitle, tmpYTitle);

                if(validTitle){
                    titleX = tmpXTitle;
                    titleY = tmpYTitle;
                }
            }
        }

        histInfo.hist = new TH1D(Form("base_%s_TH1D_c%i", histName.Data(), c), "", bins, minBin, maxBin);
        histInfo.hist -> SetStats(0);

        histInfo.hist -> GetYaxis()->SetRangeUser(0.00011, maxHeight*1.25);
        if(GetDrawFlag("yLog")){histInfo.hist->GetYaxis()->SetRangeUser(0.00011, maxHeight*12.5);}
        if(GetDrawFlag("ratio")){
            titleY = "Ratio";
            histInfo.hist->GetYaxis()->SetRangeUser(0., 1.5);
        }

        if(GetDrawFlag("norm")){
            TObjArray *tokens = titleX.Tokenize(" ");
            TString tmp = ((TObjString *)tokens->At(0))->GetString();
            titleY = Form("1/N dN/d%s", tmp.Data());
        }
        else{
            if(titleY==""){
                titleY = "Counts";
            }
        }

        histInfo.hist -> SetTitle(Form(";%s;%s", titleX.Data(), titleY.Data()));
        histInfo.hist -> SetLineWidth(0);

        int xLabelNameNum = mLabelName.size();
        for(int b=0; b<xLabelNameNum; b++){
            int idx = mLabelName[b].first+1;
            TString name = mLabelName[b].second;
            if(idx < 1 || idx > histInfo.hist->GetNbinsX()){continue;}
            histInfo.hist -> GetXaxis()->SetBinLabel(idx, name);
        }   

        histInfo.base.cIdx = c;
        histInfo.base.objIdx = -1;
        mHist.push_back(histInfo);
    }
}

void DrawingUtil::InitGraph(int& cNum)
{
    cNum = -1;
    int graphNum = mGraph.size();
    for(int i=0; i<graphNum; i++){
        if(cNum < mGraph[i].base.cIdx+1){
            cNum = mGraph[i].base.cIdx+1;
        }
    }
    if(cNum <= 0){cNum = 1;}

    for(int c=0; c<cNum; c++){
        double minX = 99999.;
        double maxX = -99999.;
        double minY = 99999.;
        double maxY = -99999.;
        TString histName;
        TString titleX;
        TString titleY;

        for(int i=0; i<graphNum; i++){
            if(mGraph[i].base.cIdx == c){
                int pointNum = mGraph[i].graph -> GetN();
                for(int p=0; p<pointNum; p++){
                    double x, y;
                    mGraph[i].graph->GetPoint(p, x, y);
                    double ex = mGraph[i].graph->GetErrorX(p);
                    double ey = mGraph[i].graph->GetErrorY(p);
                    
                    if(minX > (x-ex)){minX = (x-ex);}
                    if(maxX < (x+ex)){maxX = (x+ex);}
                    if(minY > (y-ey)){minY = (y-ey);}
                    if(maxY < (y+ey)){maxY = (y+ey);}
                }
                
                TString tmpName, tmpXTitle, tmpYTitle;
                bool validTitle = GetTitle(mGraph[i].base.title, tmpName, tmpXTitle, tmpYTitle);
                if(validTitle){
                    titleX = tmpXTitle;
                    titleY = tmpYTitle;
                }
            }
        }

        if(GetDrawFlag("ratio")){
            minY = 0.;
            maxY = 1.5;
            if(titleY == ""){
                titleY = "Ratio";
            }
        }

        InfoGraph graphInfo;
        graphInfo.graph = new TGraphErrors();
        graphInfo.graph -> SetName(Form("baseGraph_%s_%i", titleX.Data(), c));
        graphInfo.graph -> SetPoint(graphInfo.graph->GetN(), minX, minY);
        graphInfo.graph -> SetPoint(graphInfo.graph->GetN(), maxX, maxY);
        graphInfo.graph -> SetMarkerSize(0);
        graphInfo.graph -> SetLineWidth(0);
        graphInfo.graph -> SetTitle(Form(";%s;%s", titleX.Data(), titleY.Data()));
        graphInfo.base.cIdx = c;
        graphInfo.base.objIdx = -1;
        mGraph.push_back(graphInfo);
    }
}

void DrawingUtil::DrawText(int cNum)
{
    TLatex* latex = new TLatex();

    int latexNum = mLatex.size();
    double textStartX = 0.18;
    double textStartY = 0.925;
    for(int i=0; i<latexNum; i++){
        if(mLatex[i].persistency){
            TString text = mLatex[i].text;
            if(text == ""){continue;}
            double x = mLatex[i].x;
            double y = mLatex[i].y;
            if(x < 0.){x = textStartX;}
            if(y < 0.){y = textStartY;}
            if(mLatex[i].size > 0.){latex->SetTextSize(mLatex[i].size);}
            if(mLatex[i].font > 0.){latex->SetTextSize(mLatex[i].font);}
            latex -> DrawLatexNDC(textStartX, textStartY, Form("%s", text.Data()));
            textStartY -= 0.055;
        }
    }
    for(int i=0; i<latexNum; i++){
        if(!mLatex[i].persistency && mLatex[i].cIdx == -1){
            TString text = mLatex[i].text;
            if(text == ""){continue;}
            double x = mLatex[i].x;
            double y = mLatex[i].y;
            if(x < 0.){x = textStartX;}
            if(y < 0.){y = textStartY;}
            if(mLatex[i].size > 0.){latex->SetTextSize(mLatex[i].size);}
            if(mLatex[i].font > 0.){latex->SetTextSize(mLatex[i].font);}
            latex -> DrawLatexNDC(textStartX, textStartY, Form("%s", text.Data()));
            textStartY -= 0.055;
        }
    }
    for(int i=0; i<latexNum; i++){
        if(!mLatex[i].persistency && mLatex[i].cIdx == cNum){
            TString text = mLatex[i].text;
            if(text == ""){continue;}
            double x = mLatex[i].x;
            double y = mLatex[i].y;
            if(x < 0.){x = textStartX;}
            if(y < 0.){y = textStartY;}
            if(mLatex[i].size > 0.){latex->SetTextSize(mLatex[i].size);}
            if(mLatex[i].font > 0.){latex->SetTextSize(mLatex[i].font);}
            latex -> DrawLatexNDC(textStartX, textStartY, Form("%s", text.Data()));
            textStartY -= 0.055;
        }
    }
}

void DrawingUtil::DrawLegend(int cNum)
{
    TLegend* persistantLegend = new TLegend(0.64, 0.68, 0.89, 0.965);
    persistantLegend -> SetBorderSize(0);
    persistantLegend -> SetTextSize(0.035);
    int legNum = mLegend.size();
    int canvasLegNum = 0;
    for(int i=0; i<legNum; i++){
        if(mLegend[i].cIdx == -1 || mLegend[i].persistency){
            persistantLegend -> AddEntry(mLegend[i].obj, Form("%s", mLegend[i].text.Data()), mLegend[i].opt);
        }
        if(mLegend[i].cIdx == cNum && !mLegend[i].persistency){
            canvasLegNum++;
        }
    }
    if(canvasLegNum == 0){
        persistantLegend -> Draw("same");
    }
    else{
        TLegend* leg = new TLegend(0.63, 0.68, 0.9, 0.965);
        leg -> SetBorderSize(0);
        leg -> SetTextSize(0.035);
        for(int i=0; i<legNum; i++){
            if(mLegend[i].cIdx = cNum){
                leg -> AddEntry(mLegend[i].obj, Form("%s", mLegend[i].text.Data()), mLegend[i].opt);
            }
        }
        leg -> Draw("same");
    }
}

bool DrawingUtil::InitDrawOption(TString opt)
{
    opt.ToUpper();
    opt.ReplaceAll(",", " ");
    opt.ReplaceAll(".", " ");
    opt.ReplaceAll(";", " "); 
    TObjArray *tokens = opt.Tokenize(" ");
    for(int i=0; i<tokens->GetEntries(); i++){
        TString option = ((TObjString *)tokens->At(i))->GetString();
        mDrawOptArr.push_back(option);
    }
}

bool DrawingUtil::GetDrawFlag(TString opt)
{   
    opt.ToUpper();
    int num = mDrawOptArr.size();
    for(int i=0; i<num; i++){
        TString option = mDrawOptArr[i];
        if(opt == "N" || opt == "NORM"){
            if(option == "N" || option == "NORM"){return true;}
        }
        if(opt == "YLOG" && option == "YLOG"){return true;}
        if(opt == "XLOG" && option == "XLOG"){return true;}
        if(opt == "ZLOG" && option == "ZLOG"){return true;}
        if(opt == "RATIO" && option == "RATIO"){return true;}
    }
    return false;
}

TString DrawingUtil::GetDrawValue(TString opt)
{
    int num = mDrawOptArr.size();
    for(int i=0; i<num; i++){
        TString option = mDrawOptArr[i];
        if(opt == "LINE" && option.Index("LINE")){
            TObjArray *tokens = option.Tokenize("=");
            return ((TObjString *)tokens->At(tokens->GetEntries()-1))->GetString();
        }
    }
}

bool DrawingUtil::GetTitle(TString inputTitle, TString& name, TString& xTitle, TString& yTitle)
{
    if(inputTitle == "" || inputTitle.Index(";")==-1){return 0;}

    TObjArray *tokens = inputTitle.Tokenize(";");
    int entries = tokens -> GetEntries();
    if(entries==0 && inputTitle != ""){
        xTitle = inputTitle;
        return 1;
    }
    if(entries==1){
        xTitle = ((TObjString *)tokens->At(0))->GetString();
        return 1;
    }
    if(entries==2){
        xTitle = ((TObjString *)tokens->At(0))->GetString();
        yTitle = ((TObjString *)tokens->At(1))->GetString();
        return 1;
    }
    else if(entries==3){
        name = ((TObjString *)tokens->At(0))->GetString();
        xTitle = ((TObjString *)tokens->At(1))->GetString();
        yTitle = ((TObjString *)tokens->At(2))->GetString();
        return 1;
    }
    return 0;
}

void DrawingUtil::NormalizedTH1D(TH1D* hist, double scale)
{
    double entries = hist -> GetEntries();
    if(entries < 1){return;}
    for(int bin=0; bin<hist->GetNbinsX(); bin++){
        double entry = hist->GetBinContent(bin+1);
        double binWidth = hist->GetXaxis()->GetBinWidth(bin+1);
        hist->SetBinContent(bin+1, entry/binWidth/entries);
    }
}
