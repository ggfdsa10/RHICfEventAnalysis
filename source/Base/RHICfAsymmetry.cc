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
}

void RHICfAsymmetry::InitData()
{
    for(int fill=0; fill<kFillNum; fill++){
        for(int type=0; type<kTypeNum; type++){
            for(int dle=0; dle<kDLENum; dle++){
                int runIdx = mOptContainer -> GetFillToRunIdx(fill);
                int ptNum = mBinning -> GetPtBinNum(runIdx, type, dle);
                int xfNum = mBinning -> GetXfBinNum(runIdx, type, dle);

                if(ptNum > 0 && xfNum > 0){
                    mRightAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mLeftAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mBkgRightAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));
                    mBkgLeftAsymmetry[fill][type][dle].resize(ptNum, vector<double>(xfNum));   
                }
            }
        }
    }
}

void RHICfAsymmetry::SetBinning(RHICfBinning* binning)
{
    mBinning = binning;
}

void RHICfAsymmetry::SetDilution(RHICfDilutionFactor* dilution)
{
    mDilution = dilution;
}

void RHICfAsymmetry::FillPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand)
{
    if(isRightHand){
        mRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
    }
    else{
        mLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
    }
}

void RHICfAsymmetry::FillBkgPolarization(int fillIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, bool isRightHand)
{
    if(isRightHand){
        mBkgRightAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
    }
    else{
        mBkgLeftAsymmetry[fillIdx][typeIdx][dleIdx][ptIdx][xfIdx] += 1.;
    }
}

void RHICfAsymmetry::Calculate()
{

}

void RHICfAsymmetry::SaveAsymmetryData()
{

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

