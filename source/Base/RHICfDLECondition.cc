#include "RHICfDLECondition.hh"

RHICfDLECondition::RHICfDLECondition() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
}

RHICfDLECondition::~RHICfDLECondition()
{
}

void RHICfDLECondition::Init()
{
    cout << "RHICfDLECondition::Init() -- Done." << endl;
}

void RHICfDLECondition::SetBTofMult(int mult)
{
    mBTofMult = mult;
}

void RHICfDLECondition::SetBBCADC(int smallE, int smallW, int largeE, int largeW)
{
    mBBCSmallEast = smallE;
    mBBCSmallWest = smallW;
    mBBCLargeEast = largeE;
    mBBCLargeWest = largeW;
}

bool RHICfDLECondition::isSDLE()
{
    if(mBTofMult < mBTofMult_thr){
        if(mBBCSmallEast < mBBCSmallEast_thr && mBBCLargeEast < mBBCLargeEast_thr){
            return true;
        }
    }
    return false;
}

bool RHICfDLECondition::isDDLE()
{
    if(mBTofMult < mBTofMult_thr){
        if(mBBCSmallEast >= mBBCSmallEast_thr && mBBCLargeEast >= mBBCLargeEast_thr){
            return true;
        }
    }
    return false;
}

bool RHICfDLECondition::isNDLE()
{
    if(mBTofMult >= mBTofMult_thr){
        if(mBBCSmallEast >= mBBCSmallEast_thr && mBBCLargeEast >= mBBCLargeEast_thr){
            if(mBBCSmallWest >= mBBCSmallWest_thr && mBBCLargeWest >= mBBCLargeWest_thr){
                return true;
            }
        }
    }
    return false;
}

int RHICfDLECondition::GetDLEIdx()
{
    if(isSDLE()){return kSDLE;}
    if(isDDLE()){return kDDLE;}
    if(isNDLE()){return kNDLE;}
    return kALLDLE;
}

int RHICfDLECondition::GetSimEventType(int model, int id)
{
    if(model == rPythia8){
        if(id == 101){return kND;}
        if(id == 103 || id == 104){return kSD;}
        if(id == 105){return kDD;}
    }
    else{
        if(id == 0){return -1;}
        id = abs(id);
        if(id >= 10){id -= 10;}
        if(id == 4){return kSD;}
        if(id == 1){return kND;}
        if(id == 2){return kDD;}
    }
    return -1;
}
