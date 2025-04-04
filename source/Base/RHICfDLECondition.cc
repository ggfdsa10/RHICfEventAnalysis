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
