#ifndef RHICfDLECondition_hh
#define RHICfDLECondition_hh

#include "RHICfOptContainer.hh"

class RHICfDLECondition
{
    public:
        RHICfDLECondition();
        ~RHICfDLECondition();

        void Init();

        void SetBTofMult(int mult);
        void SetBBCADC(int smallE, int smallW, int largeE, int largeW);

        bool isSDLE();
        bool isDDLE();
        bool isNDLE();

        int GetDLEIdx();

        int GetSimEventType(int model, int id); // convert from id to event type

    private:
        RHICfOptContainer* mOptContainer;

        int mBTofMult;
        int mBBCSmallEast;
        int mBBCSmallWest;
        int mBBCLargeEast;
        int mBBCLargeWest;

        static const int mBTofMult_thr = 1;
        static const int mBBCSmallEast_thr = 40;
        static const int mBBCSmallWest_thr = 40;
        static const int mBBCLargeEast_thr = 6;
        static const int mBBCLargeWest_thr = 4;

};

#endif
