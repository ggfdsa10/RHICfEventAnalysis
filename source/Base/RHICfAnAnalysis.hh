#ifndef RHICfAnAnalysis_hh
#define RHICfAnAnalysis_hh

#include "TString.h"

class RHICfAnAnalysis
{
    public:
        RHICfAnAnalysis();
        ~RHICfAnAnalysis();

        int Init();
        int Calculate();
        int SaveData();

        // 
        // void Ontest();

    private:

};

#endif
