
#include "RHICfSimAnalysis.hh"

using namespace std;

int main(int argc, char** argv)
{
    RHICfSimAnalysis* sim = new RHICfSimAnalysis();
    // sim -> SetRunType("TS");
    sim -> CalculatePi0();  
    sim -> Init();
    sim -> Calculate();
    sim -> Finish();

    cout << "Terminated." << endl;
    return 0;
}
