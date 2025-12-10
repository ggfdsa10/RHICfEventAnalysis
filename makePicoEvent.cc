// RHICf Pico Event generator
#include "RHICfPicoEventMaker.hh"

using namespace std;

int main(int argc, char** argv)
{
    RHICfPicoEventMaker* picoEventMaker = new RHICfPicoEventMaker();
    picoEventMaker -> SetInputDataPath("/gpfs01/star/pwg/slee5/RHICfEvent/rhicf"); 
    picoEventMaker -> SetRunType("TOP"); 
    // picoEventMaker -> SetExecuteFileNum(5); 

    // picoEventMaker -> MakeRHICfStream();
    // picoEventMaker -> MakePhysicsStream();
    picoEventMaker -> MakeRHICfParticleEvent();

    picoEventMaker -> CalculateNeutron();  
    // picoEventMaker -> CalculatePi0();  

    picoEventMaker -> Init();
    picoEventMaker -> Calculate();
    picoEventMaker -> Finish();

    cout << "Terminated." << endl;
    
    return 0;
}
