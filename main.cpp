#include "test.h"
#include "coproto/Socket/AsioSocket.h"
#include "coproto/Socket/Socket.h"
#include "libOTe/config.h"
#include "Coeff128.h"
#include "libOTe/Triple/Foleage/fft/FoleageFft.h"
#include <fstream>
#include <ctime>
#include <omp.h>
#include "libOTe/Tools/CoeffCtx.h"
#include "OLE.h"
#include "F4.h"
#include <bitset>
//#define trial 2
//#define num_thread 1
using namespace osuCrypto;


void printUsage() {
    printf("Usage: ./main [OPTIONS] n\n");
    printf("Options:\n");
    printf("  --QA_Syndrome\tTests syndrome encoding of QA code.\n");
    printf("  --EA_Syndrome\tTests syndrome encoding of EA code.\n");
    printf("  --EC_Syndrome\tTests syndrome encoding of EC code.\n");
    printf("  --QA_VOLE\tTests VOLE based on QA code.\n");
    printf("  --EA_VOLE\tTests VOLE based on EA code.\n");
    printf("  --EC_VOLE\tTests VOLE based on EC code..\n");
    printf("  --OLE\tTests OLE based on QA code.\n");
    printf("  --QA_PCS\tTests  encoding of QA code for PCS.\n");
    printf("  --EA_PCS\tTests  encoding of EA code for PCS.\n");
    printf("  --RAA_PCS\tTests  encoding of RAA code for PCS.\n");
}

   int main(int argc, char **argv)
{
        if (argc <3)
        {
        printUsage();
        }
        else
        {
        int num_var=std::atoi(argv[2]);
        int n=ipow(3,num_var);
        if (strcmp(argv[1], "--QA_Syndrome") == 0) {
            QA_encode_test(num_var,n);
        } else if (strcmp(argv[1], "--EA_Syndrome") == 0) {
             EA_encode_test<block, CoeffCtxGF128>(n, 5, 21);
        } else if (strcmp(argv[1], "--EC_Syndrome")==0) {
            EC_encode_test<block, CoeffCtxGF128>(n,2,7,24,true); 
        } else if (strcmp(argv[1], "--QA_VOLE") == 0) {
            VOLE_QASD<block, CoeffCtxGF128_F4>(n);
        } else if (strcmp(argv[1], "--EC_VOLE") == 0) {
             Vole_LPN<block, block, CoeffCtxGF128>(n, osuCrypto::MultType::ExConv7x24, false, false, false);
  
        } else if (strcmp(argv[1], "--EA_VOLE") ==0) {
            Vole_LPN<block, block,CoeffCtxGF128>(n, osuCrypto::MultType::ExAcc21, false, false, false);
        } else if (strcmp(argv[1], "--OLE") ==0) {
           Extension_Field_OLE<block, CoeffCtxGF128_F4>(num_var, 4, 2);
           Extension_Field_OLE<block, CoeffCtxGF128_F4>(num_var, 3, 3);
        }
        else if(strcmp(argv[1], "--QA_PCS") == 0) {
             QA_encode_test_PCS(num_var,n,2);
             QA_encode_test_PCS(num_var,n,3);
             QA_encode_test_PCS(num_var,n,4);
             QA_encode_test_PCS(num_var,n,5);
        }
        else if(strcmp(argv[1], "--EA_PCS") == 0) {
            int logn=log2ceil(n);
            int inverse_rate=2;
            int nonzerocolumn=(logn*18)/inverse_rate;
            EA_encode_test_PCS<block, CoeffCtxGF128>(n, nonzerocolumn, inverse_rate);
        }
        else if(strcmp(argv[1], "--RAA_PCS") == 0) {
            RAA_encode_test_PCS<block, CoeffCtxGF128>(n, 2);
            RAA_encode_test_PCS<block, CoeffCtxGF128>(n, 4);
            RAA_encode_test_PCS<block, CoeffCtxGF128>(n, 8);
        }
            

        else
        {
            printUsage();
        }
    }

     
}

  