# Extension Field PCG
The PCG implementation based on FFT transform over F2^128.

## Installing libOTe
Our code is developed based on libOTe [libOTe](https://github.com/osu-crypto/libOTe). To install libOTe library, run the following command.

```
cd libOTe
python build.py --all --install --boost --sodium --openssl -D FETCH_SODIUM=OFF -D SODIUM_MONTGOMERY=FALSE
```
If installation requires sudo access, run the following command

```
cd libOTe
python build.py --all --install --sudo --boost --sodium --openssl -D FETCH_SODIUM=OFF -D SODIUM_MONTGOMERY=FALSE
```
Please manually delete it if a previous version of libOTe is installed in your computer.

## Running benchmarks
Benchmarks:

```
cmake .
make
./main   --QA_Syndrome n Tests syndrome encoding of QA code of length 3^n
./main   --EA_Syndrome n Tests syndrome encoding of EA code of length 3^n
./main   --EC_Syndrome n Tests syndrome encoding of EC code of length 3^n
./main   --QA_VOLE n Tests VOLE boased on QA code of length 3^n
./main   --EA_VOLE n Tests VOLE boased on EA code of length 3^n
./main   --EC_VOLE n Tests VOLE boased on EC code of length 3^n
./main   --OLE n Tests OLE boased on QA code of length 3^n
./main   --QA_PCS n Tests encoding of QA code for PCS of msssage length 3^n
./main    --RAA_PCS n Tests encoding of RAA code for PCS of msssage length 3^n
./main    --EA_PCS n Tests encoding of EA code for PCS of msssage length 3^n
```
