# OGG Vorbis Decoder
Author: b11902051 巫俋霆

## Overview

This is a simple decoder based on Ogg Vorbis I spec.

The output will be 16-bit PCM stream, which can used by simple example encoder from Vorbis.

## Remark

In Huffman decode, the lookup time complexity is $\mathcal O(n)$, just brute force.

In MDCT, I use the easiest formula to implement, which has $\mathcal O(n^2)$ time complexity.

The buttleneck of program should be MDCT module, but I have no enough time to implement fast algorithm.

## Compile

Using Makefile:
```bash=
make
```

or compile by `g++` command:
```bash=
g++ -std=c++17 -O3 -Wall -o main main.cpp
```

## Usage

```bash=
./main [input file] [outputfile]
```

## For those who want to test the result

1. Clone [Vorbis library: encoder example](https://github.com/xiph/vorbis/blob/master/examples/encoder_example.c) to your local directory.

2. Change line 127 of `encoder_example.c` to the suitable value. If you don't know the information of the original audio file, install `oggz` then run:
```bash=
oggz info [inputfile]
```

3. Compile `encoder_example.c`, which can be done by:
```base=
gcc encoder_example.c -o encoder_example -O3 -logg -lvorbis -lvorbisenc
```

4. The usage of sample encoder is:
```bash=
./encoder_example < [inputfile] > [outputfile]
```

5. Then you can use `ogg123` to play:
```bash=
ogg123 [inputfile]
```
