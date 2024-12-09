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

1. Clone [Vorbis library](https://github.com/xiph/vorbis/tree/master) to your local directory.

2. Build the library.

3. Find the directory called `examples`

4. Modify Makefile to add the flag `-logg`

5. Change line 127 to the suitable value. If you don't know the information of the original audio file, install `oggz` then run:
```bash=
oggz info [inputfile]
```

6. Compile `encoder_example.c`

7. The usage of sample encoder is:
```bash=
./encoder_example < [inputfile] > [outputfile]
```

8. Then use `ogg123` to play:
```bash=
ogg123 [inputfile]
```
