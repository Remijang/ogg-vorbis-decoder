#ifndef UTIL_HPP
#define UTIL_HPP

#include <math.h>

int ilog(int x) {
	int ret = 0;
	while(x > 0) {
		ret++;
		x >>= 1;
	}
	return ret;
}

int ilog(unsigned int x) {
	int ret = 0;
	while(x > 0) {
		ret++;
		x >>= 1;
	}
	return ret;
}

float float32_unpack(unsigned int x) {
	int mantissa = x & 0x1fffff;
	int sign = x & 0x80000000;
	int exponent = (x & 0x7fe00000) >> 21;
	if(sign != 0) mantissa = ~mantissa;
	return mantissa * pow(2, exponent - 788);
}

unsigned int lookup1_values(unsigned int entries, unsigned int dimensions) {
	unsigned int ret = 1;
	while(1){
		unsigned long long tmp = 1;
		for(unsigned int i = 0; i < dimensions; ++i){
			tmp *= ret;
			if(tmp > 0xffffffff) break;
		}
		if(tmp > entries) break;
		ret++;
	}
	return ret - 1;
}

double bark(double x) {
	return 13.1 * atan(0.0074 * x) + 2.24 * atan(0.0000000185 * x * x) + 0.0001 * x;
}

#endif
