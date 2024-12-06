#ifndef IO_HPP
#define IO_HPP

#include <stdlib.h>
#include <vector>

using namespace std;

struct io_buf {
	unsigned char* now;
	int ptr;
	int offset;
	int len;

	io_buf (unsigned char *in, int n): now(in), len(n) {
		ptr = offset = -1;
	}

	void padding() {
		offset = -1;
	}

	int getbit() {
		offset++;
		if((offset & 7) == 0) {
			ptr++;
			return now[ptr] & 1;
		}
		else return now[ptr] >> (offset & 7) & 1;
	}

	unsigned int read_u(int n) {
		unsigned int ans = 0;
		for(int i = 0; i < n; ++i) ans |= (getbit() << i);
		return ans;
	}

	int read_s(int n) {
		int ans = 0;
		for(int i = 0; i < n - 1; ++i) ans |= (getbit() << i);
		int msb = getbit();
		if(n < 31 && msb == 1){
			ans = (1 << n) - (ans | (msb << (n - 1)));
		}
		else ans |= (msb << (n - 1));
		return ans;
	}

	vector<unsigned char> read_utf() {
		vector<unsigned char> v;
		unsigned char c;
		c = read_u(8);
		v.insert(v.begin(), c);
		if((c & 0x80) == 0) return v;
		c = read_u(8);
		v.insert(v.begin(), c);
		if((c & 0xE0) == 0xC0) return v;
		c = read_u(8);
		v.insert(v.begin(), c);
		if((c & 0xF0) == 0xE0) return v;
		c = read_u(8);
		v.insert(v.begin(), c);
		if((c & 0xF8) == 0xF0) return v;
		exit(-1);
	}

};


#endif
