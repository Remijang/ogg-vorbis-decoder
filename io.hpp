#ifndef IO_HPP
#define IO_HPP

#include <stdlib.h>
#include <vector>
#include <queue>

using namespace std;

struct io_buf {
	unsigned char* buf;
	int ptr;
	int offset;
	int len;
	queue<unsigned int> seg;
	queue<unsigned char> now;
	queue<unsigned int> qq;
	int quota;
	int end;
	int count;
	int end_p;

	io_buf (unsigned char *in, int n): buf(in), len(n) {
		ptr = offset = -1;
		quota = count = end = 0;
	}

	void parse() {
		ptr = 0;
		now.push(0);
		while(end != 1) {
			end = end == 1 ? 1 : buf[ptr + 5] >> 2;
			ptr += 26;
			int jmp = 0;
			int count = buf[ptr++];
			for(int i = 0; i < count; ++i) {
				int tmp = buf[ptr++];
				seg.push(tmp);
				jmp += tmp;
			}
			for(int i = 0; i < jmp; ++i) now.push(buf[ptr++]);
		}
		while(!seg.empty()) {
			int sum = 0, tmp;
			do {
				tmp = seg.front();
				seg.pop();
				sum += tmp;
			} while(!seg.empty() && tmp == 255);
			qq.push(sum);
		}
	}

	void padding() {
		offset = -1;
		while(quota > 0) {
			now.pop();
			quota--;
		}
	}

	void new_p() {
		offset = -1;
		end_p = 0;
		count = 0;
		quota = qq.front();
		qq.pop();
	}

	int getbit() {
		if(quota < 0) return 0;
		offset++;
		if((offset & 7) == 0) {
			quota--;
			if(quota < 0) {
				end_p = 1;
				return 0;
			}
			count++;
			now.pop();
			return now.front() & 1;
		}
		else return now.front() >> (offset & 7) & 1;
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
