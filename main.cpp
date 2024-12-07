#include <bits/stdc++.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "io.hpp"
#include "header.hpp"

using namespace std;

void parse(io_buf &in, FILE* fout) {
	identification id;
	comment		   cc;
	setup		   ss;
	int h_count = 0;
	in.parse();

	while(h_count < 3) {
		in.new_p();
		unsigned int type = in.read_u(8); // header type
		//cout << type << endl;
		for(int i = 0; i < 6; ++i) in.read_u(8); // vorbis
		if(type == 1) id.init(in);
		if(type == 3) cc.init(in);
		if(type == 5) ss.init(in, id);
		h_count++;
		in.padding();
	}
	packet p(id);
	int count = 0;
	int total = 0;
	while(1) {
		in.new_p();
		bool type = in.read_u(1); // header type
		if(type != 0) exit(-4);
		p.decode(in, id, ss);
		for(int i = 0; i < id.audio_channels; ++i) {
			for(int j = 0; j < p.Y_ret[0].size(); ++j) { 
				float ff = p.Y_ret[i][j] * 32767.f + .5f;
				int s = floor(ff);
				if(s > 32767) s = 32767;
				if(s < -32768) s = -32768;
				short s2 = s;
				fwrite(&s2, sizeof(short), 1, fout);
			}
		}
		total += p.Y_ret[0].size();
		//printf("%d: size%d\nread%d\n", count + 3, p.Y_ret[0].size(), in.count);
		count++;
		in.padding();
		if(in.qq.empty()) break;
	}
	/*
	cout << total << endl;
	cout << count << endl;
	cout << in.ptr << " " << in.len << endl;
	cout << in.now.size() << endl;
	*/
}

int main(int argc, char *argv[]) {
	/*
	 * open file
	 */
	int fin = open(argv[1], O_RDWR, (mode_t) 0644);
	FILE *fout = fopen(argv[2], "wb");

	/*
	 * set io buffer
	 */
	struct stat st;
	stat(argv[1], &st);
	int fin_siz = st.st_size;
	unsigned char *in_map = (unsigned char*) mmap(NULL, fin_siz, PROT_READ, MAP_SHARED, fin, 0);
	io_buf input(in_map, fin_siz);
	parse(input, fout);
	munmap(in_map, fin_siz);
}
