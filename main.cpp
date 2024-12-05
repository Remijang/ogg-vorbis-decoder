#include <bits/stdc++.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "io.hpp"
#include "header.hpp"

using namespace std;

void parse(io_buf &in) {
	vector<unsigned char> q;
	identification id;
	comment		   cc;
	setup		   ss;
	for(int i = 0; i < 7; ++i) q.push_back(in.read_u(8));
	while(1) {
		if(q[1] == 'v' && q[2] == 'o' && q[3] == 'r' && q[4] == 'b' && q[5] == 'i' && q[6] == 's'){
			cout << "abc" << " " << (int)q[0] << endl;
			if(q[0] == 1) id.init(in);
			if(q[0] == 3) cc.init(in);
			if(q[0] == 5) ss.init(in, id);
			q.clear();
			for(int i = 0; i < 7; ++i) q.push_back(in.read_u(8));
			continue;
		}
		q.push_back(in.read_u(8));
		q.erase(q.begin());
		if(in.ptr >= 20000) break;
	}
	cout << cc.vender_string.data() << endl;
	for(unsigned int i = 0; i < cc.user_comment_list_length; ++i){
		cout << cc.user_comment_list[i].data() << endl;
	}
	packet p;
	int count = 0;
	while(in.ptr < in.len) {
		cout << count << endl;
		p.decode(in, id, ss);
		count++;
	}
	cout << in.ptr << " " << in.len << endl;
	cout << in.offset << endl;
}

int main() {
	/*
	 * open file
	 */
	int fin = open("test.ogg", O_RDWR, (mode_t) 0644);

	/*
	 * set io buffer
	 */
	struct stat st;
	stat("test.ogg", &st);
	int fin_siz = st.st_size;
	unsigned char *in_map = (unsigned char*) mmap(NULL, fin_siz, PROT_READ | PROT_WRITE, MAP_SHARED, fin, 0);
	io_buf input(in_map, fin_siz);

	parse(input);

}
