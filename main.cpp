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
	int h_count = 0;
	while(h_count < 3) {
		if(q[1] == 'v' && q[2] == 'o' && q[3] == 'r' && q[4] == 'b' && q[5] == 'i' && q[6] == 's'){
			cout << "abc" << " " << (int)q[0] << endl;
			if(q[0] == 1) id.init(in);
			if(q[0] == 3) cc.init(in);
			if(q[0] == 5) ss.init(in, id);
			h_count++;
			if(h_count < 3) {
				q.clear();
				for(int i = 0; i < 7; ++i) q.push_back(in.read_u(8));
				continue;
			}
			else break;
		}
		q.push_back(in.read_u(8));
		q.erase(q.begin());
	}
	/*
	cout << cc.vender_string.data() << endl;
	for(unsigned int i = 0; i < cc.user_comment_list_length; ++i){
		cout << cc.user_comment_list[i].data() << endl;
	}
	*/
	packet p(id);
	int count = 0;
	while(in.ptr < in.len) {
		p.decode(in, id, ss);
		count++;
	}
	cout << count << endl;
	cout << in.ptr << " " << in.len << endl;
	cout << in.offset << endl;
}

int main(int argc, char *argv[]) {
	/*
	 * open file
	 */
	int fin = open(argv[1], O_RDWR, (mode_t) 0644);

	/*
	 * set io buffer
	 */
	struct stat st;
	stat(argv[1], &st);
	int fin_siz = st.st_size;
	unsigned char *in_map = (unsigned char*) mmap(NULL, fin_siz, PROT_READ, MAP_SHARED, fin, 0);
	{
		io_buf input(in_map, fin_siz);

		parse(input);
	}
	munmap(in_map, fin_siz);
}
