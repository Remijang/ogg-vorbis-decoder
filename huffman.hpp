#ifndef HUFFMAN_HPP
#define HUFFMAN_HPP

#include <queue>
#include <vector>

struct node {
	int val;
	node *left, *right;

	node(int n): val(n) {
		left = right = NULL;
	}

	node(int n, node *a, node *b): val(n), left(a), right(b) {}
};

struct cmp {
	bool operator() (node *a, node *b) {
		return a->val > b->val;
	}
};

void huffman(int n, int *length, long long *codeword) {
	priority_queue <node *, vector<node *>, cmp> p[33];
	for(int i = 0; i < n; ++i) if(length[i] > 0) {
		node* tmp = new node(i);
		p[length[i]].push(tmp);
	}
	for(int l = 32; l > 0; --l) {
		while(!p[l].empty()) {
			if(p[l].size() % 2 != 0) exit(-1);
			node *n1 = p[l].top();
			p[l].pop();
			node *n2 = p[l].top();
			p[l].pop();
			if(n1 == NULL || n2 == NULL) exit(-2);
			node *tmp = new node(n1->val, n1, n2);
			p[l - 1].push(tmp);
		}
	}

	if(p[0].size() != 1) exit(-3);

	auto dfs = [&] (auto&& self, node *cur, long long s) -> void {
		if(cur->left == NULL && cur->right == NULL) {
			codeword[cur->val] = s;
			delete cur;
			return;
		}
		self(self, cur->left, s << 1);
		self(self, cur->right, (s << 1) + 1);
		delete cur;
	};

	dfs(dfs, p[0].top(), 0);
}


#endif
