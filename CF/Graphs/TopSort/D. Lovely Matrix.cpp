#include <bits/stdc++.h>
using namespace std;

#ifndef ONLINE_JUDGE
#include "algo/debug.h"
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif

#define int long long

/*
Super Hard and intresting problem, So given an arrays sort each array in ascending order, you can swap all columns only, like choose index i and j
and swap them. find if it is possible and plausible sort.

for all bunch of 1 and bunch of 2, we know that bunch of 1 should be before 2, so just make an edge from all 1s to a dummy node and make an edge from 
dummy node to bunch of 2.

This ensures 1 -> dummy -> 2 

So we can create a graph, Now top sort is the order so that all the nodes x should come before y,  x --> (z1..zk) --> y {essentially top sort defination}
*/
void solve() {
	int n, m;
	cin >> n >> m;
	int dum = m + 1;
	vector<vector<int>> g(m + 1);
	vector<int> in(m + 1, 0ll);
	for (int i = 0; i < n; i++) {
		map<int, vector<int>> mp;
		for (int j = 0, x; j < m; j++) {
			cin >> x;
			if (x == -1) {
				continue;
			}
			mp[x].push_back(j + 1);
		}

		for (auto it = mp.begin(); it != mp.end(); it++) {
			auto it2 = it;
			it2++;
			if (it2 == mp.end()) {
				break;
			}

			g.push_back(vector<int> {});
			in.push_back(0);

			for (auto id : it->second) {
				g[id].push_back(dum);
				in[dum]++;
			}

			for (auto id : it2->second) {
				g[dum].push_back(id);
				in[id]++;
			}
			dum++;
		}
	}
	debug(g, in.size(), dum);
	queue<int> q;
	for (int i = 1; i < dum; i++) {
		if (in[i] == 0) {
			q.push(i);
		}
	}
	vector<int> tp;
	while (q.size()) {
		int node = q.front();
		q.pop();
		tp.push_back(node);

		for (auto ch : g[node]) {
			in[ch]--;
			if (in[ch] == 0) {
				q.push(ch);
			}
		}
	}
	if (tp.size() != dum - 1) {
		cout << "-1\n";
	} else {
		for (auto x : tp) {
			if (x > m) {
				continue;
			}
			cout << x << " ";
		}
	}
}

signed main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#else
#endif
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);

	int t = 1;
	//cin >> t;

	while (t--) {
		solve();
	}

	// cerr << "Time elapsed: " << ((long double)clock() / CLOCKS_PER_SEC) << " s.\n";
}

