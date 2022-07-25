#include <bits/stdc++.h>
using namespace std;
/*
Given a DAG G, delete some edges such that indegree and outdegree of each node decreases or remain 0, now find the largest set S, 
such that there is a path from a to b or b to a for each a, b in the set.

Now this hints towards top sort, after top sort, drawing graph in straight line, such that each edge goes from left to right, 

Now we can prove that this largest set S, is nothing but a simple path from a-> b -> c -> d...
*/
int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#else
#endif
	ios::sync_with_stdio(0); cin.tie(0);

	int n, m;
	cin >> n >> m;

	vector<vector<int>> g(n + 1), rg(n + 1);
	vector<int> in(n + 1, 0ll), out(n + 1, 0ll);
	for (int i = 0, a, b; i < m; i++) {
		cin >> a >> b;
		g[a].push_back(b);
		rg[b].push_back(a);
		in[b]++;
		out[a]++;
	}
	queue<int> q;
	for (int i = 1; i <= n; i++) {
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
			if (!in[ch]) {
				q.push(ch);

			}
		}
	}

	vector<int> dp(n + 1, 0);
	for (int i = 0; i < n; i++) {
		vector<int> pos;
		bool ok = 0;
		for (auto x : rg[tp[i]]) {
			if (out[x] == 1) {
				ok = 1;
			} else {
				pos.push_back(dp[x]);
			}
		}
		if (ok) {
            //If ok we have already deleted a in degree wala edge
            // So just take the maximum
			pos.push_back(0);
			dp[tp[i]] = 1 + *max_element(pos.begin(), pos.end());
		} else {
			if (pos.size() > 1) {
                // if pos.size > 1 then just take the max
				dp[tp[i]] = 1 + *max_element(pos.begin(), pos.end());
			} else {
                //o.w its just 1
				dp[tp[i]] = 1;
			}
		}
	}

	cout << *max_element(dp.begin(), dp.end()) << "\n";
}
