#include <bits/stdc++.h>
using namespace std;

#ifndef ONLINE_JUDGE
#include "algo/debug.h"
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif

#define int long long

const int LOGN = 22;
/*
Revise When superFree.
*/
void solve() {
	int n;
	cin >> n;
	vector<vector<int>> g(n + 1);

	for (int i = 1, u, v; i < n; i++) {
		cin >> u >> v;
		g[u].push_back(v);
		g[v].push_back(u);
	}

	vector<vector<int>> dp(n + 1, vector<int>(LOGN, 0ll));
	vector<int> depth(n + 1, -1ll), in(n + 1, 0ll); // to normalise with 0
	int timer = 0;

	function<void (int, int)> dfs = [&](int node, int par = 0) {
		depth[node] = depth[par] + 1;
		dp[node][0] = par;
		in[node] = timer++;

		//node-> 2^0 par
		for (int i = 1; i < LOGN; i++) {
			dp[node][i] = dp[dp[node][i - 1]][i - 1];
		}

		for (auto child : g[node]) {
			if (child == par) {
				continue;
			}

			dfs(child, node);
		}
	};

	function<int (int, int)> ancestor = [&](int node, int k) {
		//par k distnace up !
		for (int i = 0; i < LOGN && node; i++) {
			if (k & (1ll << i)) {
				node = dp[node][i];
			}
		}
		return (node == 0 ? -1 : node);
	};

	function<int (int, int)> lca = [&](int node1, int node2) {
		if (depth[node1] < depth[node2]) {
			swap(node2, node1);
		}

		node1 = ancestor(node1, depth[node1] - depth[node2]);

		if (node1 == node2)return node1;


		for (int i = LOGN - 1; i >= 0; i--) {
			if (dp[node1][i] != dp[node2][i]) {
				node1 = dp[node1][i];
				node2 = dp[node2][i];
			}
		}

		return dp[node1][0];
	};
	dfs(1, 0);

	int q;
	cin >> q;

	vector<vector<pair<int, int>>> vt(n + 1);//virtual-tree

	vector<pair<int, int>> dp_ans(n + 1, {0ll, 0ll});
	vector<bool> special(n + 1, false);
	bool ok = 1;

	function<void(int, int, int)> dfs_ans = [&](int node, int par, int ww) {
		if (!ok) {
			return;
		}
		for (auto [ch, w] : vt[node]) {
			if (ch ^ par) {
				dfs_ans(ch, node, w);
			}
		}

		if (special[node]) {
			dp_ans[node].second = 1;
			for (auto [ch, w] : vt[node]) {
				if (ch ^ par) {
					if (dp_ans[ch].second) {
						//there is a special node in the cut
						if (w == 1) {
							ok = false;
							return;
						} else {
							dp_ans[node].first += dp_ans[ch].first + 1;
						}
					} else {
						dp_ans[node].first += dp_ans[ch].first;
					}
				}
			}
		} else {
			int cnt = 0, sm = 0;
			for (auto [ch, w] : vt[node]) {
				if (ch ^ par) {
					if (dp_ans[ch].second) {
						cnt++;
					}
					sm += dp_ans[ch].first;
				}
			}
			if (cnt >= 2) {
				dp_ans[node] = {1, 0};
			} else {
				dp_ans[node] = {0, (cnt ? 1 : 0)};
				if (ww == 1) {
					if (special[par]) {
						dp_ans[node] = {(cnt > 0), 0};
					}
				}
			}

			dp_ans[node].first += sm;
		}

	};


	while (q--) {
		int k;
		cin >> k;
		ok = 1;
		vector<int> c(k);
		for (auto &x : c) {
			cin >> x;
			special[x] = 1;
		}

		vector<pair<int, int>> ne;

		for (auto &x : c) {
			ne.push_back({in[x], x});
		}
		sort(ne.begin(), ne.end());



		for (int i = 0; i + 1 < k; i++) {
			int z = lca(ne[i].second, ne[i + 1].second);
			ne.push_back({in[z], z});
		}



		set<pair<int, int>> s(ne.begin(), ne.end());


		ne.clear();
		for (auto &x : s) {
			ne.push_back(x);
		}


		vector<int> cur;
		cur.push_back(1);

		for (int i = (ne.front().second == 1 ? 1 : 0); i < ne.size(); i++) {
			while (1) {
				auto x = cur.back();
				if (lca(x, ne[i].second) == x) {
					break;
				} else {
					cur.pop_back();
				}
			}
			auto x = cur.back();
			int w = depth[ne[i].second] - depth[x];
			vt[x].push_back({ne[i].second, w});
			vt[ne[i].second].push_back({x, w});

			cur.push_back(ne[i].second);
		}


		dfs_ans(1, 0, 0);

		if (!ok) {
			cout << -1 << "\n";
		} else {
			cout << dp_ans[1].first << "\n";
		}
		dp_ans[1] = {0, 0};

		vt[1].clear();

		for (auto [inT, node] : ne) {
			dp_ans[node] = {0, 0};
			special[node] = false;
			vt[node].clear();
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
