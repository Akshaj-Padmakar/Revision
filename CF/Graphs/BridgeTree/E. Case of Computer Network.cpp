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
Given a Graph G, and M instructions, i.e orient the edges from path u to v or v to u.
So Is it possible or not ?
*/

const int LOGN = 22;
struct dsu {
public:
	dsu(int _n) {
		n = _n;
		par.assign(n, 0);
		sz.assign(n, 1);
		for (int i = 0; i < n; i++) {
			par[i] = i;
		}
	}

	int find(int v) {
		assert(v >= 0 && v < n);
		if (v == par[v]) {
			return v;
		}
		return par[v] = find(par[v]);
	}

	void merge(int u, int v) {
		assert(u >= 0 && u < n && v >= 0 && v < n);

		u = find(u), v = find(v);
		if (u == v) {
			return;
		}
		if (sz[u] > sz[v]) {
			swap(u, v);
		}
		par[u] = v;
		sz[v] += sz[u];
	}

	int size(int u) {
		assert(u >= 0 && u < n);
		u = find(u);
		return sz[u];
	}

	bool same(int u, int v) {
		assert(u >= 0 && u < n && v >= 0 && v < n);

		return find(u) == find(v);
	}

	void groups(vector<vector<int>> &grp) {
		vector<vector<int>> temp(n);
		for (int i = 0; i < n; i++ ) {
			temp[find(i)].push_back(i);
		}

		for (int i = 0; i < n; i++) {
			if (temp[i].size() >= 1) {
				grp.push_back(temp[i]);
			}
		}
	}

private:
	int n;
	vector<int> par, sz;
};
void solve() {
	int n, m, q;
	cin >> n >> m >> q;

	vector<vector<int>> g(n + 1);
	vector<pair<int, int>> e;
	dsu dsg(n + 1);
	for (int i = 0, a, b; i < m; i++) {
		cin >> a >> b;
		g[a].push_back(b);
		g[b].push_back(a);
		e.push_back({a, b});
		dsg.merge(a, b);
	}

	vector<int> dpIN(n + 1, 0ll), in(n + 1, 0ll), par(n + 1, 0ll), isBridge(n + 1, 0ll), vis(n + 1, 0ll);
	//dpIN[i] -> min in time node i can reach from going down in the subtree of node i
	//isBridge[node] -> 1 node <-> par[node] is a bridge.

	int timer = 0;
	function<void(int, int)> dfs1 = [&](int node, int p) {
		in[node] = dpIN[node] = timer++;
		par[node] = p;
		vis[node] = 1;
		bool cc = true;
		for (auto ch : g[node]) {
			if (ch == p && cc) {
				cc = false;//multiple edges from node <-> par are possible
				continue;
			}
			if (!vis[ch]) {
				dfs1(ch, node);
			}
			dpIN[node] = min(dpIN[node], dpIN[ch]);
		}

		if (p) {
			isBridge[node] = (dpIN[node] > in[p]);//Bridge only if I go in the subtree and I cant go to lower in time node.
		}
	};

	for (int i = 1; i <= n; i++) {
		if (!vis[i]) {
			dfs1(i, 0);
		}
	}



	for (auto &[a, b] : e) {
		if (par[a] == b) {
			swap(a, b);
		}
		//par,child format
	}
	dsu ds(n + 1);

	for (auto &[p, node] : e) {
		if (par[node] == p && isBridge[node]) {
		} else {
			ds.merge(p, node);//these forms a Non-Bridge Connected componenet.
		}
	}

	vector<vector<int>> grp;
	ds.groups(grp);
	vector<int> id(n + 1, 0ll);// Giving each of these Non-Bridge CC a identity.

	for (int i = 1; i < grp.size(); i++) {
		for (auto x : grp[i]) {
			id[x] = i;
		}
	}
	vector<vector<int>> bridgeTree(grp.size() + 1);
	for (auto &[p, node] : e) {
		if (par[node] == p && isBridge[node]) {
			bridgeTree[id[p]].push_back(id[node]);
			bridgeTree[id[node]].push_back(id[p]);// Creating a bridge from 1 cc to other.
		}
	}
	n = grp.size() - 1;
	//Standard LCA stuff.
	vector<vector<int>> dp(n + 1, vector<int>(LOGN, 0ll));
	vector<int> depth(n + 1, -1ll);// to normalise with 0

	function<void (int, int)> dfs = [&](int node, int par = 0) {
		depth[node] = depth[par] + 1;
		dp[node][0] = par;

		//node-> 2^0 par
		for (int i = 1; i < LOGN; i++) {
			dp[node][i] = dp[dp[node][i - 1]][i - 1];
		}

		for (auto child : bridgeTree[node]) {
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

	for (int i = 1; i <= n; i++) {
		if (depth[i] == -1) {
			dfs(i, 0);
		}
	}

	vector<int> up(n + 1, 0), down(n + 1, 0);
	while (q--) {
		int u, v;
		cin >> u >> v;
		if (!dsg.same(u, v)) {
			cout << "No\n";// Not same cc, I just cant reach here.
			return;
		}
		if (ds.same(u, v)) {
			//ok
			// If same cc Non-bridge CC I can move from any edge to any other. I can form a SCC.
		} else {
			int l = lca(id[u], id[v]);

			// Now sort of prefix sum, standard dp on tree sort of.
			//up from u to lca, down from v to lca, lca pe we cancel the effect.
			// Just some Corner Cases. Likely to miss these.

			if (l == id[u]) {
				down[id[v]]++, down[l]--;
			} else if (l == id[v]) {
				up[id[u]]++, up[l]--;
			} else {
				up[id[u]]++, down[id[v]]++, up[l]--, down[l]--;
			}

		}
	}
	bool ok = 1;
	vis.assign(n, false);
	function<void(int, int)> dfsBT = [&](int node, int par) {
		vis[node] = 1;
		for (auto ch : bridgeTree[node]) {
			if (ch ^ par) {
				dfsBT(ch, node);
				up[node] += up[ch];
				down[node] += down[ch];
			}
		}
		if (par) {
			if (up[node] && down[node]) {
				ok = false;// Both Up and Down orientation is not possible.
			}
		}
	};

	for (int i = 1; i <= n; i++) {
		if (!vis[i]) {
			dfsBT(i, 0);
			if (!ok) {
				cout << "No\n";
				return;
			}
		}
	}

	cout << "Yes\n";

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

