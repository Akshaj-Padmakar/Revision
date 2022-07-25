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
Given A weighted tree G, for each pair of node a,b -> the cost is defined as xr of all values on path frm a to b. % 2 of the on bits in the result.

If some of the edge weights have gone missing, But m new cost are given, determine if the new m gives a plausible tree, print the weights.
*/
void solve() {
	int n, m;
	cin >> n >> m;
	vector<vector<pair<int, int>>> t(n + 1), g(n + 1);
	vector<tuple<int, int, int>> e;
	for (int i = 0, a, b, w; i < n - 1; i++) {
		cin >> a >> b >> w;

		e.push_back(make_tuple(a, b, w));

		if (w != -1) {
			w = __builtin_popcount(w) % 2;
			g[a].push_back({b, w});
			g[b].push_back({a, w});
		}
	}

	for (int i = 0, a, b, w; i < m; i++) {
		cin >> a >> b >> w;
		w = __builtin_popcount(w) % 2;
		g[a].push_back({b, w});
		g[b].push_back({a, w});
	}

	vector<int> val(n + 1, -1ll);
	bool ok = true;
	function<void(int)> dfs = [&](int node) {
		for (auto [ch, w] : g[node]) {
			if (val[ch] == -1) {
				val[ch] = val[node] ^ w;
				dfs(ch);
			} else {
				if (val[ch] != (val[node]^w)) {
					ok = false;
				}
			}
		}
	};
	for (int i = 1; i <= n; i++) {
		if (val[i] == -1) {
			val[i] = 0;
			dfs(i);
			if (!ok) {
				cout << "NO\n";
				return;
			}
		}
	}
	cout << "YES\n";

	for (auto &[a, b, w] : e) {
		if (w == -1) {
			w = val[a]^val[b]);
		}
		cout << a << " " << b << " " << w << "\n";
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
	cin >> t;

	while (t--) {
		solve();
	}

	// cerr << "Time elapsed: " << ((long double)clock() / CLOCKS_PER_SEC) << " s.\n";
}

