#include <bits/stdc++.h>
using namespace std;

#ifndef ONLINE_JUDGE
#include "algo/debug.h"
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif

const int LOGN = 22;
const int INF = 1e8;

void solve() {
	int n;
	cin >> n;
	vector<int> v(n), hsh(1ll << LOGN, 0ll);
	for (auto &x : v) {
		cin >> x;
		hsh[x]++;
	}
	vector<int> val(1ll << LOGN , -INF);
	/*
	val[msk] -> stores an index 'i' such that v[i] is a subset of msk

	*/

	for (int i = 0; i < 1ll << LOGN; i++) {
		val[i] = (hsh[i] == 0 ? -INF : i);
	}
	for (int j = 1; j <= LOGN; j++) {
		vector<int> nxtval(1ll << LOGN, -INF);
		for (int i = 0; i < 1ll << LOGN; i++) {
			if (i & 1ll << (j - 1)) {
				nxtval[i] = max(nxtval[i], max(val[i ^ (1ll << (j - 1))], val[i]));
				//i have the jth but on -> and first j bits are allowed to change , i has an edge from i^1ll<<j and i
			} else {
				nxtval[i] = max(nxtval[i], val[i]);
				//if its switched off cant do anything
			}
		}
		swap(nxtval, val);
	}

	for (int i = 0; i < n; i++) {
		int ans = val[v[i] ^ ((1ll << LOGN) - 1)];
		// subset of v[i]^{all bits on}

		cout << (ans == -INF ? -1 : ans) << " ";
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
