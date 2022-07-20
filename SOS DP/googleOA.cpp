#include <bits/stdc++.h>
using namespace std;

#ifndef ONLINE_JUDGE
#include "algo/debug.h"
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif
/*
Problem -1
We have an array A of N elements. 
We will be asked Q queries,
Each query is in the form of a single integer, X, and we have to tell whether there exists an index i in the array such that the bitwise AND of A[i] and X equals 0. 
If such an index exists print YES, otherwise print NO.

Constraints : 1<=N<=1e5, 1<=Q<=1e5, 1<=A[i]<=1e5
*/
const int LOGN = 20;
const int INF = 1e8;

void solve() {
	int n;
	cin >> n;
	vector<int> a(n);
	vector<int> hsh(1ll << LOGN, 0ll);
	for (auto &x : a) {
		cin >> x;
		hsh[x]++;
	}

	vector<int> val(1ll << LOGN, -INF);
	for (int i = 0; i < 1ll << LOGN; i++) {
		val[i] = (hsh[i] == 0 ? -INF : i);
	}

	for (int i = 1; i <= LOGN; i++) {
		vector<int> nxt(1ll << LOGN, -INF);
		for (int msk = 0; msk < 1ll << LOGN; msk++) {
			if (msk & 1ll << (i - 1)) {
				nxt[msk] = max(val[msk], val[msk ^ (1ll << (i - 1))]);
			} else {
				nxt[msk] = val[msk];
			}
		}
		swap(nxt, val);
	}
	int q;
	cin >> q;

	while (q--) {
		int x;
		cin >> x;

		int ans = val[(((1ll << LOGN) - 1)^x)];
		cout << (ans == -INF ? "NO\n" : "YES\n");
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

