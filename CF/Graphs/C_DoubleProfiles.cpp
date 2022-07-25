#include <bits/stdc++.h>
using namespace std;

#ifndef ONLINE_JUDGE
#include "algo/debug.h"
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif

#define int long long

const int MOD = 1e9 + 7;
const int FMOD = 998244353;
int add(int a, int b, int M = MOD) {
	return ((a + b) % M + M) % M;
}
int sub(int a, int b, int M = MOD) {
	return ((a - b) % M + M) % M;
}
/*
Given a Graph G, count the number of pairs i,j such that either i,j have a edge and rest all the nodes in the adj list are same, o.w i,j dont have a edge and all nodes are same in its adj list.
The benift of using two hash ?

So, If I use two hash we have 1/M1*M2 probablity of collision for a single comparsion, and using M1 ~ M2 ~ 10 ^9 decreases collision probablity by a lot.

https://codeforces.com/contest/154/problem/C
*/
void solve() {
	int n, m;
	cin >> n >> m;
	vector<vector<int>> g(n + 1);
	vector<pair<int, int>> e;
	for (int i = 0, a, b; i < m; i++) {
		cin >> a >> b;
		g[a].push_back(b);
		g[b].push_back(a);
		e.push_back({a, b});
	}

	vector<pair<int, int>> val(n + 1);
	map<pair<int, int>, int> dp;

	vector<int> pw(n + 1, 0ll), pw1(n + 1, 0ll);
	pw[0] = 1, pw1[0] = 1;

	for (int i = 1; i <= n; i++) {
		pw[i] = pw[i - 1] * 31; pw[i] %= MOD;
		pw1[i] = pw1[i - 1] * 31; pw1[i] %= FMOD;
	}

	for (int i = 1; i <= n; i++) {
		int sm = 0, sm1 = 0;
		for (auto x : g[i]) {
			sm += pw[x];
			sm1 += pw1[x];
			sm %= MOD, sm1 %= FMOD;
		}

		val[i] = {sm, sm1};
		dp[ {sm, sm1}]++;
	}
	int ans = 0;
	for (int i = 1; i <= n; i++) {
		ans += dp[val[i]] - 1;
	}
	ans /= 2;

	for (auto [a, b] : e) {
		if (sub(val[a].first, pw[b], MOD) == sub(val[b].first, pw[a], MOD) && sub(val[a].second, pw1[b], FMOD) == sub(val[b].second, pw1[a], FMOD)) {
			ans++;
		}
	}

	cout << ans << "\n";
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

