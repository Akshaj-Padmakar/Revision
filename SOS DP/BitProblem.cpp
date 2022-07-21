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
find the count of y such that, x|y = x -> x should be a subset of x -> 
find the count of y such that, 
*/

const int LOGN = 20;
const int MAXN = (1ll << LOGN) + 5;

void solve() {
	int n;
	cin >> n;
	vector<int> v(n);
	vector<int> hsh(MAXN, 0ll);
	for (auto &x : v) {
		cin >> x;
		hsh[x]++;
	}
 
	vector<int> orr(1ll << LOGN, 0ll); // allow first i bits to vary then find the count.
	for (int i = 0; i < 1ll << LOGN; i++) {
		orr[i] = hsh[i];
	}
	for (int i = 1; i <= LOGN; i++) {
		vector<int> nxt(1ll << LOGN, 0ll); 
 
		for (int msk = 0; msk < 1ll << LOGN; msk++) {
 
			if (msk & (1ll << (i - 1))) {
				nxt[msk] = orr[msk] + orr[msk ^ (1ll << (i - 1))]; // if this bit is on -> we have 2 options to switch this on or not
			} else {
				nxt[msk] = orr[msk]; // just 1 option
			}
		}
		swap(orr, nxt);
	}
 
	vector<int> andd(1ll << LOGN, 0); // allow first i bits to vary, andd[msk] -> vary first i bits, and if its on I have two options
	for (int i = 0; i < 1ll << LOGN; i++) {
		andd[i] = hsh[i];
	}
	for (int i = 1; i <= LOGN; i++) {
		vector<int> nxt(1ll << LOGN, 0ll);
		for (int msk = (1ll << LOGN) - 1; msk >= 0; msk--) {
			if (msk & (1ll << (i - 1))) {
				nxt[msk] = andd[msk]; 
			} else {
				nxt[msk] = andd[msk] + andd[msk ^ (1ll << (i - 1))];
			}
		}
		swap(nxt, andd);
	}
 
 
	for (int i = 0; i < n; i++) {
		cout << orr[v[i]] << " " << andd[v[i]] << " " <<n - orr[((1ll << LOGN) - 1)^v[i]] << "\n"
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
