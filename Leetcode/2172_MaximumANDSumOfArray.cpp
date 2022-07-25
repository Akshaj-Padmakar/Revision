/*
2172. Maximum AND Sum of Array
*/

class Solution {
public:
    int maximumANDSum(vector<int>& nums, int numSlots) {
	while (nums.size() < numSlots * 2) {
		nums.push_back(0);
	}
	int n = nums.size();
	vector<int> dp(1ll << n, 0ll);

	for (int msk = 0; msk < 1ll << n; msk++) {
		int cnt = __builtin_popcount(msk);
		if (cnt % 2) {
			continue;
		}
		for (int i = 0; i < n; i++) {
			if (!(msk & (1ll << i))) {
				for (int j = i + 1; j < n; j++) {
					if (!(msk & (1ll << j))) {
						dp[msk ^ (1ll << i) ^ (1ll << j)] = max(dp[msk ^ (1ll << i) ^ (1ll << j)], dp[msk] + ((nums[i] & (cnt / 2 + 1)) + (nums[j] & (cnt / 2 + 1))));
					}
				}
			}
		}
	}

	return dp[(1ll << n) - 1];
}
};
