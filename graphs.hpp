/*
    LARGE STACK SIZE --> !!
    g++ -std=c++17 -Wl,-stack_size,0x10000000 A.cpp -o A
    ./A <input.txt >output.txt
*/
auto check = [&](int a) {
	// when returning a value
};

function<void (int, int)> dfs = [&](int node, int p) {
	// for writing a recursion
};

1
// ))(???)(?(
__builtin_popcount
// -------------------------------<ADJACENCY_LIST>-------------------------------

//1. > Graphs representation


int node, edges;

cin >> nodes >> edges;

rep(i, edges) {
	int a, b;
	cin >> a >> b;
	// for undirected graphs;
	v[a].pb(b);
	v[b].pb(a);
}

// -------------------------------<ADJACENCY_LIST>-------------------------------



// -----------------------------------<DFS>----------------------------------

vi vis(MAXN);
vi v[MAXN];

void dfs(int node) {
	vis[node] = 1;


	// for-each loop;
	for (int child : v[node]) {
		if (vis[child] == 0) {
			dfs(child);
		}
	}
}

//time complexity-- > O(node + edges);

// -----------------------------------</DFS>----------------------------------


// calucalting connected components in a graph;

// strongly connencted components--> only single path for directed graphs

// weakly connected components--> making directed paths undirected;


int main() {
	int nodes, edges;
	cin >> nodes >> edges;

	int cc = 0;
	for (int i = 1; i <= nodes; i++) {
		if (vis[nodes] == 0) {
			dfs(nodes);
			cc++;
		}
	}
	cout << cc << "\n";
}



// -----------------------------------<SSSP>----------------------------------

//Single Sources Shortest Path (on tree)
// shortest distance between a source node and everyother node in the graph

vi v[MAXN];
vi vis(MAXN);
vi dis(MAXN);

void dfs(int node, int par int distance) {

	//vis[node] = 1;
	dis[node] = distance;
	for (int child : v[node]) {
		if (child != par) {
			dfs(child, node , distance + 1);
		}
	}
}

function<void(int, int, int)> dfs = [&](int node, int p, int d) {
	dis[node] = d;
	for (auto ch : g[node]) {
		if (ch ^ p) {
			dfs(ch, node, d + 1);
		}
	}
}

// -----------------------------------</SSSP>----------------------------------


// -----------------------------------<Bipartite_Graph>----------------------------------

// if there is a edge then the color of the vertices joining this edge should be diff
//each edge should have a color

vi v[MAXN];
vi vis(MAXN);
vector<bool> col(MAXN);
bool dfs(int node, int c) {
	vis[node] = 1;

	col[node] = c;

	for (int child : v[node]) {
		if (vis[child] == 0) {
			if (!dfs(child, c ^ 1)) {
				return false;
			}
		} else if (col[child] == col[node]) {
			return false;
		}
	}

	return true;
}

function<bool(int, int)> dfs = [&](int node, int c) {
	vis[node] = 1;
	col[node] = c;
	for (auto ch : g[node]) {
		if (!vis[ch]) {
			if (!dfs(ch, c ^ 1)) {
				return false;
			}
		} else {
			if (col[node] == col[ch]) {
				return false;
			}
		}
	}
	return true;
};


// -----------------------------------</Bipartite_Graph>----------------------------------



// -----------------------------------<Cycle_detection>----------------------------------

//any edge that connects it to its ancestor which is not its parent is called a "BACK-EDGE"
//BACK-EDGE indicates a cycle in the graph

// For Undirected Graph !!
vi v[MAXN];
vi vis(MAXN);

bool dfs(int node, int par) {
	vis[node] = 1;

	for (int child : v[node]) {
		if (vis[child] == 0) {
			if (dfs(child, node) == true) {
				return true;
			}
		} else {
			if (child != par) {
				return true;
			}
		}
	}
	return false;
}

function<bool(int, int)> dfs = [&](int node, int p) {
	vis[node] = 1;
	for (auto ch : g[node]) {
		if (vis[ch]) {
			if (ch != p) {
				return true;
			}
		} else {
			if (dfs(ch, node)) {
				return true;
			}
		}
	}
	return false;
}
// returns true for a cycle





vi v[MAXN];
vi vis(MAXN);
vi par(MAXN);
bool dfs(int node, int p) {
	vis[node] = 1;
	par[node] = par;
	for (auto child : v[node]) {
		if (!vis[child]) {
			if (dfs(child, node)) {
				return 1;
			}
		} else if (vis[node] == 1) {
			s = node;
			e = child;
			// x= s---> e
			return 1;
		}
	}
	return 0;
}


// For A Directed Graph


int n;
vector<vector<int>> adj;
vector<char> color;
vector<int> parent;
int cycle_start, cycle_end;

bool dfs(int v) {
	color[v] = 1;
	for (int u : adj[v]) {
		if (color[u] == 0) {
			parent[u] = v;
			if (dfs(u))
				return true;
		} else if (color[u] == 1) {
			cycle_end = v;
			cycle_start = u;
			return true;
		}
	}
	color[v] = 2;
	return false;
}

void find_cycle() {
	color.assign(n, 0);
	parent.assign(n, -1);
	cycle_start = -1;

	for (int v = 0; v < n; v++) {
		if (color[v] == 0 && dfs(v))
			break;
	}

	if (cycle_start == -1) {
		cout << "Acyclic" << endl;
	} else {
		vector<int> cycle;
		cycle.push_back(cycle_start);
		for (int v = cycle_end; v != cycle_start; v = parent[v])
			cycle.push_back(v);
		cycle.push_back(cycle_start);
		reverse(cycle.begin(), cycle.end());

		cout << "Cycle found: ";
		for (int v : cycle)
			cout << v << " ";
		cout << endl;
	}
}

// -----------------------------------</Cycle_detection>----------------------------------


// -----------------------------------<In_Out_Time>----------------------------------

//
int timer = 1;
vi in[(int)1e5 + 5];
vi out[(int)1e5 + 5];

void dfs(int node) {
	vis[node] = 1;

	in[v] = timer++;
	// timer increased after saving in array.

	for (int child : v[node]) {
		if (vis[child] == 0) {
			dfs(child);
		}
	}

	out[v] = timer++;
}


function<void(int)> dfs = [&](int node) {
	vis[node] = 1;
	in[node] = timer++;

	for (auto ch : g[node]) {
		if (!vis[ch]) {
			dfs(ch);
		}
	}
	out[node] = timer++;
};

function<bool(int, int)> isAncestor = [&](int ancestor, int child) {
	return in[ancestor] <= in[child] && out[ancestor] >= out[child];
};

// if a node(x) lies in a subtree of another node(y) then the in[x] > in[y] && out[x]<out[y]

// 	otherwise in[x]>in[y] && out[x] > out[y] then x , y doesnt lie in the subtree of other node

// -----------------------------------</In_Out_Time>----------------------------------




// -----------------------------------<Diameter of tree>----------------------------------

// spoj --> longest path in a tree

void dfs(int n, int d) {
	vis[node] = 1;
	if (d > max_d) {
		max_d = d;
		max_node = node;
	}

	for (int child : v[node]) {
		if (vis[child] == 0) {
			dfs(child, d + 1);
		}
	}
}

void solvethetestcase() {
	// reading the tree;

	max_node = -1;

	dfs(1, 0);
	// started frm 1 with distnace ? 1

	rep(i, n + 3) {
		vis[i] = 0;
	}

	max_d = -1;
	dfs(max_node, 0);
	cout << max_d;
}

function<int(int, int)> diameter = [&]() {
	SSSP_Tree(1, 0, 0);
	int mx_node = max_element(dis.begin(), dis.end()) - dis.begin();

	dfs(mx_node, 0, 0);

	return max_element(dis.begin, dis.end());
};

// -----------------------------------</Diameter of tree>----------------------------------


// -----------------------------------<CENTROID>----------------------------------


vi v[MAXN];
vi sz(MAXN);
int dfs(int node, int p) {
	sz[node] = 1;
	for (auto child : v[node]) {
		if (child ^ p) {
			sz[node] += dfs(child, node);
		}
	}
	return sz[node];
}
int n;
int dfs_centroid(int node, int p) {
	for (auto child : v[node]) {
		if (child ^ p) {
			if (sz[child] * 2 > n) {
				return dfs_centroid(child, node);
			}
		}
	}
	return node;
}
void solvethetestcase() {
	cin >> n;

	rep(i, n - 1) {
		int a, b;
		cin >> a >> b;
		v[a].pb(b);
		v[b].pb(a);
	}

	dfs(1, -1);
	cout << dfs_centroid(1, -1) << "\n";
}





function<vector<int>()> Centroid = [&]() {
	vector<int> centroid;
	vector<int> sz(n + 1);
	function<void (int, int)> dfs = [&](int u, int par) {
		sz[u] = 1;
		bool is_centroid = true;
		for (auto v : v[u])
			if (v ^ par) {
				dfs(v, u);
				sz[u] += sz[v];
				if (sz[v] > n / 2) {
					is_centroid = false;
				}
			}
		if (n - sz[u] > n / 2) is_centroid = false;
		if (is_centroid) centroid.push_back(u);
	};
	dfs(1, -1);

	return centroid;
};

// -----------------------------------</CENTROID>----------------------------------


// -----------------------------------<Subtree_Size_Of_Tree>----------------------------------

vi subsize[N];

// stores the subtree size of node i
// dfs returns the subtree size of node

int dfs(int node) {
	v[node] = 1;
	int curr_size = 1;

	for (int child : v[node]) {
		if (!vis[child]) {
			curr_size += dfs(child);
		}
	}

	subsize[node] = curr_size;

	return curr_size;
}

function<void(int, int)> dfs = [&](int node, int par) {
	sz[node] = 1;
	for (auto ch : g[node]) {
		if (ch ^ par) {
			dfs(ch, node);
			sz[node] += sz[ch];
		}
	}
};

// -----------------------------------</Subtree_Size_Of_Tree>----------------------------------



// -----------------------------------------<BFS>-------------------------------------

//https://atcoder.jp/contests/abc168/tasks/abc168_d

vi v[MAXN];
vi dis(MAXN, INF);
vi vis(MAXN, 0ll);
vi par(MAXN);
// stores where I came from.
void bfs(int source) {
	queue<int> q;
	q.push(1);
	dis[1] = 0; vis[1] = 1; par[1] = -1;

	while (q.size()) {
		int x = q.front();
		q.pop();
		for (auto child : v[x]) {
			if (!vis[child]) {
				dis[child] = dis[x] + 1;
				vis[child] = 1;
				q.push(child);
				par[child] = x;
			}
		}
	}
}



function<void(int)> bfs = [&](int source) {
	queue<int> q;
	q.push(source);
	dis[source] = 0;

	while (q.size()) {
		auto node = q.front();
		q.pop();

		for (auto ch : g[node]) {
			if (dis[ch] > dis[node] + 1) {
				dis[ch] = dis[node] + 1;
				q.push(ch);
			}
		}
	}
};

queue<int> q;

//source
q.push(path.back());
vis[path.back()] = 1;
dis[path.back()] = 0;
while (!q.empty()) {
	int node = q.front();
	q.pop();

	for (auto child : g2[node]) {
		if (!vis[child]) {
			q.push(child);
			vis[child] = 1;
		}
		//dis[child] = (dis[child] > dis[node] + 1) ? (dis[node] + 1 ) : dis[child];
		if (dis[child] > dis[node] + 1) {
			dis[child] = dis[node] + 1;

		}
		debug(dis[child], child, node, dis[node]);
	}
}

// -----------------------------------------</BFS>-------------------------------------


// -----------------------------------------<BRIDGE>-------------------------------------

// A BRIDGE is a edge which when removed increses the connceted component of the graph by one

// A EDGE can be -->

// BACK EDGE - connects a node to its ancestor which is not its parent

// FORWARD EDGE - edges which are traversed during a dfs call !
//                forms the basis for DFS tree of a graph


vi v[MAXN];
vi in(MAXN);
vi low(MAXN);
vi vis (MAXN);
int timer;

// in time  --> time you enter the node in dfs
// low time --> least time taken by a node's ancestors is low time


void dfs(int node, int parent) {
	vis[node] = 1;
	in[node] = low[node] = timer;

	timer++;

	for (auto child : v[node]) {
		if (child == par) {
			continue;
		}

		if (vis[child]) {
			// this is a back edge {child - node}

			// minimize the low time of node

			low[node] = min(low[node], in[child]);

		}

		if (vis[child] == 0) {
			// this denotes a forward edge

			dfs(child, node);

			if (low[child] > in[node]) {
				// this a bridge
				// child - node

				// why ??
				// bc otherwise node is connected to other node whoese in time is even smaller
				//than the in time of child
				// hence cann't be a bridge
			}


			low[node] = min(low[node], low[child]);
		}


	}
}


// -----------------------------------------</BRIDGE>-------------------------------------


// -------------------------------------<ARTICULATION_POINT>-------------------------------

//nodes which when removed increase the cc of the graph
// has no relation with bridges

vi v[MAXN];
vi vis(MAXN);
vi low(MAXN);
vi in(MAXN);

int timer;

void dfs(int node, int parent = -1) {
	vis[node] = 1;
	in[node] = low[node] = timer;
	timer++;
	int children = 0;
	for (auto child : v[node]) {
		if (child == parent) {
			continue;
		} else if (vis[node] == 1) {
			//backedge
			low[node] = min(low[node], in[child]);

		} else {
			dfs(child, node);
			low[child] = min(low[child], low[node]);
			if (low[node] > in[child] && p != -1) {
				// node is a articulation point
			}
			children++;

		}

	}
	if (p == -1 && children > 1) {
		// root is a articulation point
	}

}




// -------------------------------------</ARTICULATION_POINT>-------------------------------



// -----------------------------------------<KAHN'S ALGO>-----------------------------------
// for directed acyclic graph (DAG)
vi v[MAXN];
vi indegree(MAXN);
vi ans;



void kahn(int n) {
	queue<int> q;
	ford(i, 1, n + 1) {
		if (indegree[i] == 0) {
			q.push(i);
		}
	}
	while (!q.empty()) {
		int x = q.front();
		q.pop();
		ans.pb(x);
		for (auto child : v[x]) {
			indegree[child]--;
			if (indegree[child] == 0) {
				q.push(child);
			}
		}
	}
}

function<void(int)> kahn = [&]() {
	queue<int> q;
	for (int i = 1; i <= n; i++) {
		if (in[i] == 0) {
			q.push(i);
		}
	}

	while (q.size()) {
		int node = q.front();
		q.pop();
		tpsrt.push_back(node);

		for (auto ch : g[node]) {
			in[ch]--;
			if (!in[ch]) {
				q.push(ch);
			}
		}
	}
	assert(ans.size() == n);
	//Check for a DAG
}


void solvethetestcase() {
	int n, m;
	cin >> n >> m;

	rep(i, m) {
		int a, b;
		cin >> a >> b;
		v[a].pb(b);
		indegree[b]++;
		//a->b
	}

	kahn(int n);

	// ans stores our topologically sorted nodes

}





int n; // number of vertices
vector<vector<int>> adj; // adjacency list of graph
vector<bool> visited;
vector<int> ans;

void dfs(int v) {
	visited[v] = true;
	for (int u : adj[v]) {
		if (!visited[u])
			dfs(u);
	}
	ans.push_back(v);
}

void topological_sort() {
	visited.assign(n, false);
	ans.clear();
	for (int i = 0; i < n; ++i) {
		if (!visited[i])
			dfs(i);
	}
	reverse(ans.begin(), ans.end());
}





// -----------------------------------------</KAHN'S ALGO>-----------------------------------





// -----------------------------------------<DIJSKTRA'S ALGORITHM>-----------------------------------
// min weighted path from a source

vpii v[MAXN];
vi dis(MAXN, INF);


void solvethetestcase() {
	int n, m;
	cin >> n >> m;
	rep(i, m) {
		int a, b, w; cin >> a >> b >> w;
		v[a].pb(mp(b, w));
		v[b].pb(mp(a, w));
	}

	PQD(pii) p;
	// dist, node
	p.push(mp(0, 1));
	dis[1] = 0;

	while (!p.empty()) {
		auto x = p.top();
		p.pop();

		if (x.S > dis[x])continue;

		for (auto xyz : v[x]) {
			if (x.F + xyz.S < dis[xyz.F]) {
				p.push(mp(xyz.S, xyz.F));
				dis[xyz.F] = x.F + xyz.S;
			}
		}
	}
}
}
//0 1
//





vector<pair<int, int>> dr1 = {{ -1, 0}, {1, 0}, {0, 1}, {0, -1}};
vector<pair<int, int>> dr2 = {{1, 0}, { -1, 0}, {0, -1}, {0, 1}, {2, 0}, { -2, 0}, {0, 2}, {0, -2}, {1, 2}, {2, 1}, { -1, 2}, {2, -1}, {1, -2}, { -2, 1}, { -1, -2}, { -2, -1}};




function<void(int)> dijktra = [&](int source) {
	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > pq;
	//dis,node

	pq.push({source, 0});


	while (pq.size()) {
		auto [d, node] = pq.top();
		pq.front();

		for (auto [ch, w] : g[node]) {
			if (dis[ch] > dis[node] + w) {
				dis[ch] = dis[node] + w;
				pq.push(dis[ch], ch);
			}
		}
	}
};



// -----------------------------------------</DIJSKTRA'S ALGORITHM>-----------------------------------





// -----------------------------------------<MST>-----------------------------------

// -----------------------------------------<Kruskal's Algo.>-----------------------------------
// A greedy Algo
// can be proved by contradiction
void solvethetestcase() {
	int n, m;
	cin >> n >> m;

	dsu ds(n);

	vector<pair<int, pii>> G;
	rep(i, m) {
		int a, b, w;
		cin >> a >> b >> w;
		a--, b--;
		G.pb({w, {a, b}});
	}

	sort(all(G));
	// on the basis of weight

	int ans = 0;
	rep(i, m) {
		int a = G[i].S.F;
		int b = G[i].S.S;
		int w = G[i].F;

		if (ds.same(a, b)) {
			//ok

		} else {
			ans += w;
			ds.merge(a, b);
		}
	}
	if (ds.groups().size() > 1) {
		cout << -1 << "\n";
		return;
		// we had dis-connected graph intially!
	}

	cout << ans << "\n";
}





// -----------------------------------------</Kruskal's Algo.>-----------------------------------



// -----------------------------------------<LCA>-----------------------------------


vi v[MAXN];
int dp[MAXN][LOGN];
vi depth(MAXN);

void dfs(int node, int par = 0) {
	depth[node] = depth[par] + 1;
	dp[node][0] = par;

	//node-> 2^0 par
	for (int i = 1; i < LOGN; i++) {
		dp[node][i] = dp[dp[node][i - 1]][i - 1];
	}

	for (auto child : v[node]) {
		if (child == par)continue;

		dfs(child, node);
	}
}


int ancestor(int node, int k) {
	//par k distnace up !

	for (int i = 0; i < LOGN && node; i++) {
		if (k & (1ll << i)) {
			node = dp[node][i];
		}
	}
	return (node == 0 ? -1 : node);
}


int lca(int node1, int node2) {
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
}





//---------------------<Inside main!>-----------------------------

vector<vector<int>> dp(n + 1, vector<int>(LOGN, 0ll));
vector<int> depth(n + 1, -1ll);// to normalise with 0

function<void (int, int)> dfs = [&](int node, int par = 0) {
	depth[node] = depth[par] + 1;
	dp[node][0] = par;

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
function<int (int, int)> jump = [&](int node, int val) {
	//jumps from node to node x such that age[x]<=val but age[par[x]]>val
	for (int i = LOGN - 1; ~i; --i) {
		if (dp[node][i] && age[dp[node][i]] <= val) {
			node = dp[node][i];
		}
	}
	return node;
};
// -----------------------------------------</LCA>-----------------------------------



















































