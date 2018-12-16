#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <set>
#include <iterator>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <queue>
using namespace std;

class Timer {
	private:
		time_t startTime, stopTime;
	public:
		double seconds();
		unsigned int milliSeconds();
		void resetTimer();
		void startTimer();
		void stopTimer();
};

struct Node {
	string templateName;
	set<Node*> adjacentTemplates; // a list of structurally_similar_templates
	Node(string t_) : templateName(t_) {};
	Node(string t_, set<Node*> a_) : templateName(t_), adjacentTemplates(a_) {};
	void print();
};

class Graph {
	private:
		set<Node*> g;
		string ppiName;
		vector<string> templatesList;
	public:
		void addNode(string t_);
		void addEdge(string t1_, string t2_);
		set<Node*> getUniverse();
		void print();
		string getPpiName();
		void setPpiName(string name);
		vector<string> getTmplatesList();
		void setTmplatesList(vector<string> tempList);
};

/*
class MaxHeap {
private:
	int *heapArr; // a pointer to the heap array
	int capacity;
	int hSize;
public:
	MaxHeap(vector<pair<string, int> > nameValuePairs);
	void heapify(int i);
	int parent(int i);
	int leftChild(int i);
	int rightChild(int i);
	int getMax();
	void insert(int value);
	int extractMax();
	void updateValue(int i, int newValue);
	void deleteValue(int i);
};
*/

void printCluster(set<Node*> s);
set<Node*> setUnion(set<Node*> a, set<Node*> b);
set<Node*> setIntersection(set<Node*> a, set<Node*> b);
set<Node*> setDifference(set<Node*> a, set<Node*> b);
set<int> setDifference(set<int> a, set<int> b);
vector < set<Node*> > bronKerbosch_v1(set<Node*> R, set<Node*> P, set<Node*> X, vector<set<Node*> > cliques);
vector < set<Node*> > bronKerbosch_v2(set<Node*> R, set<Node*> P, set<Node*> X, vector<set<Node*> > cliques,  ofstream& out);
vector<Graph> parseGraphs(const char *filename);
vector<string> split(string s, string delimiter);
Node* pickPivot(set<Node*> P);
set<string> pickRepresentative(vector<string> templates, vector<set<Node*> > cliques);
void printCompnents(vector<vector<int> > ppiGraph, vector<string> templates, vector<vector<int> > connectedCompnents);
vector<vector<int> > calculateConnectedComponents(vector<vector<int> > adj, vector<string> nodeNames);
vector<Graph> constructSubGraphsFromCC(vector<vector<int> > cc, vector<vector<int> > adj, vector<string> vNames, string gName);
auto compare = [](pair<int, string> *lhs, pair<int, string> *rhs) {
	return (*lhs).first < (*rhs).first;
};

void printQueue(priority_queue<pair<int, string> > q, unordered_map<string, set<int> > templateCliques);

double const THRESHOLD = 0.80;

int main() {
	ifstream in("adjMatrixDiscardZero.txt");
	if (!in) {
		std::cout << "Error opening input file" << endl;
		return 1;
	}
	string line;
	vector<vector<int> > connectedCompnents;
	vector<Graph> subgraphs;
	vector<double> execTime;
	int totalCliquesCount = 0, toalRepCount = 0;
	vector<string> allRepresentatives;

	while (!in.eof()) {
		getline(in, line);
		vector<string> lineEntries;
		lineEntries = split(line, "\t");
		string ppiName;
		vector<string> ppiTemplates;
		vector<vector<int> > ppiGraph;
		if (!lineEntries.empty() && lineEntries[0] == "entry:") {
			ppiName = lineEntries[1];
			std::cout << ppiName << ": " << endl;
			getline(in, line);
			lineEntries = split(line, "\t");
			for (int i = 1; i < lineEntries.size(); i++) {
				ppiTemplates.push_back(lineEntries[i]);
			}
			for (int i = 0; i < ppiTemplates.size(); i++) {
				getline(in, line);
				lineEntries = split(line, "\t");
				vector<int> adj;
				for (int j = 0; j < lineEntries.size(); j++) {
					double similarityScore = atof(lineEntries[j].c_str());
					if (similarityScore >= THRESHOLD) {
						adj.push_back(j);
					}
				}
				ppiGraph.push_back(adj);
			}

			connectedCompnents = calculateConnectedComponents(ppiGraph, ppiTemplates);
			printCompnents(ppiGraph, ppiTemplates, connectedCompnents);
			subgraphs = constructSubGraphsFromCC(connectedCompnents, ppiGraph, ppiTemplates, ppiName);
			/*cout << endl << subgraphs.size() << endl;
			for (auto sg : subgraphs) {
				printCluster(sg.getUniverse());
			}*/
			for (auto sg : subgraphs) {
				ofstream outFile("cliques.txt");
				Timer t;
				set<Node*> R, P, X;
				vector<set<Node*> > cliques;
				//sg.print();
				P = sg.getUniverse();
				t.startTimer();
				cliques = bronKerbosch_v2(R, P, X, cliques, outFile);
				t.stopTimer();
				/*for (set<Node*> cl : cliques) {
					cout << "clique: " << endl;
					printCluster(cl);
				}*/
				outFile.close();
				execTime.push_back(t.milliSeconds());

				t.resetTimer();

				// pick cluster representatives
				vector<string> templates = sg.getTmplatesList();
				std::cout << "Picking Reprentatives for : ";
				printCluster(sg.getUniverse());
				set<string> rep = pickRepresentative(templates, cliques);
				allRepresentatives.insert(allRepresentatives.end(), rep.begin(), rep.end());
				totalCliquesCount += cliques.size();
				toalRepCount += rep.size();
				//std::cout << "picked " << to_string(rep.size()) << " representatives, from " << to_string(cliques.size()) << " cliques." << endl;

			}

			ppiGraph.clear();
			std::cout << endl;
		}
	}
	std::cout << "Average time per graph: " << (accumulate(execTime.begin(), execTime.end(), 0.0) / execTime.size()) / 1000.0 << " seconds." <<endl;
	std::cout << "In total picked " << to_string(toalRepCount) << " representatives, for " << to_string(totalCliquesCount) << " cliques" << endl;
	
	ofstream repFile("representatives.txt");
	for (string r : allRepresentatives) {
		const char* temp = (r + "\n").c_str();
		repFile.write(temp, sizeof(temp));
	}
	repFile.close();
	
	
	return 0;
}


vector < set<Node*> > bronKerbosch_v1(set<Node*> R, set<Node*> P, set<Node*> X, vector<set<Node*> > cliques) {
	if (P.empty() && X.empty()) {
		cliques.push_back(R);
		//cout << "clique found!" << endl;
	}
	set<Node*>::iterator v = P.begin();
	while (!P.empty() && v != P.end()) {
		set<Node*> singleton = { (*v) };
		cliques = bronKerbosch_v1(setUnion(R, singleton), setIntersection(P, (*v)->adjacentTemplates), setIntersection(X, (*v)->adjacentTemplates), cliques);
		P = setDifference(P, singleton);
		X = setUnion(X, singleton);
		if (!P.empty())
			v = P.begin();
	}
	return cliques;
}

vector < set<Node*> > bronKerbosch_v2(set<Node*> R, set<Node*> P, set<Node*> X, vector<set<Node*> > cliques, ofstream& out) {
	if (P.empty() && X.empty()) {
		cliques.push_back(R);
		//out.write("clique,", 7);
	}
	else {
		// pick the pivot as the node with the highest number of neighbors included in P
		Node* u;
		u = pickPivot(setUnion(P,X));
		set<Node*> ext_u = setDifference(P, u->adjacentTemplates);
		//cout << "P-N(u) = ";
		//printCluster(ext_u);

		set<Node*>::iterator v;
		for (v = ext_u.begin(); v != ext_u.end(); ++v) {
			//out.write(((*v)->templateName + ",").c_str(), 7);
			set<Node*> singleton = { (*v) };
			cliques = bronKerbosch_v2(setUnion(R, singleton), setIntersection(P, (*v)->adjacentTemplates), setIntersection(X, (*v)->adjacentTemplates), cliques, out);
			P = setDifference(P, singleton);
			X = setUnion(X, singleton);
			//out.write("back,", 5);
		}
	}
	return cliques;
}

/*
	This function takes a list of connected components, the original graph's adjacency lists and vertex names in
	the original graph. The function returns a list of subgraphs of type (Graph) each representing one of the
	connected components. This list can then be input to bronKerbosch_v1 or bronKerbosch_v2 functions to enumerate
	all maximal cliques
*/
vector<Graph> constructSubGraphsFromCC(vector<vector<int> > cc, vector<vector<int> > adj, vector<string> vNames, string gName) {
	vector<Graph> subgraphs;
	for (vector<int> c : cc) {
		Graph g;
		g.setPpiName(gName);
		vector<string> nodesList;

		for (int v : c) {
			g.addNode(vNames[v]);
			nodesList.push_back(vNames[v]);
		}
		for (int v : c) {
			for (int u : adj[v]) {
				g.addEdge(vNames[v], vNames[u]);
			}
		}
		g.setTmplatesList(nodesList);
		subgraphs.push_back(g);
		//g.print();
	}
	return subgraphs;
}

Node * pickPivot(set<Node*> P)
{
	set<Node*>::iterator v;
	Node* pivot = nullptr;
	int maxNeighbors = -1;
	for (v = P.begin(); v != P.end(); v++) {
		set<Node*> neighbors = (*v)->adjacentTemplates;
		int neighborsCount = (setIntersection(neighbors, P)).size();
		//cout << "P: " << (*v)->templateName << " has " << to_string(neighborsCount) << " neighbors" << endl;
		if (neighborsCount > maxNeighbors) {
			maxNeighbors = neighborsCount;
			pivot = *v;
		}
	}
	//cout << "Pivot " << pivot->templateName << " picked with: " << to_string(maxNeighbors) << " neighbors in P" << endl;
	return pivot;
}

set<string> pickRepresentative(vector<string> templates, vector<set<Node*> > cliques) {
	unordered_map<string, set<int> > templateCliques;
	unordered_map<string, int>  nodesCoversPerTemplate;
	int cliqueIndex = 0;
	set<string> representatives;
	unordered_map<string, Node *> templateNodes;

	for (auto clique : cliques) {
		for (auto temp : clique) {
			int s = cliques[cliqueIndex].size() - 1;
			if (nodesCoversPerTemplate.find(temp->templateName) != nodesCoversPerTemplate.end()) {
				nodesCoversPerTemplate[temp->templateName] += 1;
			} else {
				nodesCoversPerTemplate[temp->templateName] = 1 ;
			}
			
			string key = temp->templateName;
			if (templateCliques.find(key) != templateCliques.end()) {
				templateCliques[key].insert(cliqueIndex);
			} else {
				set<int> cl;
				cl.insert(cliqueIndex);
				templateCliques[key] = cl;
				templateNodes[key] = temp;
			}
		}
		cliqueIndex++;
	}
	priority_queue<pair<int, Node *> > pQueue;
	for (auto t : templates) {
		//pQueue.push(make_pair(cliquesCountPerTemplate[t], t));
		pQueue.push(make_pair(nodesCoversPerTemplate[t], templateNodes[t]));
	}
	set<Node*> nodesCovered;

	while (!pQueue.empty()) {
		set<Node*> singleton = { (pQueue.top().second) };
		if (pQueue.top().second != nullptr && setIntersection(singleton, nodesCovered).empty()) {
			std::cout << "Rep: " << pQueue.top().second->templateName << " For " << to_string(pQueue.top().first) << " clique" << endl;
			representatives.insert(pQueue.top().second->templateName);
			set<int> clustersCovered = templateCliques[pQueue.top().second->templateName];
			for (auto c : clustersCovered) {
				nodesCovered = setUnion(nodesCovered, cliques[c]);
			}
			//printCluster(nodesCovered);
			//cout << "--------------" << endl;
		}
		pQueue.pop();
		
	}
	//printQueue(pQueue, templateCliques);
	
	return representatives;
}

vector<vector<int> > calculateConnectedComponents(vector<vector<int> > adj, vector<string> nodeNames)
{
	std::cout << "Calculating Connected Components ... " << endl;
	// Mark all the vertices as not visited 
	const int V = adj.size();
	vector<bool> visited(V);
	fill(visited.begin(), visited.end(), false);

	queue<int> q;
	vector<vector<int> > connectedComp;
	for (int v = 0; v < V; v++) {
		if (visited[v] == false) {
			vector<int> comp;
			q.push(v);
			//comp.push_back(v);
			visited[v] = true;
			while (!q.empty()) {
				comp.push_back(q.front());
				int t = q.front();
				q.pop();
				for (auto n : adj[t]) {
					if (visited[n] == false) {
						visited[n] = true;
						q.push(n);
					}
				}

			}
			connectedComp.push_back(comp);
			comp.clear();
		}

	}
	return connectedComp;
}

/*
	Util functions
*/
vector<string> split(string s, string delimiter) {
	vector<string> list;
	size_t pos = 0;
	string token;
	while ((pos = s.find(delimiter)) != string::npos) {
		token = s.substr(0, pos);
		list.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	list.push_back(s);
	return list;
}


/*
	Graph methods
*/
void Graph::addNode(string t_) {
	Node* n = new Node(t_);
	g.insert(n);
}

void Graph::addEdge(string t1_, string t2_) {
	Node* n1 = NULL, *n2 = NULL;
	for (auto t : g) {
		if (t->templateName == t1_) { n1 = t; }
		else if (t->templateName == t2_) { n2 = t; }
	}
	if (n1 && n2)
		n1->adjacentTemplates.insert(n2), n2->adjacentTemplates.insert(n1);
}

set<Node*> Graph::getUniverse() {
	return g;
}

string Graph::getPpiName() {
	return ppiName;
}

void Graph::setPpiName(string name) {
	ppiName = name;
}

vector<string> Graph::getTmplatesList() {
	return templatesList;
}

void Graph::setTmplatesList(vector<string> temp_list) {
	templatesList = temp_list;
}


/*
	Set Operations 
*/
set<Node*> setUnion(set<Node*> a, set<Node*> b) {
	set<Node*> c;
	set_union(a.begin(), a.end(), b.begin(), b.end(), inserter(c, c.end()));
	return c;
}

set<Node*> setIntersection(set<Node*> a, set<Node*> b) {
	set<Node*> c;
	set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(c, c.end()));
	return c;
}

set<Node*> setDifference(set<Node*> a, set<Node*> b) {
	set<Node*> c;
	set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(c, c.end()));
	return c;
}

set<int> setDifference(set<int> a, set<int> b) {
	set<int> c;
	set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(c, c.end()));
	return c;
}


/*
	This function was used to parse graphs from the input file and construct them to be input
	directly for bronKerbosch algorithm, without dividing them into their connected components first.
*/
vector<Graph> parseGraphs(const char *filename) {
	vector<Graph> graphs;
	std::cout << "Started parsing adjacency matrices.. " << endl;
	Timer t;
	t.startTimer();
	ifstream infile(filename);
	if (!infile) {
		std::cout << "Error opening input file" << endl;
		return graphs;
	}
	string line;
	vector<string> lineEntries;
	while (!infile.eof()) {
		Graph g;
		getline(infile, line);
		lineEntries = split(line, "\t");
		int numOfNodes = 0;
		if (lineEntries[0] == "entry:") {
			//cout << "Parsing Adjacency Matrix for : " << lineEntries[1] << endl;
			g.setPpiName(lineEntries[1]);

			getline(infile, line);
			lineEntries = split(line, "\t");

			vector<string> nodesList;
			for (int i = 1; i < lineEntries.size(); i++) {
				g.addNode(lineEntries[i]);
				numOfNodes++;
				nodesList.push_back(lineEntries[i]);
			}
			//cout << nodesList.size() << endl;
			g.setTmplatesList(nodesList);
			int i = 1;
			// parse the similarity matrix
			for (int row = 1; row <= numOfNodes; row++) {
				getline(infile, line);
				lineEntries = split(line, "\t");

				for (int column = row + 1; column <= numOfNodes; column++) {
					double similarityScore = ::atof(lineEntries[column - 1].c_str());
					if (similarityScore >= THRESHOLD) {
						g.addEdge(nodesList[row - 1], nodesList[column - 1]);
					}
				}
			}
		}
		graphs.push_back(g);
	}
	t.stopTimer();
	std::cout << "Finished Parsing adjacency matrices in " << to_string(t.seconds()) << " seconds." << endl;
	infile.close();

	return graphs;
}


/**
	For debugging and testing
*/
void printQueue(priority_queue<pair<int, string> > q, unordered_map<string, set<int> > templateCliques) {
	while (!q.empty()) {
		std::cout << q.top().second << " - " << q.top().first << endl;
		for (auto t : templateCliques[q.top().second])
			std::cout << t << ", ";
		std::cout << endl;
		q.pop();
	}
}

void Node::print() {
	std::cout << templateName << " (" << this << ") is structurally similar to : { ";
	for (auto a : adjacentTemplates)
		std::cout << a->templateName << " ";
	std::cout << "}\n";
}

void printCluster(set<Node*> s) {
	std::cout << "[ ";
	for (auto temp : s) {
		std::cout << temp->templateName << " ";
		//for (Node* tt : temp->adjacentTemplates) {
		//	cout << tt->templateName << " ";
		//}
		//cout << endl;
	}
	std::cout << "]\n";
}

void printCompnents(vector<vector<int> > ppiGraph, vector<string> templates, vector<vector<int> > connectedCompnents) {
	sort(connectedCompnents.begin(), connectedCompnents.end(), [](const vector<int> & a, const vector<int> & b) { return a.size() > b.size(); });
	//cout << connectedCompnents.size() << endl;
	for (vector<int> comp : connectedCompnents) {
		std::cout << "templates: ";
		for (int node : comp) {
			std::cout << templates[node] << "\t";
		}

		std::cout << endl;
	}

}

void Graph::print() {
	for (auto t : g)
		(*t).print();
}

void Timer::startTimer() {
	startTime = clock();
}

void Timer::stopTimer() {
	stopTime = clock();
}

void Timer::resetTimer() {
	startTime = stopTime = 0;
}

double Timer::seconds() {
	return double(stopTime - startTime) / CLOCKS_PER_SEC;
}

unsigned int Timer::milliSeconds() {
	return (double(stopTime - startTime) / CLOCKS_PER_SEC) * 1000;
}