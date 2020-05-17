#include <iostream>
#include <string>
#include <random>
#include <iomanip>
#include "Graph.h";


using namespace std;

class Locality {
private:
    string name;
    size_t population;
public:
    Locality() : name(0), population(0) {}

    Locality(const string& str) : name(str), population(0) {}

    Locality(const string& str, size_t pop) : name(str), population(pop) {}

    Locality(const Locality& obj) = delete;

    Locality& operator=(const Locality& obj) = delete;

    ~Locality() {
        name.clear();
        population = 0;
    }

    Locality& operator=(Locality&& obj) {
        name = obj.name;
        obj.name.clear();
        population = obj.population;
        obj.population = 0;
        return *this;
    }

    string GetName() {
        return name;
    }

    size_t GetPopulation() {
        return population;
    }

    bool operator==(const Locality& rhs) const {
        return name == rhs.name && population == rhs.population;
    }

    bool operator!=(const Locality& rhs) const {
        return !(rhs == *this);
    }
};

void printGraph(const Graph<Locality*, double>& graph) {
    for (auto i = 0; i < graph.GetVertexCount(); ++i) {
        cout << graph.GetVertex(i)->GetName() << endl;
        for (auto j = 0; j < graph.GetEdgeCount(graph.GetVertex(i)); ++j) {
            cout << std::setw(5) << "->" << graph.GetEdge(graph.GetVertex(i), j).dstVertex->GetName() << "(" << graph.GetEdge(graph.GetVertex(i), j).edge << ")" << endl;
        }
        cout << endl;
    }
}

void createRandomEdges(Graph<Locality*, double>* _graph) {
    auto count = _graph->GetVertexCount();

    size_t v1 = 0;
    size_t v2 = 0;
    srand(2);
    for (auto i = 0; i < count; ++i) {
        v1 = rand() % count;
        v2 = rand() % count;
        if (v1 != v2) {
            _graph->AddEdge(_graph->GetVertex(v1), _graph->GetVertex(v2), rand() % 7000);
        }
        else
            --i;
    }
}


int main() {
    Graph<Locality*, double> graph;

    const size_t MAX_POP = 5000000;
    const size_t MAX_WAY = 1500;
    const size_t size = 10;

    auto visited = new bool[size];
    auto stack = new std::size_t[size];

    srand(time(0));

    Locality Samara("Samara", rand() % MAX_POP);
    Locality Moscow("Moscow", rand() % MAX_POP);
    Locality Sochi("Sochi", rand() % MAX_POP);
    Locality Peter("Saint-Petersburg", rand() % MAX_POP);
    Locality Saratov("Saratov", rand() % MAX_POP);
    Locality Bugulma("Bugulma", rand() % MAX_POP);
    Locality Krasnoyarsk("Krasnoyarsk", rand() % MAX_POP);
    Locality Omsk("Omsk", rand() % MAX_POP);
    Locality Rostov("Rostov-on-Don", rand() % MAX_POP);
    Locality Magadan("Magadan", rand() % MAX_POP);

    graph.AddVertex(&Samara);
    graph.AddVertex(&Moscow);
    graph.AddVertex(&Sochi);
    graph.AddVertex(&Peter);
    graph.AddVertex(&Saratov);
    graph.AddVertex(&Bugulma);
    graph.AddVertex(&Krasnoyarsk);
    graph.AddVertex(&Omsk);
    graph.AddVertex(&Rostov);
    graph.AddVertex(&Magadan);

    createRandomEdges(&graph);

    printGraph(graph);

    
    /*cout << "Depth: ";
    DepthFirstSearch(graph, 'a', [](auto vertex) {
        cout << vertex << ' ';
    });
    cout << endl;

    cout << "Breadth: ";
    BreadthFirstSearch(graph, 'a', [](auto vertex) {
        cout << vertex << ' ';
    });
    cout << endl;

    const auto src = 'a';
    const auto dst = 'd';

    {
        char path[11];
        size_t length = 0;
        cout << Dijkstra(graph, src, dst, path, &length) << ": ";
        for (size_t i = 0; i < length; ++i) {
            cout << path[i] << ' ';
        }
        cout << endl;
    }

    {
        char path[11];
        size_t length = 0;
        cout << BellmanFord(graph, src, dst, path, &length) << ": ";
        for (size_t i = 0; i < length; ++i) {
            cout << path[i] << ' ';
        }
        cout << endl;
    }*/
    getchar();
}
