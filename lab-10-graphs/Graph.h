#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <iostream>

// 1)TVertex - some info which store in Vertex
// 2)TEdge - some info which store in edge
// 3)TVertex, TEdge - lightweight
template<typename TVertex, typename TEdge>
class Graph {
public:
    struct Edge {
        TVertex dstVertex;//where should go
        TEdge edge;//some info which store in edge

        Edge(const TVertex &dstVertex, const TEdge edge) : dstVertex(dstVertex), edge(edge) {}
    };

private:
    struct EdgeNode {//node of list
        Edge edge;
        EdgeNode *next;

        EdgeNode(const Edge &edge) : edge(edge), next(nullptr) {}

        EdgeNode(Edge edge, EdgeNode *next) : edge(edge), next(next) {}

    };

    struct Vertex {
        TVertex vertex;//info which store in vertex
        EdgeNode *edges;//list

        Vertex() : edges(nullptr) {}

        Vertex(const TVertex &vertex) : vertex(vertex), edges(nullptr) {}

    };

    Vertex *_graph;
    std::size_t _count;

public:
    Graph() : _graph(nullptr), _count(0) {}

    std::size_t FindVertexIndex(TVertex vertex) const {
        for (std::size_t i = 0; i < _count; ++i)if (_graph[i].vertex == vertex)return i;
        return _count;
    }

    std::size_t FindVertexIndexOrThrow(TVertex vertex) const {
        auto index = FindVertexIndex(vertex);
        if (index == _count)throw std::invalid_argument("Vertex not found!");
        return index;
    }

private:
    //we can store more than one edge between two vertex
    /*static EdgeNode *getEdges(Graph graph) {
        return graph->edges;
    }*/

public:
    ~Graph() {
        clear();
    }

    void clear() {
        /* for (std::size_t i = 0; i < _count; ++i) {
             for (std::size_t j = 0; j < getEdgeCount(_graph[i].vertex); ++j) {
                 this->DeleteEdge(_graph[i].vertex,this->GetEdge(_graph[i].vertex,j).dstVertex,GetEdge(_graph[i].vertex,j).edge);
             }
             _graph[i].edges= nullptr;
         }
         _count = 0;*/
        delete[] _graph;
        _graph = nullptr;
        _count = 0;
    }

    void AddEdge(TVertex srcVertex, TVertex dstVertex, TEdge edge) {
        const auto srcIndex = FindVertexIndexOrThrow(srcVertex);
        //   const auto dstIndex = FindVertexIndexOrThrow(dstVertex);
        Edge e(dstVertex, edge);
        _graph[srcIndex].edges = new EdgeNode(e, _graph[srcIndex].edges);//head insert
    }


    [[nodiscard]]std::size_t GetVertexCount() const {
        return _count;
    }

    [[nodiscard]]TVertex GetVertex(std::size_t index) const {
        return _graph[index].vertex;
    }

    std::size_t getEdgeCount(TVertex vertex) const {
        std::size_t size = 0;
        auto index = FindVertexIndexOrThrow(vertex);
        auto tmp = _graph[index].edges;
        while (tmp) {
            ++size;
            tmp = tmp->next;
        }
        return size;
    }

    Edge GetEdge(TVertex srcVertex, std::size_t index) const {//?
        auto node = _graph[FindVertexIndexOrThrow(srcVertex)].edges;
        for (size_t i = 0; i < index; ++i)
            node = node->next;
        return node->edge;
    }

    /*Edge GetEdge(TVertex srcVertex, std::size_t index) {//?
        auto node = _graph[FindVertexIndexOrThrow(srcVertex)].edges;
        for (size_t i = 0; i < index; ++i)
            node = node->next;
        return node->edge;
    }*/


    void AddVertex(TVertex vertex) {
        const auto index = FindVertexIndex(vertex);
        if (index != _count)
            throw std::invalid_argument("Vertex already exists!");
        auto graph = new Vertex[_count + 1];
        for (std::size_t i = 0; i < _count; ++i)graph[i] = _graph[i];
        graph[index] = vertex;
        delete[] _graph;
        ++_count;
        _graph = graph;
    }

    void DeleteEdge(TVertex srcVertex, TVertex dstVertex) {
        const auto srcIndex = FindVertexIndexOrThrow(srcVertex);
        EdgeNode *edgeNode = _graph[srcIndex].edges;
        if (!edgeNode)throw "edge not found";
        if (edgeNode->edge.dstVertex == dstVertex) {
            _graph[srcIndex].edges = _graph[srcIndex].edges->next;
            delete edgeNode;
            return;
        }
        while (edgeNode->next && edgeNode->next->edge.dstVertex != dstVertex)edgeNode = edgeNode->next;// 8 9 88 99
        if (edgeNode->next == nullptr)throw "edge not found";
        auto tmp = edgeNode->next;
        edgeNode->next = edgeNode->next->next;
        delete tmp;
    }

    void DeleteVertex(TVertex vertex) {
        const auto index = FindVertexIndexOrThrow(vertex);
        for (std::size_t i = 0; i < _count; ++i) {
            if (i != index) {
                for (std::size_t j = 0; j < getEdgeCount(_graph[i].vertex); ++j) {
                    auto edge = this->GetEdge(_graph[i].vertex, j);
                    if (edge.dstVertex == vertex)DeleteEdge(_graph[i].vertex, vertex);
                }
            }
        }
        auto graph = new Vertex[_count - 1];
        for (std::size_t i = 0; i < index; ++i) {
            graph[i] = _graph[i];
        }
        for (std::size_t i = index; i < _count - 1; ++i) {
            graph[i] = _graph[i + 1];
        }
        delete[] _graph;
        --_count;
        _graph = graph;
    }

    void EditEdge(TVertex vertex,TVertex dstVertex, TEdge oldEdge, TEdge newEdge) {
        auto vertexIndex = FindVertexIndexOrThrow(vertex);
        for (std::size_t i = 0; i < this->getEdgeCount(_graph[vertexIndex].vertex); ++i) {
            auto EdgeNode = _graph[vertexIndex].edges;
            if (EdgeNode->edge.edge == oldEdge&&EdgeNode->edge.dstVertex==dstVertex) {
                EdgeNode->edge.edge = newEdge;
                return;
            } else EdgeNode = EdgeNode->next;
        }
        throw "Edge not Find";
    }
    void EditEdge(TVertex vertex, TVertex OldDstVertex,TEdge oldEdge,TVertex newDstVertex, TEdge newEdge) {
        auto vertexIndex = FindVertexIndexOrThrow(vertex);
        auto newDstIndex = FindVertexIndexOrThrow(newDstVertex);
        for (std::size_t i = 0; i < this->getEdgeCount(_graph[vertexIndex].vertex); ++i) {
            auto EdgeNode = _graph[vertexIndex].edges;
            if (EdgeNode->edge.edge == oldEdge&&EdgeNode->edge.dstVertex==OldDstVertex) {
                EdgeNode->edge.edge = newEdge;
                EdgeNode->edge.dstVertex = newDstVertex;
                return;
            } else EdgeNode = EdgeNode->next;
        }
        throw "Edge not Find";
    }

    void EditVertex(TVertex oldVertex, TVertex newVertex) {
        auto vertexIndex = FindVertexIndexOrThrow(oldVertex);
        auto newVertexIndex = FindVertexIndex(newVertex);
        if(newVertexIndex!=_count)throw"vertex already exist";
        *oldVertex = *newVertex;
        return;
    }
/* void EditEdge(TVertex srcVertex, TVertex dstVertex, TVertex newDstVertex,TEdge edge) {
    DeleteEdge(srcVertex, dstVertex);
     AddEdge(srcVertex, newDstVertex, edge);
 }
 void EditEdge(TVertex srcVertex,TVertex dstVertex,TEdge edge){
     DeleteEdge(srcVertex,dstVertex);
     AddEdge(srcVertex,dstVertex,edge);
 }*/

/*void EditVertex(TVertex oldVertex, TVertex newVertex) {
    auto index = FindVertexIndexOrThrow(oldVertex);
    _graph[index].vertex = newVertex;
}*/

/*  bool DelVertex(TVertex vertex) {
      auto vertexIndex = this->FindVertexIndexOrThrow(vertex);
              auto graph = new Vertex[_count - 1];
              for (std::size_t i = 0; i < vertexIndex; ++i) {
                  graph[i] = _graph[i];
              }
              for (std::size_t i = vertexIndex; i < _count - 1; ++i) {
                  graph[i] = _graph[i + 1];
              }
*//*for (auto j = 0; j < _count; ++j) {
DelEdge(_graph[i].vertex, vertex)
}*//*
                delete[] _graph;
                _graph = graph;
                --_count;
                return true;
            }
        return false;
    }
*/

};

template<typename TVertex, typename TEdge, typename TFunctional>
void dfs(const Graph<TVertex, TEdge> &graph, TVertex begin, bool *used, std::size_t *stack, std::size_t &stackSize,
         TFunctional f) {

    auto beginIndex = graph.FindVertexIndexOrThrow(begin);
    used[beginIndex] = true;
    f(graph.GetVertex(beginIndex));
    const TVertex vertex = graph.GetVertex(beginIndex);
    //auto tmp= getEdges(graph[beginIndex]);
    for (std::size_t j = 0; j < graph.getEdgeCount(vertex); ++j) {
        beginIndex = graph.FindVertexIndexOrThrow(graph.GetEdge(vertex, j).dstVertex);
        if (!used[beginIndex]) {
            stack[stackSize++] = beginIndex;
            used[beginIndex] = true;
        }
    }
    while (stackSize > 0) {
        auto vertexIndex = stack[--stackSize];
        dfs(graph, graph.GetVertex(vertexIndex), used, stack, stackSize, f);
    }

}

template<typename TVertex, typename TEdge, typename TFunctional>
void DepthFirstSearch(const Graph<TVertex, TEdge> &graph, TVertex begin, TFunctional f) {
    auto count = graph.GetVertexCount();
    auto stackSize = 0;
    auto used = new bool[count];
    for (std::size_t i = 0; i < count; ++i) {
        used[i] = false;
    }
    auto stack = new std::size_t[count];
    try {
        auto vertexIndex = graph.FindVertexIndexOrThrow(begin);
        used[vertexIndex] = true;
        stack[stackSize++] = vertexIndex;
        while (stackSize > 0) {
            auto vertexInd = stack[--stackSize];
            f(graph.GetVertex(vertexInd));
            auto vertex = graph.GetVertex(vertexInd);
            for (std::size_t i = 0; i < graph.getEdgeCount(vertex); ++i) {
                vertexInd = graph.FindVertexIndexOrThrow(graph.GetEdge(vertex, i).dstVertex);
                if (!used[vertexInd]) {
                    stack[stackSize++] = vertexInd;
                    used[vertexInd] = true;
                }
            }
        }
    }
    catch (...) {
        delete[] used;
        delete[] stack;
        throw;
    }
    delete[] used;
    delete[] stack;
}

template<typename TVertex, typename TEdge, typename TFunctional>
void BreadthFirstSearch(const Graph<TVertex, TEdge> &graph, TVertex begin, TFunctional f) {
    auto count = graph.GetVertexCount();
    auto used = new bool[count];
    auto queueSize = 0;
    for (std::size_t i = 0; i < count; ++i) {
        used[i] = false;
    }
    auto queue = new std::size_t[count];
    try {
        auto beginIndex = graph.FindVertexIndexOrThrow(begin);
        used[beginIndex] = true;
        queue[queueSize++] = beginIndex;
        while (queueSize > 0) {
            auto vertexIndex = queue[0];
            for (std::size_t i = 0; i < queueSize - 1; ++i) {
                queue[i] = queue[i + 1];
            }
            --queueSize;
            auto vertex = graph.GetVertex(vertexIndex);
            f(vertex);
            for (std::size_t i = 0; i < graph.getEdgeCount(vertex); ++i) {
                vertexIndex = graph.FindVertexIndexOrThrow(graph.GetEdge(vertex, i).dstVertex);
                if (!used[vertexIndex]) {
                    queue[queueSize++] = vertexIndex;
                    used[vertexIndex] = true;
                }
            }
        }
    }
    catch (...) {
        delete[] used;
        delete[] queue;
        throw;
    }
    delete[] used;
    delete[] queue;
}

const auto INFINITY_DISTANCE = std::numeric_limits<float>::infinity();

template<typename TVertex, typename TEdge>
double Dijkstra(const Graph<TVertex, TEdge> &graph, TVertex begin, TVertex end, TVertex *path, std::size_t *pathLengh) {
    const auto count = graph.GetVertexCount();
    auto u = new bool[count];//used
    auto d = new TEdge[count];//distance
    auto p = new std::size_t[count];
    for (std::size_t i = 0; i < count; ++i) {
        u[i] = false;
        d[i] = INFINITY_DISTANCE;
        p[i] = count;
    }
    try {
        auto beginIndex = graph.FindVertexIndexOrThrow(begin);
        d[beginIndex] = 0;
        for (std::size_t i = 0; i < count; ++i) {// STEP 1: Select an unmarked vertex with minimum distance
            auto selectDistance = INFINITY_DISTANCE;
            auto selectedIndex = count;
            for (std::size_t j = 0; j < count; ++j) {
                if (!u[j] && d[j] < selectDistance) {
                    selectedIndex = j;
                    selectDistance = d[j];
                }
            }
            if (selectedIndex == count)break;
            u[selectedIndex] = true;
            auto selectedVertex = graph.GetVertex(selectedIndex);
            auto edgeCount = graph.getEdgeCount(selectedVertex);
            for (std::size_t j = 0; j < edgeCount; ++j) {//relaxation
                auto edge = graph.GetEdge(selectedVertex, j);
                auto dstIndex = graph.FindVertexIndexOrThrow(edge.dstVertex);
                if (d[selectedIndex] + edge.edge < d[dstIndex]) {
                    d[dstIndex] = d[selectedIndex] + edge.edge;
                    p[dstIndex] = selectedIndex;
                }
            }
        }
        auto endIndex = graph.FindVertexIndexOrThrow(end);
        double distance = d[endIndex];
        *pathLengh = 0;
        if (distance != INFINITY_DISTANCE && path != nullptr) {
            auto currentIndex = endIndex;
            while (currentIndex != beginIndex) {
                path[(*pathLengh)++] = graph.GetVertex(currentIndex);
                currentIndex = p[currentIndex];
            }
            path[(*pathLengh)++] = begin;
            std::reverse(path, path + *pathLengh);
        }
        delete[] u;
        delete[] p;
        delete[] d;
        return distance;
    }
    catch (...) {
        delete[] u;
        delete[] p;
        delete[] d;
        throw;
    }
}

template<typename TVertex, typename TEdge>
double
BellmanFord(const Graph<TVertex, TEdge> &graph, TVertex begin, TVertex end, TVertex *path, std::size_t *pathLengh) {
    const auto count = graph.GetVertexCount();
    auto d = new double[count];//distance
    auto p = new std::size_t[count];
    for (std::size_t i = 0; i < count; ++i) {
        d[i] = INFINITY_DISTANCE;
        p[i] = count;
    }
    try {
        auto beginIndex = graph.FindVertexIndexOrThrow(begin);
        d[beginIndex] = 0;
        for (std::size_t i = 0; i < count - 1; ++i) {
            for (std::size_t vertexIndex = 0; vertexIndex < count; ++vertexIndex) {
                auto vertex = graph.GetVertex(vertexIndex);
                auto edgeCount = graph.getEdgeCount(vertex);
                for (std::size_t edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
                    auto edge = graph.GetEdge(vertex, edgeIndex);
                    auto dstIndex = graph.FindVertexIndexOrThrow(edge.dstVertex);
                    if (edge.edge + d[vertexIndex] < d[dstIndex]) {
                        d[dstIndex] = edge.edge + d[vertexIndex];
                        p[dstIndex] = vertexIndex;
                    }
                }
            }
        }
        auto endIndex = graph.FindVertexIndexOrThrow(end);
        double distance = d[endIndex];
        *pathLengh = 0;
        if (distance != INFINITY_DISTANCE && path != nullptr) {
            auto currentIndex = endIndex;
            while (currentIndex != beginIndex) {
                path[(*pathLengh)++] = graph.GetVertex(currentIndex);
                currentIndex = p[currentIndex];
            }
            path[(*pathLengh)++] = begin;
            std::reverse(path, path + *pathLengh);
        }
        delete[] p;
        delete[] d;
        return distance;
    }
    catch (...) {
        delete[] d;
        delete[] p;
        throw;
    }

}

