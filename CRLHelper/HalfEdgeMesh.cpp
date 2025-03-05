#include "CRLHelper/HalfEdgeMesh.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>

void HalfEdgeMesh::build(const MatrixXF &V, const MatrixXI &F) {
    int N = V.rows();
    int M = F.rows();

    vertices.resize(N);
    faces.resize(M);

    std::unordered_map<std::pair<int, int>, int, PairHash> edgeMap;

    // Add vertices
    for (int i = 0; i < N; i++) {
        Vector3F pos = V.row(i);
        vertices[i] = {i, pos, {}};
    }

    // Add faces and half-edges
    int edgeCount = 0;
    for (int i = 0; i < M; i++) {
        Vector3I face = F.row(i);
        int v0 = face(0);
        int v1 = face(1);
        int v2 = face(2);

        int e0 = edgeCount++, e1 = edgeCount++, e2 = edgeCount++;
        faces[i] = {i, e0};  // Store index in Face struct

        halfEdges.push_back({e0, v0, -1, e1, i});
        halfEdges.push_back({e1, v1, -1, e2, i});
        halfEdges.push_back({e2, v2, -1, e0, i});

        vertices[v0].halfEdges.push_back(e0);
        vertices[v1].halfEdges.push_back(e1);
        vertices[v2].halfEdges.push_back(e2);

        edgeMap[{v0, v1}] = e0;
        edgeMap[{v1, v2}] = e1;
        edgeMap[{v2, v0}] = e2;
    }

    // Set twin edges
    for (const auto &[key, heIndex] : edgeMap) {
        int vA = key.first, vB = key.second;
        auto twinIt = edgeMap.find({vB, vA});
        if (twinIt != edgeMap.end()) {
            int twinIndex = twinIt->second;
            halfEdges[heIndex].twin = twinIndex;
            halfEdges[twinIndex].twin = heIndex;
        }
    }
}

const Vertex &HalfEdgeMesh::vertexInFace(int faceIndex, int vertexIndexInFace) const {
    int heIndex = faces[faceIndex].halfEdge;
    for (int i = 0; i < vertexIndexInFace; i++) {
        heIndex = halfEdges[heIndex].next;
    }
    return vertices[halfEdges[heIndex].vertex];
}

Vector2I HalfEdgeMesh::edgeVertexIndices(int edgeIndex) const {
    return {halfEdges[edgeIndex].vertex, halfEdges[halfEdges[edgeIndex].twin].vertex};
}
