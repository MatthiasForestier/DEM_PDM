#pragma once

#include "CRLHelper/VecMatDef.h"

struct HalfEdge {
    int index;      // Index in the halfEdges vector
    int vertex;     // Destination vertex index
    int twin = -1;  // Opposite half-edge index
    int next;       // Next half-edge in the face
    int face;       // Face to which this half-edge belongs
};

struct Vertex {
    int index;                   // Index in the vertices vector
    Vector3F pos;                // Position
    std::vector<int> halfEdges;  // Incident half-edges
};

struct Face {
    int index;     // Index in the faces vector
    int halfEdge;  // One of the half-edges in this face
};

class HalfEdgeMesh {
   public:
    std::vector<Vertex> vertices;
    std::vector<HalfEdge> halfEdges;
    std::vector<Face> faces;

    HalfEdgeMesh() {}
    HalfEdgeMesh(const MatrixXF &V, const MatrixXI &F) { build(V, F); }

    void build(const MatrixXF &V, const MatrixXI &F);

    const Vertex &vertexInFace(int faceIndex, int vertexIndexInFace) const;
    Vector2I edgeVertexIndices(int edgeIndex) const;

    std::vector<int> incidentEdges(int vertexIndex) const { return vertices[vertexIndex].halfEdges; }

    std::vector<int> incidentFaces(int vertexIndex) {
        std::vector<int> faces;
        for (int he : vertices[vertexIndex].halfEdges) {
            faces.push_back(halfEdges[he].face);
        }
        return faces;
    }

    int oppositeFace(int faceIndex, int edgeIndex) {
        int twinEdge = halfEdges[edgeIndex].twin;
        return (twinEdge != -1) ? halfEdges[twinEdge].face : -1;
    }
    int uniqueEdgeIndex(int edgeIndex) {
        int twin = halfEdges[edgeIndex].twin;
        return (twin < edgeIndex) ? twin : edgeIndex;
    }

   private:
    struct PairHash {
        size_t operator()(const std::pair<int, int> &p) const {
            return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
        }
    };
};
