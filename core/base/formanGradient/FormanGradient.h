/// \ingroup base
/// \class ttk::FormangGradient
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %formangGradient processing package.
///
/// %FormangGradient is a TTK processing package that takes a scalar field on
/// the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkFormangGradient.cpp %for a usage example.

#pragma once

// ttk common includes
#include <Debug.h>
#include <MemoryUsage.h>
#include <Triangulation.h>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/optional/optional.hpp>

using namespace std;
using namespace ttk;

typedef pair<short int, SimplexId> Simplex;
typedef set<Simplex, boost::function<bool(const Simplex &, const Simplex &)>>
  SimplexesSet;

namespace ttk {

  class FormanGradient : virtual public Debug {

  public:
    FormanGradient();

    ~FormanGradient();

    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeGradient(dataType *outputData,
                        const dataType *inputData,
                        const triangulationType *triangulation) {

      double maxF = 0;
      SimplexId vertexNum = triangulation->getNumberOfVertices();

      vector<pair<dataType, SimplexId>> thepairs(vertexNum);
      for(int i = 0; i < vertexNum; i++) {
        maxF = std::max(maxF, (double)inputData[i]);
        thepairs[i] = pair<dataType, SimplexId>(inputData[i], i);
      }

      sort(thepairs.begin(), thepairs.end());

      // building the indexing
      indexing_ = vector<SimplexId>(vertexNum);
      for(SimplexId i = 0; i < vertexNum; i++) {
        indexing_[thepairs[i].second] = i;
      }

      // prepare gradient
      gradient_.push_back(vector<vector<SimplexId>>(
        triangulation->getNumberOfVertices(), vector<SimplexId>(1, -1)));
      gradient_.push_back(vector<vector<SimplexId>>(
        triangulation->getNumberOfEdges(), vector<SimplexId>(2, -1)));

      if(dimensionality_ == 2) {
        gradient_.push_back(vector<vector<SimplexId>>(
          triangulation->getNumberOfTriangles(), vector<SimplexId>(1, -1)));
      } else {
        gradient_.push_back(vector<vector<SimplexId>>(
          triangulation->getNumberOfTriangles(), vector<SimplexId>(2, -1)));
        gradient_.push_back(vector<vector<SimplexId>>(
          triangulation->getNumberOfCells(), vector<SimplexId>(1, -1)));
      }

      // compute indexing on the data
      SimplexId vertexNumber = triangulation->getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < vertexNumber; i++) {
        // compute homotopy expansion
        // logFile_ << "i" << std::endl;
        // logFile_ << "-----------------------------------------" << std::endl;
        homotopyExpansion(i, triangulation);
      }

      return 1; // success
    }

    int preconditionTriangulation(ttk::AbstractTriangulation *triangulation) {
      dimensionality_ = triangulation->getDimensionality();

      if(triangulation) {
        triangulation->preconditionVertexStars();
        triangulation->preconditionVertexEdges();
        triangulation->preconditionEdges();

        if(dimensionality_ >= 2) {
          triangulation->preconditionTriangles();
          triangulation->preconditionEdgeTriangles();
          triangulation->preconditionTriangleEdges();
        }

        if(dimensionality_ == 3) {
          triangulation->preconditionCellTriangles();
          triangulation->preconditionTriangleStars();
          triangulation->preconditionVertexTriangles();
        }
      }
      return 0;
    };

    template <typename triangulationType>
    int extractBoundary(Simplex simpl,
                        int dimension_b,
                        int index_b,
                        const triangulationType &triangulation) {
      SimplexId id = -1;
      switch(simpl.first) {
        case 0:
          return id;
        case 1:
          triangulation->getEdgeVertex(simpl.second, index_b, id);
          break;
        case 2:
          if(dimension_b == 0)
            triangulation->getTriangleVertex(simpl.second, index_b, id);
          else {
            triangulation->getTriangleEdge(simpl.second, index_b, id);
          }
          break;
        case 3:
          if(dimension_b == 0)
            triangulation->getCellVertex(simpl.second, index_b, id);
          else if(dimension_b == 0)
            triangulation->getCellEdge(simpl.second, index_b, id);
          else
            triangulation->getCellTriangle(simpl.second, index_b, id);
          break;
      }

      return id;
    }

    template <typename triangulationType>
    int extractCoboundary(Simplex simpl,
                          int dimension_cb,
                          int index_cb,
                          const triangulationType &triangulation) {
      SimplexId id = -1;
      // logFile_ << "dim: " << simpl.first << ", id: " << simpl.second << std::endl;

      switch(simpl.first) {
        case 0:
          if(dimension_cb == 1)
            triangulation->getVertexEdge(simpl.second, index_cb, id);
          else if(dimension_cb == 2)
            triangulation->getVertexTriangle(simpl.second, index_cb, id);
          else
            triangulation->getVertexStar(simpl.second, index_cb, id);
          break;

        case 1:
          if(dimension_cb == 2)
            triangulation->getEdgeTriangle(simpl.second, index_cb, id);
          else
            triangulation->getEdgeStar(simpl.second, index_cb, id);
          break;

        case 2:
          triangulation->getTriangleStar(simpl.second, index_cb, id);
          break;
        case 3:
          break;
      }

      return id;
    }

    template <typename triangulationType>
    int sizeCofacets(Simplex simpl,
                     int dimension_cb,
                     const triangulationType &triangulation) {
      switch(simpl.first) {
        case 0: {
          if(dimension_cb == 1)
            return triangulation->getVertexEdgeNumber(simpl.second);
          else if(dimension_cb == 2)
            return triangulation->getVertexTriangleNumber(simpl.second);
          else
            return triangulation->getVertexStarNumber(simpl.second);
        }

        case 1: {
          if(dimension_cb == 2)
            return triangulation->getEdgeTriangleNumber(simpl.second);
          else
            return triangulation->getEdgeStarNumber(simpl.second);
        }

        case 2:
          if(dimensionality_ > 2)
            return triangulation->getTriangleStarNumber(simpl.second);
          else
            return 0;
        case 3:
          return 0;
      }

      return 0;
    }

    template <typename triangulationType>
    void simplexToVertices(Simplex simplex,
                           vector<SimplexId> &vertices,
                           const triangulationType &triangulation) {
      vertices = vector<SimplexId>(simplex.first + 1);
      // Timer t;
      for(int i = 0; i < simplex.first + 1; i++) {
        SimplexId vertex;

        // t.reStart();
        switch(simplex.first) {
          case 0:
            vertex = simplex.second;
            break;
          case 1:
            triangulation->getEdgeVertex(simplex.second, i, vertex);
            break;
          case 2:
            triangulation->getTriangleVertex(simplex.second, i, vertex);
            break;
          case 3:
            triangulation->getCellVertex(simplex.second, i, vertex);
            break;
        }
        // relationTime_ += t.getElapsedTime();

        vertices[i] = vertex;
      }
    }

  protected:
    template <typename triangulationType>
    void homotopyExpansion(int v, const triangulationType &triangulation) {
      auto lexycographic
        = bind(&FormanGradient::cmpSimplexesLevelSets<triangulationType>, this,
               _1, _2, triangulation);
      SimplexesSet pqone = SimplexesSet(lexycographic);
      SimplexesSet pqzero = SimplexesSet(lexycographic);
      SimplexesSet critical = SimplexesSet(lexycographic);

      Simplex vertex = Simplex(0, v);

      pqzero.insert(vertex);
      // Timer t;
      int szCofacets = sizeCofacets(vertex, 1, triangulation);
      // relationTime_ += t.getElapsedTime();
      // logFile_ << "szCofacets: " << szCofacets << std::endl;

      for(int i = 0; i < szCofacets; i++) {
        // t.reStart();
        int coface = extractCoboundary(vertex, 1, i, triangulation);
        // relationTime_ += t.getElapsedTime();
        int index = getIndex(Simplex(1, coface), triangulation);
        if(index == indexing_[v])
          pqone.insert(Simplex(1, coface));
      }

      while(!pqzero.empty()) {
        while(!pqone.empty()) {
          Simplex top = *pqone.begin();
          pqone.erase(top);

          Simplex pairable;
          int numPairable
            = numPairableFaces(top, critical, pairable, triangulation);

          if(numPairable == 0) {
            pqzero.insert(top);
          } else {
            setPair(pairable, top);
            pqzero.erase(pairable);

            // t.reStart();
            int sz = sizeCofacets(pairable, pairable.first + 1, triangulation);
            // relationTime_ += t.getElapsedTime();
            for(int i = 0; i < sz; i++) {
              // t.reStart();
              int coface = extractCoboundary(
                pairable, pairable.first + 1, i, triangulation);
              // relationTime_ += t.getElapsedTime();

              Simplex newS = Simplex(pairable.first + 1, coface);
              int index = getIndex(newS, triangulation);
              if(index == indexing_[v]) {
                Simplex other;
                numPairable
                  = numPairableFaces(newS, critical, other, triangulation);
                if(numPairable == 0)
                  pqzero.insert(newS);
                else if(numPairable == 1)
                  pqone.insert(newS);
              }
            }

            // t.reStart();
            sz = sizeCofacets(top, top.first + 1, triangulation);
            // relationTime_ += t.getElapsedTime();
            for(int i = 0; i < sz; i++) {
              // t.reStart();
              int coface
                = extractCoboundary(top, top.first + 1, i, triangulation);
              // relationTime_ += t.getElapsedTime();

              Simplex newS = Simplex(top.first + 1, coface);
              int index = getIndex(newS, triangulation);
              if(index == indexing_[v]) {
                numPairable
                  = numPairableFaces(newS, critical, pairable, triangulation);
                if(numPairable == 0)
                  pqzero.insert(newS);
                else if(numPairable == 1)
                  pqone.insert(newS);
              }
            }
          }
        }

        // time for new critical
        if(!pqzero.empty()) {
          Simplex top = *pqzero.begin();
          pqzero.erase(top);

          if(isCritical(top)) {
            critical.insert(top);

            // t.reStart();
            int sz = sizeCofacets(top, top.first + 1, triangulation);
            // relationTime_ += t.getElapsedTime();
            for(int i = 0; i < sz; i++) {
              // t.reStart();
              int coface
                = extractCoboundary(top, top.first + 1, i, triangulation);
              // relationTime_ += t.getElapsedTime();

              Simplex newS = Simplex(top.first + 1, coface);
              int index = getIndex(newS, triangulation);
              if(index == indexing_[v]) {
                Simplex pairable;
                int numPairable
                  = numPairableFaces(newS, critical, pairable, triangulation);
                if(numPairable == 0)
                  pqzero.insert(newS);
                else if(numPairable == 1)
                  pqone.insert(newS);
              }
            }
          }
        }
      }
    }

    template <typename triangulationType>
    int getIndex(Simplex simplex, const triangulationType &triangulation) {
      // from the simplex extract its boundary vertices
      vector<SimplexId> vertices;
      simplexToVertices(simplex, vertices, triangulation);

      int index_val = indexing_[vertices[0]];

      for(SimplexId v = 1; v < (SimplexId)vertices.size(); v++) {
        index_val = index_val > indexing_[vertices[v]] ? index_val
                                                       : indexing_[vertices[v]];
      }

      return index_val;
    }

    void setPair(const Simplex &tail, const Simplex &head);
    bool getPair(const Simplex &simpl, Simplex &paired);
    bool isCritical(const Simplex &simpl);

    template <typename triangulationType>
    int numPairableFaces(Simplex simplex,
                         SimplexesSet &critical,
                         Simplex &chosen,
                         const triangulationType &triangulation) {
      int tot = 0;
      int index = getIndex(simplex, triangulation);

      for(int i = 0; i < simplex.first + 1; i++) {
        // Timer t;
        int face
          = extractBoundary(simplex, simplex.first - 1, i, triangulation);
        // relationTime_ += t.getElapsedTime();
        Simplex simplB = Simplex(simplex.first - 1, face);
        int indexb = getIndex(simplB, triangulation);

        if(isCritical(simplB) && critical.find(simplB) == critical.end()
           && indexb == index) {
          tot++;
          chosen = simplB;
        }
      }

      return tot;
    }

    template <typename triangulationType>
    bool cmpSimplexesLevelSets(const Simplex &s1,
                               const Simplex &s2,
                               const triangulationType &triangulation) {
      if(s1.first == s2.first) {
        vector<SimplexId> indexing1;
        vector<SimplexId> indexing2;

        getVertexIndexingOfVertices(s1, indexing1, triangulation);
        getVertexIndexingOfVertices(s2, indexing2, triangulation);

        sort(indexing1.begin(), indexing1.end(), std::greater<int>());
        sort(indexing2.begin(), indexing2.end(), std::greater<int>());

        for(int i = 0; i < (SimplexId)indexing1.size(); i++) {
          if(indexing1[i] == indexing2[i]) {
            continue;
          }

          return indexing1[i] < indexing2[i];
        }
      }
      return s1.first < s2.first;
    }

    template <typename triangulationType>
    void getVertexIndexingOfVertices(Simplex simplex,
                                     vector<SimplexId> &vec,
                                     const triangulationType &triangulation) {
      simplexToVertices(simplex, vec, triangulation);
      for(SimplexId v = 0; v < (SimplexId)vec.size(); v++) {
        vec[v] = indexing_[vec[v]];
      }
    }

    template <typename triangulationType>
    void computeBarycenter(Simplex &simplex,
                           vector<float> &coords,
                           const triangulationType &triangulation) {
      coords = vector<float>(3, 0);
      // Timer t;
      for(int i = 0; i < simplex.first + 1; i++) {
        SimplexId vertex;

        // t.reStart();
        switch(simplex.first) {
          case 0:
            vertex = simplex.second;
            break;
          case 1:
            triangulation->getEdgeVertex(simplex.second, i, vertex);
            break;
          case 2:
            triangulation->getTriangleVertex(simplex.second, i, vertex);
            break;
          case 3:
            triangulation->getCellVertex(simplex.second, i, vertex);
            break;
        }
        // relationTime_ += t.getElapsedTime();

        float x, y, z;
        triangulation->getVertexPoint(vertex, x, y, z);
        coords[0] += x;
        coords[1] += y;
        coords[2] += z;
      }

      for(int i = 0; i < 3; i++)
        coords[i] /= (float)simplex.first + 1;
    }

  protected:
    int dimensionality_;
    double relationTime_;
    vector<SimplexId> indexing_;
    vector<vector<vector<SimplexId>>> gradient_;
    // std::ofstream logFile_;
  };
} // namespace ttk