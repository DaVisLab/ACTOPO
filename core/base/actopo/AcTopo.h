
#pragma once

// base code includes
#include <AbstractTriangulation.h>
#include <CellArray.h>
#include <FlatJaggedArray.h>
#include <bitset>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <list>
#include <map>
#include <set>
#include <string>

// parallel library
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <semaphore.h>
#include <thread>

#define EDGE_ID 1
#define TRIANGLE_ID 2
#define PRINT_WAITING_TIME 1

namespace ttk {

  class ImplicitCluster {
  private:
    /* components */
    SimplexId nid;
    /* boundary cells */
    std::vector<std::vector<bool>> boundaryVertices_;
    std::vector<std::vector<bool>> boundaryEdges_;
    std::vector<std::vector<bool>> boundaryTriangles_;
    /* vertex relations */
    std::vector<FlatJaggedArray> vertexEdges_;
    std::vector<FlatJaggedArray> vertexLinks_;
    std::vector<FlatJaggedArray> vertexNeighbors_;
    std::vector<FlatJaggedArray> vertexStars_;
    std::vector<FlatJaggedArray> vertexTriangles_;
    /* edge relations */
    // edgeVertex relation can be extracted from internal edge list
    std::vector<FlatJaggedArray> edgeLinks_;
    std::vector<FlatJaggedArray> edgeStars_;
    std::vector<FlatJaggedArray> edgeTriangles_;
    /* triangle relations */
    // triangleVertex relation can be extracted from internal triangle list
    std::vector<std::vector<std::array<SimplexId, 3>>> triangleEdges_;
    std::vector<FlatJaggedArray> triangleLinks_;
    std::vector<FlatJaggedArray> triangleStars_;
    /* cell relations */
    std::vector<std::vector<std::array<SimplexId, 6>>> tetraEdges_;
    std::vector<FlatJaggedArray> cellNeighbors_;
    std::vector<std::vector<std::array<SimplexId, 4>>> tetraTriangles_;

  public:
    ImplicitCluster() {
    }
    ImplicitCluster(SimplexId id, SimplexId num) : nid(id) {
      /* boundary cells */
      boundaryVertices_ = std::vector<std::vector<bool>>(num);
      boundaryEdges_ = std::vector<std::vector<bool>>(num);
      boundaryTriangles_ = std::vector<std::vector<bool>>(num);
      /* vertex relations */
      vertexEdges_ = std::vector<FlatJaggedArray>(num);
      vertexLinks_ = std::vector<FlatJaggedArray>(num);
      vertexNeighbors_ = std::vector<FlatJaggedArray>(num);
      vertexStars_ = std::vector<FlatJaggedArray>(num);
      vertexTriangles_ = std::vector<FlatJaggedArray>(num);
      /* edge relations */
      // edgeVertex relation can be extracted from internal edge list
      edgeLinks_ = std::vector<FlatJaggedArray>(num);
      edgeStars_ = std::vector<FlatJaggedArray>(num);
      edgeTriangles_ = std::vector<FlatJaggedArray>(num);
      /* triangle relations */
      // triangleVertex relation can be extracted from internal triangle list
      triangleEdges_ = std::vector<std::vector<std::array<SimplexId, 3>>>(num);
      triangleLinks_ = std::vector<FlatJaggedArray>(num);
      triangleStars_ = std::vector<FlatJaggedArray>(num);
      /* cell relations */
      tetraEdges_ = std::vector<std::vector<std::array<SimplexId, 6>>>(num);
      cellNeighbors_ = std::vector<FlatJaggedArray>(num);
      tetraTriangles_ = std::vector<std::vector<std::array<SimplexId, 4>>>(num);
    }
    ~ImplicitCluster() {
    }

    inline void clear(const ThreadId consumerId) {
      // keep the edge and triangle lists
      /* boundary cells */
      boundaryVertices_[consumerId] = std::vector<bool>{};
      boundaryEdges_[consumerId] = std::vector<bool>{};
      // keep boundary triangles 
      // boundaryTriangles_ = std::vector<bool>{};
      /* vertex relations */
      vertexEdges_[consumerId] = FlatJaggedArray{};
      vertexLinks_[consumerId] = FlatJaggedArray{};
      vertexNeighbors_[consumerId] = FlatJaggedArray{};
      vertexStars_[consumerId] = FlatJaggedArray{};
      vertexTriangles_[consumerId] = FlatJaggedArray{};
      /* edge relations */
      // edgeVertex relation can be extracted from internal edge list
      edgeLinks_[consumerId] = FlatJaggedArray{};
      edgeStars_[consumerId] = FlatJaggedArray{};
      edgeTriangles_[consumerId] = FlatJaggedArray{};
      /* triangle relations */
      // triangleVertex relation can be extracted from internal triangle list
      triangleEdges_[consumerId] = std::vector<std::array<SimplexId, 3>>{};
      triangleLinks_[consumerId] = FlatJaggedArray{};
      triangleStars_[consumerId] = FlatJaggedArray{};
      /* cell relations */
      tetraEdges_[consumerId] = std::vector<std::array<SimplexId, 6>>{};
      cellNeighbors_[consumerId] = FlatJaggedArray{};
      tetraTriangles_[consumerId] = std::vector<std::array<SimplexId, 4>>{};
    }

    friend class AcTopo;
  };

  class AcTopo final : public AbstractTriangulation {

  public:
    #ifdef PRINT_WAITING_TIME
    // Keep track of timings
    mutable std::vector<double> waitingTimes_;
    mutable std::vector<int> missCnts_;
    #endif
    mutable std::ofstream consumerFile_, producerFile_, logFile_;

    const std::vector<std::string> RelationNames = {
      "BlankRelation", // Empty
      "VVRelation", // vertexNeighbor
      "VTRelation", // vertexStar
      "EVRelation", // edgeVertex
      "ETRelation", // edgeStar
      "FVRelation", // triangleVertex
      "FTRelation", // triangleStar
      "TVRelation", // cellVertex
      "TTRelation", // cellNeighbor
      // external maps computation
      "IEMapComp", // internalEdgeMap
      "IFMapComp", // internalTriangleMap
      "IEFMapComp", // intenralEdgeMap and internalTriangleMap
      "IFBMapComp", // internalTriangleMap and boundaryCell
      // require external maps
      "VERelation", // vertexEdge
      "VFRelation", // vertexTriangle
      "EFRelation", // edgeTriangle
      "FERelation", // triangleEdge
      "TERelation", // cellEdge
      "TFRelation", // cellTriangle
      "VLRelation", // vertexLink
      "ELRelation", // edgeLink
      "FLRelation", // triangleLink
      "BVRelation", // boundaryVertex
      "BERelation", // boundaryEdge
      "BFRelation", // boundaryTriangle
      "RelationNum" // total number of relations
    };

    AcTopo();

    AcTopo(const AcTopo &rhs);

    AcTopo &operator=(const AcTopo &rhs);

    ~AcTopo();

    /**
     * Print out the clusters in the buffer.
     */
    // inline void printBuffer(SimplexId clusterId,
    //                  std::vector<SimplexId> &removed) const {
    //   std::cout << clusterId << ", " << currCluster_ << ": [ ";
    //   for(size_t i = 0; i < sbuffer_.size(); i++) {
    //     std::cout << sbuffer_[i] << " ";
    //   }
    //   std::cout << "]" << std::endl;
    //   std::cout << "removed: [";
    //   for(SimplexId r : removed) {
    //     std::cout << r << " ";
    //   }
    //   std::cout << "]" << std::endl;
    // }

    /**
     * Print out the vector.
     */
    inline void printVector(std::vector<SimplexId> &vec, std::string name) const {
      std::cout << name << ": [ ";
      for(size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
      }
      std::cout << "]" << std::endl;
    }

    /**
     * Set up vertices from the input.
     */
    inline int setInputPoints(const SimplexId &pointNumber,
                              const void *pointSet,
                              const int *indexArray,
                              const bool &doublePrecision = false) {

      if(vertexNumber_)
        clear();

      vertexNumber_ = pointNumber;
      pointSet_ = pointSet;
      vertexIndices_ = indexArray;
      doublePrecision_ = doublePrecision;

      return 0;
    }

    /**
     * Set up cells from the input.
     */

#ifdef TTK_CELL_ARRAY_NEW
    // Layout with connectivity + offset array (new)
    int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *connectivity,
                             const LongSimplexId *offset);
#else
    // Flat layout with a single array (legacy & default one)
    int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray);
#endif

    /**
     * Reorder the input vertices.
     */
    int reorderVertices(std::vector<SimplexId> &vertexMap);

    /**
     * Reorder the input cells.
     */
    int reorderCells(const std::vector<SimplexId> &vertexMap,
                     const SimplexId &cellNumber,
                     const LongSimplexId *connectivity,
                     const LongSimplexId *offset);

    /**
     * Reorder the input cells.
     */
    int reorderCells(const std::vector<SimplexId> &vertexMap,
                     const LongSimplexId *cellArray);

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if(localEdgeId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].tetraEdges_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TERelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localEdgeId >= (int)(allClusters_[nid].tetraEdges_)[localCellId].size())
        return -2;
#endif
      edgeId = (allClusters_[nid].tetraEdges_[tid])[localCellId][localEdgeId];
      return 0;
    }

    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      (void)cellId;
      return (maxCellDim_ + 1) * maxCellDim_ / 2;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellEdgesInternal() override {
      if(cellEdgeVector_.empty()) {
        cellEdgeVector_.reserve(cellNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          if(allClusters_[nid].tetraEdges_[0].empty()) {
            getClusterCellEdges(nid, 0);
          }
          for(size_t i = 0; i < allClusters_[nid].tetraEdges_[0].size(); i++) {
            cellEdgeVector_.push_back({allClusters_[nid].tetraEdges_[0][i].begin(),
                                       allClusters_[nid].tetraEdges_[0][i].end()});
          }
        }
      }
      return &cellEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if(localNeighborId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].cellNeighbors_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TTRelation, tid);
      }
      
#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId >= allClusters_[nid].cellNeighbors_.size(localCellId))  
        return -2;
#endif
      neighborId = allClusters_[nid].cellNeighbors_[tid].get(localCellId, localNeighborId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].cellNeighbors_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TTRelation, tid);
      }

      return allClusters_[nid].cellNeighbors_[tid].size(localCellId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      cellNeighborList_.reserve(cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localCellNeighbors;
        if(allClusters_[nid].cellNeighbors_[0].empty()) {
          getClusterCellNeighbors(nid, 0);
        }
        allClusters_[nid].cellNeighbors_[0].copyTo(localCellNeighbors);
        cellNeighborList_.insert(cellNeighborList_.end(),
                                 localCellNeighbors.begin(),
                                 localCellNeighbors.end());
      }
      return &cellNeighborList_;
    }

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId >= getCellTriangleNumber(cellId)))
        return -2;
#endif

      SimplexId nid = vertexIndices_[cellArray_->getCellVertex(cellId, 0)];
      SimplexId localCellId = cellId - cellIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].tetraTriangles_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TFRelation, tid);
      }

      triangleId = (allClusters_[nid].tetraTriangles_[tid])[localCellId][localTriangleId];
      return 0;
    }

    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      (void)cellId;
      return (maxCellDim_ + 1) * maxCellDim_ * (maxCellDim_ - 1) / 6;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override {
      if(cellTriangleVector_.empty()) {
        cellTriangleVector_.reserve(cellNumber_);
        std::bitset<5> type(8);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          if(allClusters_[nid].tetraTriangles_[0].empty()) {
            getClusterCellTriangles(nid, 0);
          }
          for(size_t i = 0; i < allClusters_[nid].tetraTriangles_[0].size(); i++) {
            cellTriangleVector_.push_back(
              {allClusters_[nid].tetraTriangles_[0][i].begin(),
               allClusters_[nid].tetraTriangles_[0][i].end()});
          }
        }
      }
      return &cellTriangleVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localVertexId < 0)
         || (localVertexId >= cellArray_->getCellVertexNumber(cellId)))
        return -2;
#endif

      vertexId = cellArray_->getCellVertex(cellId, localVertexId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
#endif

      return cellArray_->getCellVertexNumber(cellId);
    }

    int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const override {
      return maxCellDim_;
    }

    inline const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override {
      return &edgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
      const SimplexId &edgeId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
      if(localLinkId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ELRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId >= allClusters_[nid].edgeLinks_[tid].size(localEdgeId))
        return -2;
#endif
      linkId = allClusters_[nid].edgeLinks_[tid].get(localEdgeId, localLinkId);

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeLinks_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ELRelation, tid);
      }
      return allClusters_[nid].edgeLinks_[tid].size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      edgeLinkList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeLinks;
        if(allClusters_[nid].edgeLinks_[0].empty()) {
          getClusterEdgeLinks(nid, 0);
        }
        allClusters_[nid].edgeLinks_[0].copyTo(localEdgeLinks);
        edgeLinkList_.insert(
          edgeLinkList_.end(), localEdgeLinks.begin(), localEdgeLinks.end());
      }
      return &edgeLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
      const SimplexId &edgeId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
      if(localStarId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;
      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ETRelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId >= allClusters_[nid].edgeStars_[tid].size(localEdgeId))
        return -2;
#endif
      starId = allClusters_[nid].edgeStars_[tid].get(localEdgeId, localStarId);

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeStars_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ETRelation, tid);
      }
      return allClusters_[nid].edgeStars_[tid].size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      edgeStarList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeStars;
        if(allClusters_[nid].edgeStars_[0].empty()) {
          getClusterEdgeStars(nid, 0);
        }
        allClusters_[nid].edgeStars_[0].copyTo(localEdgeStars);
        edgeStarList_.insert(
          edgeStarList_.end(), localEdgeStars.begin(), localEdgeStars.end());
      }
      return &edgeStarList_;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
      if(localTriangleId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeTriangles_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::EFRelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localTriangleId >= allClusters_[nid].edgeTriangles_[tid].size(localEdgeId))
        return -2;
#endif
      triangleId = allClusters_[nid].edgeTriangles_[tid].get(localEdgeId, localTriangleId);
      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeTriangles_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::EFRelation, tid);
      }

      return allClusters_[nid].edgeTriangles_[tid].size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override {
      edgeTriangleList_.reserve(edgeIntervals_.back() + 1);
      std::bitset<5> type(2);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeTriangles;
        if(allClusters_[nid].edgeTriangles_[0].empty()) {
          getClusterEdgeTriangles(nid, 0);
        }
        allClusters_[nid].edgeTriangles_[0].copyTo(localEdgeTriangles);
        edgeTriangleList_.insert(edgeTriangleList_.end(),
                                 localEdgeTriangles.begin(),
                                 localEdgeTriangles.end());
      }
      return &edgeTriangleList_;
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > (SimplexId)edgeIntervals_.back()))
        return -1;
      if((localVertexId != 0) && (localVertexId != 1))
        return -2;
#endif

      if(localVertexId) {
        vertexId = edgeList_[edgeId][1];
      } else {
        vertexId = edgeList_[edgeId][0];
      }
      return 0;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const override {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdgesInternal() const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!edgeIntervals_.size())
        return -1;
#endif
      return edgeIntervals_.back() + 1;
    }

    inline SimplexId getNumberOfTrianglesInternal() const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!triangleIntervals_.size())
        return -1;
#endif

      return triangleIntervals_.back() + 1;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const override {
      return vertexNumber_;
    }

    inline const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override {
      // if it is a triangle mesh
      if(maxCellDim_ == 2) {
        triangleList_.resize(cellNumber_, std::array<SimplexId, 3>());
        for(SimplexId cid = 0; cid < cellNumber_; cid++) {
          triangleList_[cid][0] = cellArray_->getCellVertex(cid, 0);
          triangleList_[cid][1] = cellArray_->getCellVertex(cid, 1);
          triangleList_[cid][2] = cellArray_->getCellVertex(cid, 2);
        }
      }
      return &triangleList_;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
      if((localEdgeId < 0) || (localEdgeId > 2))
        return -2;
#endif

      SimplexId nid = vertexIndices_[triangleList_[triangleId][0]];
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleEdges_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FERelation, tid);
      }

      edgeId = (allClusters_[nid].triangleEdges_[tid])[localTriangleId][localEdgeId];
      return 0;
    }

    inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      (void)triangleId;
      return 3;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() override {
      if(triangleEdgeVector_.empty()) {
        triangleEdgeVector_.reserve(triangleIntervals_.size() + 1);
        std::bitset<5> type(4);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          if(allClusters_[nid].triangleEdges_[0].empty()) {
            getClusterTriangleEdges(nid, 0);
          }
          for(size_t i = 0; i < allClusters_[nid].triangleEdges_[0].size(); i++) {
            triangleEdgeVector_.push_back(
              {allClusters_[nid].triangleEdges_[0][i].begin(),
               allClusters_[nid].triangleEdges_[0][i].end()});
          }
        }
      }
      return &triangleEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
      if(localLinkId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[triangleList_[triangleId][0]];
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleLinks_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FLRelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId >= allClusters_[nid].triangleLinks_[tid].size(localTriangleId))
        return -2;
#endif
      linkId = allClusters_[nid].triangleLinks_[tid].get(localTriangleId, localLinkId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[triangleList_[triangleId][0]];
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleLinks_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FLRelation, tid);
      }

      return allClusters_[nid].triangleLinks_[tid].size(localTriangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      triangleLinkList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localTriangleLinks;
        if(allClusters_[nid].triangleLinks_[0].empty()) {
          getClusterTriangleLinks(nid, 0);
        }
        allClusters_[nid].triangleLinks_[0].copyTo(localTriangleLinks);
        triangleLinkList_.insert(triangleLinkList_.end(),
                                 localTriangleLinks.begin(),
                                 localTriangleLinks.end());
      }
      return &triangleLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || !triangleIntervals_.size()
         || (triangleId > triangleIntervals_.back()))
        return -1;
      if(localStarId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[triangleList_[triangleId][0]];
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleStars_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FTRelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId >= allClusters_[nid].triangleStars_[tid].size(localTriangleId))
        return -2;
#endif
      starId = allClusters_[nid].triangleStars_[tid].get(localTriangleId, localStarId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || !triangleIntervals_.size()
         || (triangleId > triangleIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[triangleList_[triangleId][0]];
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleStars_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FTRelation, tid);
      }
      
      return allClusters_[nid].triangleStars_[tid].size(localTriangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
      triangleStarList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localTriangleStars;
        if(allClusters_[nid].triangleStars_[0].empty()) {
          getClusterTriangleStars(nid, 0);
        }
        allClusters_[nid].triangleStars_[0].copyTo(localTriangleStars);
        triangleStarList_.insert(triangleStarList_.end(),
                                 localTriangleStars.begin(),
                                 localTriangleStars.end());
      }
      return &triangleStarList_;
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return -1;
      if((localVertexId < 0) || (localVertexId > 2))
        return -2;
#endif

      vertexId = triangleList_[triangleId][localVertexId];
      return 0;
    }

    inline int getVertexEdgeInternal(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localEdgeId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexEdges_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VERelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localEdgeId >= allClusters_[nid].vertexEdges_[tid].size(localVertexId))
        return -2;
#endif
      edgeId = allClusters_[nid].vertexEdges_[tid].get(localVertexId, localEdgeId);

      return 0;
    }

    inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexEdges_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VERelation, tid);
      }
      return allClusters_[nid].vertexEdges_[tid].size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override {
      vertexEdgeList_.reserve(vertexNumber_);
      std::bitset<5> type(1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexEdges;
        if(allClusters_[nid].vertexEdges_[0].empty()) {
          getClusterVertexEdges(nid, 0);
        }
        allClusters_[nid].vertexEdges_[0].copyTo(localVertexEdges);
        vertexEdgeList_.insert(vertexEdgeList_.end(), localVertexEdges.begin(),
                               localVertexEdges.end());
      }

      return &vertexEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId > vertexIntervals_.back()))
        return -1;
      if(localLinkId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;
      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexLinks_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VLRelation, tid);
      }

      if(localLinkId >= allClusters_[nid].vertexLinks_[tid].size(localVertexId)) {
        linkId = -2;
      } else {
        linkId = allClusters_[nid].vertexLinks_[tid].get(localVertexId, localLinkId);
      }

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId > vertexIntervals_.back()))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexLinks_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VLRelation, tid);
      }

      return allClusters_[nid].vertexLinks_[tid].size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      vertexLinkList_.reserve(vertexIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexLinks;
        if(allClusters_[nid].vertexLinks_[0].empty()) {
          getClusterVertexLinks(nid, 0);
        }
        allClusters_[nid].vertexLinks_[0].copyTo(localVertexLinks);
        vertexLinkList_.insert(vertexLinkList_.end(), localVertexLinks.begin(),
                               localVertexLinks.end());
      }
      return &vertexLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_)) 
        return -1;
      if(localNeighborId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexNeighbors_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VVRelation, tid);
      }
      
#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId >= allClusters_[nid].vertexNeighbors_[tid].size(localVertexId)) {
        return -2;
#endif
      neighborId
        = allClusters_[nid].vertexNeighbors_[tid].get(localVertexId, localNeighborId);

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexNeighbors_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VVRelation, tid);
      }

      return allClusters_[nid].vertexNeighbors_[tid].size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
      vertexNeighborList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexNeighbors;
        if(allClusters_[nid].vertexNeighbors_[0].empty()) {
          getClusterVertexNeighbors(nid, 0);
        }
        allClusters_[nid].vertexNeighbors_[0].copyTo(localVertexNeighbors);
        vertexNeighborList_.insert(vertexNeighborList_.end(),
                                   localVertexNeighbors.begin(),
                                   localVertexNeighbors.end());
      }
      return &vertexNeighborList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
      const SimplexId &vertexId, float &x, float &y, float &z) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      Timer t;
      if(doublePrecision_) {
        x = ((double *)pointSet_)[3 * vertexId];
        y = ((double *)pointSet_)[3 * vertexId + 1];
        z = ((double *)pointSet_)[3 * vertexId + 2];
      } else {
        x = ((float *)pointSet_)[3 * vertexId];
        y = ((float *)pointSet_)[3 * vertexId + 1];
        z = ((float *)pointSet_)[3 * vertexId + 2];
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexStar)(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localStarId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexStars_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VTRelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId >= allClusters_[nid].vertexStars_[tid].size(localVertexId))
        return -2;
#endif
      starId = allClusters_[nid].vertexStars_[tid].get(localVertexId, localStarId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexStars_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VTRelation, tid);
      }

      return allClusters_[nid].vertexStars_[tid].size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      vertexStarList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexStars;
        if(allClusters_[nid].vertexStars_[0].empty()) {
          getClusterVertexStars(nid, 0);
        }
        allClusters_[nid].vertexStars_[0].copyTo(localVertexStars);
        vertexStarList_.insert(vertexStarList_.end(), localVertexStars.begin(),
                               localVertexStars.end());
      }
      return &vertexStarList_;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
      if(localTriangleId < 0)
        return -2;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexTriangles_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VFRelation, tid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localTriangleId >= allClusters_[nid].vertexTriangles_[tid].size(localVertexId)) {
        return -2;
#endif
      triangleId
        = allClusters_[nid].vertexTriangles_[tid].get(localVertexId, localTriangleId);

      return 0;
    }

    inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexTriangles_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VFRelation, tid);
      }

      return allClusters_[nid].vertexTriangles_[tid].size(localVertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override {
      vertexTriangleList_.reserve(vertexNumber_);
      std::bitset<5> type(2);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexTriangles;
        if(allClusters_[nid].vertexTriangles_[0].empty()) {
          getClusterVertexTriangles(nid, 0);
        }
        allClusters_[nid].vertexTriangles_[0].copyTo(localVertexTriangles);
        vertexTriangleList_.insert(vertexTriangleList_.end(),
                                   localVertexTriangles.begin(),
                                   localVertexTriangles.end());
      }
      return &vertexTriangleList_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId > edgeIntervals_.back()))
        return false;
#endif
      SimplexId nid = vertexIndices_[edgeList_[edgeId][0]];
      SimplexId localEdgeId = edgeId - edgeIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].boundaryEdges_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::BERelation, tid);
      }

      return (allClusters_[nid].boundaryEdges_[tid])[localEdgeId];
    }

    bool isEmpty() const override {
      return !vertexNumber_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const override {
      if(maxCellDim_ == 2)
        return false;

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId > triangleIntervals_.back()))
        return false;
#endif
      SimplexId nid = vertexIndices_[triangleList_[triangleId][0]];
      SimplexId localTriangleId = triangleId - triangleIntervals_[nid - 1] - 1;

      return (allClusters_[nid].boundaryTriangles_[0])[localTriangleId];
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return false;
#endif
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      ThreadId tid = 0;
      #ifdef TTK_ENABLE_OPENMP
      tid = omp_get_thread_num();
      #endif
      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].boundaryVertices_[tid].empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::BVRelation, tid);
      }

      return (allClusters_[nid].boundaryVertices_[tid])[localVertexId];
    }

    inline int preconditionBoundaryEdgesInternal() override {
      if(maxCellDim_ == 2) {
        Timer t;
        if(edgeList_.empty()) {
          clusterEdges_.resize(nodeNumber_ + 1);
          edgeIntervals_.resize(nodeNumber_ + 1);
          edgeIntervals_[0] = -1;

          // use the producer threads to compute the edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(int i = 0; i < threadNumber_; i++) {
          std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
          reqRelations_[i] = RelationType::IBEList;
          preconditions_[i] = 1;
          llck.unlock();
          leaderCondVars_[i].notify_one();
        }

        SimplexId startNodeId = clustersPerThread_ * (threadNumber_ * numProducers_) + 1;
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif
          for(SimplexId nid = startNodeId; nid <= nodeNumber_; nid++) {
            edgeIntervals_[nid] = buildInternalEdgeList(nid, true);
          }

          for(int i = 0; i < threadNumber_; i++)
            sem_wait(&semaphores_[i*numProducers_]);
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
            edgeIntervals_[nid]
              = edgeIntervals_[nid - 1] + edgeIntervals_[nid];
          }

          edgeList_.resize(edgeIntervals_.back() + 1);
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
            int j = edgeIntervals_[nid-1];
            for (size_t k = 0; k < clusterEdges_[nid].size(); k++) {
              edgeList_[++j] = std::move(clusterEdges_[nid][k]);
            }
          }

          clusterEdges_.clear(); clusterEdges_.shrink_to_fit();
          printMsg("Edges and boundary edges preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
        } else {
          // use the producer threads to compute boundary edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
          for(int i = 0; i < threadNumber_; i++) {
            std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
            reqRelations_[i] = RelationType::IBCList;
            preconditions_[i] = 1;
            llck.unlock();
            leaderCondVars_[i].notify_one();
          }
          
          SimplexId startNodeId = clustersPerThread_ * (threadNumber_ * numProducers_) + 1;
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif
          for(SimplexId nid = startNodeId; nid <= nodeNumber_; nid++) {
            buildBoundaryTopCellList(nid);
          }

          for(int i = 0; i < threadNumber_; i++)
            sem_wait(&semaphores_[i*numProducers_]);
          printMsg("Boundary edges preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
        }

      } else if(maxCellDim_ == 3) {
        preconditionBoundaryTriangles();
        relationVec_.push_back(RelationType::BERelation);
      } else {
        // unsupported dimension
        printErr(
          "[AcTopo] Unsupported dimension for boundary precondition.");
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() override {
      if(maxCellDim_ == 2) {
        return 0;
      } else if(maxCellDim_ == 3) {
        Timer t;
        if(triangleList_.empty()) {
          clusterTriangles_.resize(nodeNumber_ + 1);
          triangleIntervals_.resize(nodeNumber_ + 1);
          triangleIntervals_[0] = -1;

          // use the producer threads to compute the triangles
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        // use the producer threads to compute the edges
        for(int i = 0; i < threadNumber_; i++) {
          std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
          reqRelations_[i] = RelationType::IBFList;
          preconditions_[i] = 1;
          llck.unlock();
          leaderCondVars_[i].notify_one();
        }

          SimplexId startNodeId = clustersPerThread_ * (threadNumber_ * numProducers_) + 1;
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif
          for(SimplexId nid = startNodeId; nid <= nodeNumber_; nid++) {
            triangleIntervals_[nid] = buildInternalTriangleList(nid, true);
          }

          for(int i = 0; i < threadNumber_; i++) 
            sem_wait(&semaphores_[i*numProducers_]);

          for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
            triangleIntervals_[nid]
              = triangleIntervals_[nid - 1] + triangleIntervals_[nid];
          }

          triangleList_.resize(triangleIntervals_.back() + 1);
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
            int j = triangleIntervals_[nid-1];
            for (size_t k = 0; k < clusterTriangles_[nid].size(); k++) {
              triangleList_[++j] = std::move(clusterTriangles_[nid][k]);
            }
          }
          clusterTriangles_.clear(); clusterTriangles_.shrink_to_fit();
          printMsg("Triangles and boundary triangles preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
        } else {
          // use the producer threads to compute boundary triangles
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        // use the producer threads to compute the edges
        for(int i = 0; i < threadNumber_; i++) {
          std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
          reqRelations_[i] = RelationType::IBCList;
          preconditions_[i] = 1;
          llck.unlock();
          leaderCondVars_[i].notify_one();
        }
          
          SimplexId startNodeId = clustersPerThread_ * (threadNumber_ * numProducers_) + 1;
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_)
  #endif
          for(SimplexId nid = startNodeId; nid <= nodeNumber_; nid++) {
            buildBoundaryTopCellList(nid);
          }
          for(int i = 0; i < threadNumber_; i++) 
            sem_wait(&semaphores_[i*numProducers_]);
          printMsg("Boundary triangles preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
        }
      } else {
        // unsupported dimension
        printErr(
          "[AcTopo] Unsupported dimension for boundary precondition.");
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryVerticesInternal() override {
      if(maxCellDim_ == 2) {
        if(!hasPreconditionedBoundaryEdges_) {
          preconditionBoundaryEdgesInternal();
          hasPreconditionedBoundaryEdges_ = true;
        }
      } else if(maxCellDim_ == 3) {
        if(!hasPreconditionedBoundaryTriangles_) {
          preconditionBoundaryTrianglesInternal();
          hasPreconditionedBoundaryTriangles_ = true;
        }
      }
      relationVec_.push_back(RelationType::BVRelation);
      return 0;
    }

    inline int preconditionCellEdgesInternal() override {
      relationVec_.push_back(RelationType::TERelation);
      return 0;
    }

    inline int preconditionCellNeighborsInternal() override {
      relationVec_.push_back(RelationType::TTRelation);
      return 0;
    }

    inline int preconditionCellTrianglesInternal() override {
      relationVec_.push_back(RelationType::TFRelation);
      return 0;
    }

    inline int preconditionEdgesInternal() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexNumber_ <= 0)
        return -1;
      if(cellNumber_ <= 0)
        return -2;
      if(!cellArray_)
        return -3;
      if(nodeNumber_ <= 0)
        return -4;
#endif

      if(edgeList_.empty()) {
        Timer t;
        clusterEdges_.resize(nodeNumber_ + 1);
        edgeIntervals_.resize(nodeNumber_ + 1);
        edgeIntervals_[0] = -1;

        // use the producer threads to compute the edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(int i = 0; i < threadNumber_; i++) {
          std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
          reqRelations_[i] = RelationType::IEList;
          preconditions_[i] = 1;
          llck.unlock();
          leaderCondVars_[i].notify_one();
        }

        SimplexId startNodeId = clustersPerThread_ * (threadNumber_ * numProducers_) + 1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(SimplexId nid = startNodeId; nid <= nodeNumber_; nid++) {
          edgeIntervals_[nid] = buildInternalEdgeList(nid);
        }

        for(int i = 0; i < threadNumber_; i++)
          sem_wait(&semaphores_[i*numProducers_]);
        
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          edgeIntervals_[nid] = edgeIntervals_[nid - 1] + edgeIntervals_[nid];
        }
        edgeList_.reserve(edgeIntervals_.back()+10);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          for(auto &edge : clusterEdges_[nid]) {
            edgeList_.push_back(std::move(edge));
          }
        }
        clusterEdges_.clear(); clusterEdges_.shrink_to_fit();
        printMsg("Edges preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
      }

      return 0;
    }

    inline int preconditionEdgeLinksInternal() override {
      if(maxCellDim_ == 2 || maxCellDim_ == 3) {
        preconditionEdges();
        relationVec_.push_back(RelationType::ELRelation);
      } else {
        // unsupported dimension
        printErr("Unsupported dimension for edge link precondition.");
        return -1;
      }
      return 0;
    }

    inline int preconditionEdgeStarsInternal() override {
      relationVec_.push_back(RelationType::ETRelation);
      return 0;
    }

    inline int preconditionEdgeTrianglesInternal() override {
      relationVec_.push_back(RelationType::EFRelation);
      return 0;
    }

    inline int preconditionTrianglesInternal() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexNumber_ <= 0)
        return -1;
      if(cellNumber_ <= 0)
        return -2;
      if(!cellArray_)
        return -3;
      if(nodeNumber_ <= 0)
        return -4;
#endif

      // build triangle interval list
      if(triangleList_.empty()) {
        Timer t;
        clusterTriangles_.resize(nodeNumber_ + 1);
        triangleIntervals_.resize(nodeNumber_ + 1);
        triangleIntervals_[0] = -1;

        // use the producer threads to compute the triangles
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(int i = 0; i < threadNumber_; i++) {
          std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
          reqRelations_[i] = RelationType::IFList;
          preconditions_[i] = 1;
          llck.unlock();
          leaderCondVars_[i].notify_one();
        }
        
        SimplexId startNodeId = clustersPerThread_ * (threadNumber_ * numProducers_) + 1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
        for(SimplexId nid = startNodeId; nid <= nodeNumber_; nid++) {
          triangleIntervals_[nid] = buildInternalTriangleList(nid);
        }

        for(int i = 0; i < threadNumber_; i++)
          sem_wait(&semaphores_[i*numProducers_]);

        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          triangleIntervals_[nid]
            = triangleIntervals_[nid - 1] + triangleIntervals_[nid];
        }
        triangleList_.reserve(triangleIntervals_.back()+10);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          for(auto &triangle : clusterTriangles_[nid]) {
            triangleList_.push_back(std::move(triangle));
          }
        }
        clusterTriangles_.clear(); clusterTriangles_.shrink_to_fit();
        printMsg("Triangles preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
      }

      return 0;
    }

    inline int preconditionTriangleEdgesInternal() override {
      relationVec_.push_back(RelationType::FERelation);
      return 0;
    }

    inline int preconditionTriangleLinksInternal() override {
      relationVec_.push_back(RelationType::FLRelation);
      return 0;
    }

    inline int preconditionTriangleStarsInternal() override {
      relationVec_.push_back(RelationType::FTRelation);
      return 0;
    }

    inline int preconditionVertexEdgesInternal() override {
      relationVec_.push_back(RelationType::VERelation);
      return 0;
    }

    inline int preconditionVertexLinksInternal() override {
      if(maxCellDim_ == 2) {
        preconditionEdgesInternal();
      } else if(maxCellDim_ == 3) {
        preconditionTrianglesInternal();
      } else {
        // unsupported dimension
        printErr("Unsupported dimension for vertex link precondition.");
        return -1;
      }
      relationVec_.push_back(RelationType::VLRelation);
      return 0;
    }

    inline int preconditionVertexNeighborsInternal() override {
      relationVec_.push_back(RelationType::VVRelation);
      return 0;
    }

    inline int preconditionVertexStarsInternal() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!cellArray_)
        return -1;
#endif

      relationVec_.push_back(RelationType::VTRelation);
      return 0;
    }

    inline int preconditionVertexTrianglesInternal() override {
      relationVec_.push_back(RelationType::VFRelation);
      return 0;
    }

    /**
     * Initialize the buffer size per thread with the given ratio.
     */
    inline void setBufferSize(const float ratio = 0.2f) {
      bufferSize_ = ratio / threadNumber_ * nodeNumber_ + 1;
      this->printMsg("Buffer capacity per thread: " + std::to_string(bufferSize_));
    }

    /** 
     * Change the worker mode.
     */
    inline void setWorkMode(const int mode) {
      // std::unique_lock<std::mutex> wlock(workerMutex_);
      #ifdef PRINT_WAITING_TIME
      missCnts_ = std::vector<int>(threadNumber_, 0);
      waitingTimes_ = std::vector<double>(threadNumber_, 0.0);
      #endif
      workMode_ = mode;
      sharedClusterIds_ = std::vector<SimplexId>(threadNumber_);
      workerRelations_ = std::vector<RelationType>(threadNumber_);
      if (workMode_ == 2) {
        workerRelationIds_ = std::vector<size_t>(threadNumber_);
      } else if (workMode_ == 3) {
        workerClusterIds_ = std::vector<size_t>(threadNumber_);
        workerClusterVecs_ = std::vector<std::vector<SimplexId>>(threadNumber_);
        workerRelationIds_ = std::vector<size_t>(threadNumber_);
      } else if (workMode_ == 4) {
        workerClusterIds_ = std::vector<size_t>(threadNumber_);
        workerClusterVecs_ = std::vector<std::vector<SimplexId>>(threadNumber_);
        workerRelations_ = std::vector<RelationType>(threadNumber_);
      }
      // wlock.unlock();
      this->printMsg("Set the worker mode to " + std::to_string(mode));
    }

    /**
     * Initialize the number of producers. 
     */
    inline void setProducerNumber(const int num) {
      if(producers_.empty()) {
        numProducers_ = num;
        int totalProducerNum_ = threadNumber_ * numProducers_;
        clustersPerThread_ = nodeNumber_ / (totalProducerNum_ + threadNumber_);
        semaphores_ = std::vector<sem_t>(totalProducerNum_);
        producers_ = std::vector<std::thread>(totalProducerNum_);
        finished_ = std::vector<SimplexId>(totalProducerNum_, 1);
        for (int i = 0; i < threadNumber_; i++) {
          producers_[i * numProducers_] = std::thread(&AcTopo::leaderProcedure, this, i);
          if(sem_init(&semaphores_[i*numProducers_], 0, 0)) {
            this->printErr("Cannot initialize the leader semaphore vector!");
          }
        }

        this->printMsg("Total number of producers: " + std::to_string(totalProducerNum_));
      }
      else {
        this->printErr("Cannot reset the number of producers!");
      }
    }

  protected:
    int clear();

    /**
     * Find the corresponding node index given the id.
     */
    inline SimplexId findNodeIndex(SimplexId id, int idType) const {
      const std::vector<SimplexId> *intervals = nullptr;
      // determine which vector to search
      if(idType == EDGE_ID) {
        intervals = &edgeIntervals_;
      } else if(idType == TRIANGLE_ID) {
        intervals = &triangleIntervals_;
      } else {
        return -1;
      }

      std::vector<SimplexId>::const_iterator low
        = lower_bound(intervals->begin(), intervals->end(), id);
      return (low - intervals->begin());
    }

    /**
     * Build the internal edge list in the cluster.
     */
    int buildInternalEdgeList(const SimplexId &clusterId, bool buildBoundary = false);

    /**
     * Build the internal triangle list in the cluster.
     */
    int buildInternalTriangleList(const SimplexId &clusterId, bool buildBoundary = false);

    /**
     * Build the boundary top cell list in the node. 
     */
    int buildBoundaryTopCellList(const SimplexId &clusterId) const; 

    /**
     * Get the cell edges for all cells in a given node.
     * Check if the tetraEdges_ is NULL before calling the function.
     */
    int getClusterCellEdges(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the cell neighbors for all cells in a given node.
     */
    int getClusterCellNeighbors(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the cell triangles for all cells in a given node.
     */
    int getClusterCellTriangles(const SimplexId &clusterId, const ThreadId &threadId) const;
    /**
     * Get the edge links for all the edges in a given node.
     */
    int getClusterEdgeLinks(const SimplexId &clusterId, const ThreadId &threadId) const;
    /**
     * Get the edge stars for all the edges in a given node.
     */
    int getClusterEdgeStars(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the edge triangles for all the edges in a given node.
     */
    int getClusterEdgeTriangles(const SimplexId &clusterId, const ThreadId &threadId) const;
    /**
     * Get the triangle edges for all the triangles in a given node.
     */
    int getClusterTriangleEdges(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the triangle links for all the triangles in a given node.
     */
    int getClusterTriangleLinks(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the triangle stars for all the triangles in a given node.
     */
    int getClusterTriangleStars(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the vertex edges for all the vertices in a given node.
     */
    int getClusterVertexEdges(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the vertex links for all the vertices in a given node.
     */
    int getClusterVertexLinks(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the vertex neighbors for all the vertices in a given node.
     */
    int getClusterVertexNeighbors(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the vertex stars for all the vertices in a given node.
     * The function is similar as getVertexCells().
     */
    int getClusterVertexStars(const SimplexId &clusterId, const ThreadId &threadId) const;
    int getClusterVertexStarsVector(const SimplexId &clusterId, 
                std::vector<std::vector<SimplexId>> &vertexStarVec) const;

    /**
     * Get the vertex triangles for all the vertices in a given node.
     */
    int getClusterVertexTriangles(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the boundary vertices in a given node.
     */
    int getClusterBoundaryVertices(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Get the boundary edges in a given node.
     */
    int getClusterBoundaryEdges(const SimplexId &clusterId, const ThreadId &threadId) const;

    /**
     * Compute the required topological relation with the given cluster id and
     * relation type.
     */
    void computeRelation(const SimplexId &clusterId,
                         const RelationType &relation,
                         const ThreadId &threadId);

    /**
     * Wait on a topolgoical relation to be computed.
     */
    inline void waitForRelation(const SimplexId &clusterId,  const RelationType &relationId, const int &consumerId = 0) const {
      #ifdef PRINT_WAITING_TIME
      double start = omp_get_wtime();
      std::unique_lock<std::mutex> llck(leaderMutexes_[consumerId]);
      reqClusters_[consumerId] = clusterId;
      reqRelations_[consumerId] = relationId;
      waitings_[consumerId] = 1;
      llck.unlock();
      leaderCondVars_[consumerId].notify_one();
      sem_wait(&semaphores_[consumerId * numProducers_]);
      double end = omp_get_wtime();
      waitings_[consumerId] = 0;
      waitingTimes_[consumerId] += (end - start);
      missCnts_[consumerId]++;
      #else 
      std::unique_lock<std::mutex> llck(leaderMutexes_[consumerId]);
      reqClusters_[consumerId] = clusterId;
      reqRelations_[consumerId] = relationId;
      waitings_[consumerId] = 1;
      llck.unlock();
      leaderCondVars_[consumerId].notify_one();
      sem_wait(&semaphores_[consumerId * numProducers_]);
      waitings_[consumerId] = 0;
      #endif
    }

    /**
     * The procedure run by the head prodcuer, which needs to:
     * - Communicate with the consumer thread;
     * - Compute the relation for the current request; 
     * - Coordinate other workers;
     * - Clean the memory space.
     */
    void leaderProcedure(const int &consumerId);

    /**
     * The procedure run by the worker producer.
     */
    void preconditionFunc(const int &workerId, const int &consumerId);
    void workerProcedure(const int &workerId, const int &consumerId);

    /**
     * Protected class variables.
     */
    bool doublePrecision_;
    int maxCellDim_;
    mutable SimplexId prevClusterId_;
    SimplexId cellNumber_, vertexNumber_, nodeNumber_;
    const void *pointSet_;
    const int *vertexIndices_;
    std::vector<SimplexId> vertexIntervals_;
    std::vector<SimplexId> edgeIntervals_;
    std::vector<SimplexId> triangleIntervals_;
    std::vector<SimplexId> cellIntervals_;
    std::shared_ptr<CellArray> cellArray_;
    std::vector<std::vector<SimplexId>> externalCells_;
    std::vector<std::vector<std::array<SimplexId, 2>>> clusterEdges_;
    std::vector<std::vector<std::array<SimplexId, 3>>> clusterTriangles_;
    mutable std::vector<std::mutex> clusterMutexes_;
    mutable std::vector<ImplicitCluster> allClusters_;

    // Buffer system
    size_t bufferSize_;
    std::vector<std::vector<SimplexId>> connectivity_;
    std::vector<std::list<SimplexId>> sbuffers_;
    std::vector<boost::unordered_set<SimplexId>> bufferSets_;
    std::vector<std::mutex> bufferMutexes_;

    // Multithreading support
    int workMode_;
    int numProducers_;
    int clustersPerThread_;
    std::vector<int> preconditions_;
    mutable std::vector<std::mutex> leaderMutexes_;
    mutable std::vector<std::thread> producers_;
    mutable std::vector<int> waitings_;
    mutable std::vector<sem_t> semaphores_;
    mutable std::vector<SimplexId> reqClusters_;
    mutable std::vector<RelationType> reqRelations_;
    mutable std::vector<std::condition_variable> leaderCondVars_;

    // related to worker producers
    std::vector<int> changed_;
    std::vector<std::condition_variable> workerCondVars_;
    std::vector<std::mutex> workerMutexes_;
    std::vector<SimplexId> sharedClusterIds_;
    std::vector<size_t> workerClusterIds_;
    std::vector<std::vector<SimplexId>> workerClusterVecs_;
    std::vector<RelationType> workerRelations_;
    std::vector<size_t> workerRelationIds_;
    std::vector<int> finished_;
    mutable std::vector<RelationType> relationVec_;
  };
} // namespace ttk
