#pragma once

// base code includes
#include <AbstractTriangulation.h>
#include <CellArray.h>
#include <FlatJaggedArray.h>
#include <bitset>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <list>
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

namespace ttk {

  class ImplicitCluster {
  private:
    /* components */
    SimplexId nid;
    /* boundary cells */
    std::vector<bool> boundaryVertices_;
    std::vector<bool> boundaryEdges_;
    std::vector<bool> boundaryTriangles_;
    /* vertex relations */
    FlatJaggedArray vertexEdges_;
    FlatJaggedArray vertexLinks_;
    FlatJaggedArray vertexNeighbors_;
    FlatJaggedArray vertexStars_;
    FlatJaggedArray vertexTriangles_;
    /* edge relations */
    // edgeVertex relation can be extracted from internal edge list
    FlatJaggedArray edgeLinks_;
    FlatJaggedArray edgeStars_;
    FlatJaggedArray edgeTriangles_;
    /* triangle relations */
    // triangleVertex relation can be extracted from internal triangle list
    std::vector<std::array<SimplexId, 3>> triangleEdges_;
    FlatJaggedArray triangleLinks_;
    FlatJaggedArray triangleStars_;
    /* cell relations */
    std::vector<std::array<SimplexId, 6>> tetraEdges_;
    FlatJaggedArray cellNeighbors_;
    std::vector<std::array<SimplexId, 4>> tetraTriangles_;

  public:
    ImplicitCluster() {
    }
    ImplicitCluster(SimplexId id) : nid(id) {
    }
    ~ImplicitCluster() {
    }
    inline void clear() {
      // keep the edge and triangle lists
      /* boundary cells */
      boundaryVertices_ = std::vector<bool>{};
      boundaryEdges_ = std::vector<bool>{};
      // keep boundary triangles 
      // boundaryTriangles_ = std::vector<bool>{};
      /* vertex relations */
      vertexEdges_ = FlatJaggedArray{};
      vertexLinks_ = FlatJaggedArray{};
      vertexNeighbors_ = FlatJaggedArray{};
      vertexStars_ = FlatJaggedArray{};
      vertexTriangles_ = FlatJaggedArray{};
      /* edge relations */
      // edgeVertex relation can be extracted from internal edge list
      edgeLinks_ = FlatJaggedArray{};
      edgeStars_ = FlatJaggedArray{};
      edgeTriangles_ = FlatJaggedArray{};
      /* triangle relations */
      // triangleVertex relation can be extracted from internal triangle list
      triangleEdges_ = std::vector<std::array<SimplexId, 3>>{};
      triangleLinks_ = FlatJaggedArray{};
      triangleStars_ = FlatJaggedArray{};
      /* cell relations */
      tetraEdges_ = std::vector<std::array<SimplexId, 6>>{};
      cellNeighbors_ = FlatJaggedArray{};
      tetraTriangles_ = std::vector<std::array<SimplexId, 4>>{};
    }

    friend class AcTopo;
  };

  class AcTopo final : public AbstractTriangulation {

  public:
    // Keep track of timings
    mutable double relationTime_, waitingTime_;
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].tetraEdges_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TERelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localEdgeId >= (int)(allClusters_[nid].tetraEdges_)[localCellId].size())
        return -2;
#endif
      edgeId = (allClusters_[nid].tetraEdges_)[localCellId][localEdgeId];

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
          if(allClusters_[nid].tetraEdges_.empty()) {
            getClusterCellEdges(nid);
          }
          for(size_t i = 0; i < allClusters_[nid].tetraEdges_.size(); i++) {
            cellEdgeVector_.push_back({allClusters_[nid].tetraEdges_[i].begin(),
                                       allClusters_[nid].tetraEdges_[i].end()});
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].cellNeighbors_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TTRelation);
      }
      
#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId >= allClusters_[nid].cellNeighbors_.size(localCellId))  
        return -2;
#endif
      neighborId = allClusters_[nid].cellNeighbors_.get(localCellId, localNeighborId);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].cellNeighbors_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TTRelation);
      }
      return allClusters_[nid].cellNeighbors_.size(localCellId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      cellNeighborList_.reserve(cellNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localCellNeighbors;
        if(allClusters_[nid].cellNeighbors_.empty()) {
          getClusterCellNeighbors(nid);
        }
        allClusters_[nid].cellNeighbors_.copyTo(localCellNeighbors);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].tetraTriangles_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::TFRelation);
      }

      triangleId = (allClusters_[nid].tetraTriangles_)[localCellId][localTriangleId];
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
          if(allClusters_[nid].tetraTriangles_.empty()) {
            getClusterCellTriangles(nid);
          }
          for(size_t i = 0; i < allClusters_[nid].tetraTriangles_.size(); i++) {
            cellTriangleVector_.push_back(
              {allClusters_[nid].tetraTriangles_[i].begin(),
               allClusters_[nid].tetraTriangles_[i].end()});
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

      // Timer t;
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ELRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId >= allClusters_[nid].edgeLinks_.size(localEdgeId))
        return -2;
#endif
      linkId = allClusters_[nid].edgeLinks_.get(localEdgeId, localLinkId);

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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ELRelation);
      }
      return allClusters_[nid].edgeLinks_.size(localEdgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      edgeLinkList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeLinks;
        if(allClusters_[nid].edgeLinks_.empty()) {
          getClusterEdgeLinks(nid);
        }
        allClusters_[nid].edgeLinks_.copyTo(localEdgeLinks);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ETRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId >= allClusters_[nid].edgeStars_.size(localEdgeId))
        return -2;
#endif
      starId = allClusters_[nid].edgeStars_.get(localEdgeId, localStarId);

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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::ETRelation);
      }
      SimplexId number = allClusters_[nid].edgeStars_.size(localEdgeId);

      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      edgeStarList_.reserve(edgeIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeStars;
        if(allClusters_[nid].edgeStars_.empty()) {
          getClusterEdgeStars(nid);
        }
        allClusters_[nid].edgeStars_.copyTo(localEdgeStars);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeTriangles_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::EFRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localTriangleId >= allClusters_[nid].edgeTriangles_.size(localEdgeId))
        return -2;
#endif
      triangleId = allClusters_[nid].edgeTriangles_.get(localEdgeId, localTriangleId);

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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].edgeTriangles_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::EFRelation);
      }

      SimplexId number = allClusters_[nid].edgeTriangles_.size(localEdgeId);
      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override {
      edgeTriangleList_.reserve(edgeIntervals_.back() + 1);
      std::bitset<5> type(2);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localEdgeTriangles;
        if(allClusters_[nid].edgeTriangles_.empty()) {
          getClusterEdgeTriangles(nid);
        }
        allClusters_[nid].edgeTriangles_.copyTo(localEdgeTriangles);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleEdges_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FERelation);
      }

      edgeId = (allClusters_[nid].triangleEdges_)[localTriangleId][localEdgeId];
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
          if(allClusters_[nid].triangleEdges_.empty()) {
            getClusterTriangleEdges(nid);
          }
          for(size_t i = 0; i < allClusters_[nid].triangleEdges_.size(); i++) {
            triangleEdgeVector_.push_back(
              {allClusters_[nid].triangleEdges_[i].begin(),
               allClusters_[nid].triangleEdges_[i].end()});
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FLRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId >= allClusters_[nid].triangleLinks_.size(localTriangleId))
        return -2;
#endif
      linkId = allClusters_[nid].triangleLinks_.get(localTriangleId, localLinkId);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FLRelation);
      }

      return allClusters_[nid].triangleLinks_.size(localTriangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      triangleLinkList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localTriangleLinks;
        if(allClusters_[nid].triangleLinks_.empty()) {
          getClusterTriangleLinks(nid);
        }
        allClusters_[nid].triangleLinks_.copyTo(localTriangleLinks);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FTRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId >= allClusters_[nid].triangleStars_.size(localTriangleId))
        return -2;
#endif
      starId = allClusters_[nid].triangleStars_.get(localTriangleId, localStarId);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].triangleStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::FTRelation);
      }
      
      return allClusters_[nid].triangleStars_.size(localTriangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
      triangleStarList_.reserve(triangleIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localTriangleStars;
        if(allClusters_[nid].triangleStars_.empty()) {
          getClusterTriangleStars(nid);
        }
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

      // Timer t;
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexEdges_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VERelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localEdgeId >= allClusters_[nid].vertexEdges_.size(localVertexId))
        return -2;
#endif
      edgeId = allClusters_[nid].vertexEdges_.get(localVertexId, localEdgeId);

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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexEdges_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VERelation);
      }
      SimplexId number = allClusters_[nid].vertexEdges_.size(localVertexId);

      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override {
      vertexEdgeList_.reserve(vertexNumber_);
      std::bitset<5> type(1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexEdges;
        if(allClusters_[nid].vertexEdges_.empty()) {
          getClusterVertexEdges(nid);
        }
        allClusters_[nid].vertexEdges_.copyTo(localVertexEdges);
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

      // Timer t;
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VLRelation);
      }

      if(localLinkId >= allClusters_[nid].vertexLinks_.size(localVertexId)) {
        linkId = -2;
      } else {
        linkId = allClusters_[nid].vertexLinks_.get(localVertexId, localLinkId);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexLinks_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VLRelation);
      }

      SimplexId number = allClusters_[nid].vertexLinks_.size(localVertexId);
      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      vertexLinkList_.reserve(vertexIntervals_.back() + 1);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexLinks;
        if(allClusters_[nid].vertexLinks_.empty()) {
          getClusterVertexLinks(nid);
        }
        allClusters_[nid].vertexLinks_.copyTo(localVertexLinks);
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

      // Timer t;
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexNeighbors_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VVRelation);
      }
      
#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId >= allClusters_[nid].vertexNeighbors_.size(localVertexId)) {
        return -2;
#endif
      neighborId
        = allClusters_[nid].vertexNeighbors_.get(localVertexId, localNeighborId);

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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexNeighbors_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VVRelation);
      }

      SimplexId number = allClusters_[nid].vertexNeighbors_.size(localVertexId);
      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
      vertexNeighborList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexNeighbors;
        if(allClusters_[nid].vertexNeighbors_.empty()) {
          getClusterVertexNeighbors(nid);
        }
        allClusters_[nid].vertexNeighbors_.copyTo(localVertexNeighbors);
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

      // Timer t;
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VTRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId >= allClusters_[nid].vertexStars_.size(localVertexId))
        return -2;
#endif
      starId = allClusters_[nid].vertexStars_.get(localVertexId, localStarId);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexStars_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VTRelation);
      }

      SimplexId number = allClusters_[nid].vertexStars_.size(localVertexId);
      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      vertexStarList_.reserve(vertexNumber_);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexStars;
        if(allClusters_[nid].vertexStars_.empty()) {
          getClusterVertexStars(nid);
        }
        allClusters_[nid].vertexStars_.copyTo(localVertexStars);
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

      // Timer t;
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      if(allClusters_[nid].vertexTriangles_.empty()) {
        // getClusterVertexTriangles(nid);
        waitForRelation(nid, RelationType::VFRelation);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(localTriangleId >= allClusters_[nid].vertexTriangles_.size(localVertexId)) {
        return -2;
#endif
      triangleId
        = allClusters_[nid].vertexTriangles_.get(localVertexId, localTriangleId);

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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].vertexTriangles_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::VFRelation);
      }

      SimplexId number = allClusters_[nid].vertexTriangles_.size(localVertexId);

      return number;
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override {
      vertexTriangleList_.reserve(vertexNumber_);
      std::bitset<5> type(2);
      for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
        std::vector<std::vector<SimplexId>> localVertexTriangles;
        if(allClusters_[nid].vertexTriangles_.empty()) {
          getClusterVertexTriangles(nid);
        }
        allClusters_[nid].vertexTriangles_.copyTo(localVertexTriangles);
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

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].boundaryEdges_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::BERelation);
      }

      return (allClusters_[nid].boundaryEdges_)[localEdgeId];
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

      return (allClusters_[nid].boundaryTriangles_)[localTriangleId];
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return false;
#endif
      SimplexId nid = vertexIndices_[vertexId];
      SimplexId localVertexId = vertexId - vertexIntervals_[nid - 1] - 1;

      std::unique_lock<std::mutex> clck(clusterMutexes_[nid]);
      if(allClusters_[nid].boundaryVertices_.empty()) {
        clck.unlock();
        waitForRelation(nid, RelationType::BVRelation);
      }

      return (allClusters_[nid].boundaryVertices_)[localVertexId];
    }

    inline int preconditionBoundaryEdgesInternal() override {
      if(maxCellDim_ == 2) {
        Timer t;
        currRelation_ = RelationType::IBCList;
        std::unique_lock<std::mutex> llck(leaderMutex_);
        precondition_ = true;
        llck.unlock();
        leaderCondVar_.notify_one();
        sem_wait(&semaphores_[0]);
        printMsg("Boundary edges preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
      } else if(maxCellDim_ == 3) {
        preconditionBoundaryTriangles();
        relationVec_.push_back(RelationType::BERelation);
      } else {
        // unsupported dimension
        printErr(
          "[AcTopo] Unsupported dimension for boundary precondtion.");
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() override {
      if(maxCellDim_ == 2) {
        return 0;
      } else if(maxCellDim_ == 3) {
        Timer t;
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
        // use the producer threads to compute the edges
        std::unique_lock<std::mutex> llck(leaderMutex_);
        currRelation_ = RelationType::IBCList;
        precondition_ = true;
        llck.unlock();
        leaderCondVar_.notify_one();
        sem_wait(&semaphores_[0]);
        printMsg("Boundary triangles preconditioned in " + std::to_string(t.getElapsedTime()) + "s.");
      } else {
        // unsupported dimension
        printErr(
          "[AcTopo] Unsupported dimension for boundary precondtion.");
        return -1;
      }
      return 0;
    }

    inline int preconditionBoundaryVerticesInternal() override {
      if(maxCellDim_ == 2) {
        preconditionBoundaryEdges();
      } else if(maxCellDim_ == 3) {
        preconditionBoundaryTriangles();
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
        std::unique_lock<std::mutex> llck(leaderMutex_);
        currRelation_ = RelationType::IEList;
        precondition_ = true;
        llck.unlock();
        leaderCondVar_.notify_one();
        sem_wait(&semaphores_[0]);

        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          edgeIntervals_[nid] = edgeIntervals_[nid - 1] + edgeIntervals_[nid];
        }
        edgeList_.reserve(edgeIntervals_.back()+1);
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
        printErr("Unsupported dimension for edge link precondtion.");
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
        std::unique_lock<std::mutex> llck(leaderMutex_);
        currRelation_ = RelationType::IFList;
        precondition_ = true;
        llck.unlock();
        leaderCondVar_.notify_one();
        sem_wait(&semaphores_[0]);

        for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
          triangleIntervals_[nid]
            = triangleIntervals_[nid - 1] + triangleIntervals_[nid];
        }

        triangleList_.reserve(triangleIntervals_.back()+1);
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
        printErr("Unsupported dimension for vertex link precondtion.");
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
     * Initialize the buffer with the ratio.
     */
    inline void setBufferSize(const float ratio = 0.2f) {
      bufferSize_ = nodeNumber_ * ratio + 1;
      this->printMsg("Buffer capacity: " + std::to_string(bufferSize_));
    }

    /** 
     * Change the worker mode.
     */
    inline void setWorkMode(const int mode) {
      // std::unique_lock<std::mutex> wlock(workerMutex_);
      workMode_ = mode;
      // wlock.unlock();
      this->printMsg("Set the worker mode to " + std::to_string(mode));
    }

    /**
     * Initialize the number of producers. 
     */
    inline void setProducerNumber(const int num) {
      if(producers_.empty()) {
        numProducers_ = num;
        semaphores_ = std::vector<sem_t>(numProducers_+1);
        producers_ = std::vector<std::thread>(numProducers_);
        producers_[0] = std::thread(&AcTopo::leaderProcedure, this);

        for(int i = 0; i <= numProducers_; i++) {
          if(sem_init(&semaphores_[i], 0, 0)) {
            this->printErr("Cannot initialize the semaphore vector!");
          }
        }

        this->printMsg("Number of producers: " + std::to_string(numProducers_));
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
     * Build the internal edge list in the node.
     */
    int buildInternalEdgeList(const SimplexId &clusterId);

    /**
     * Build the internal triangle list in the node.
     */
    int buildInternalTriangleList(const SimplexId &clusterId);

    /**
     * Build the boundary top cell list in the node. 
     */
    int buildBoundaryTopCellList(const SimplexId &clusterId) const; 

    /**
     * Get the cell edges for all cells in a given node.
     * Check if the tetraEdges_ is NULL before calling the function.
     */
    int getClusterCellEdges(const SimplexId &clusterId) const;

    /**
     * Get the cell neighbors for all cells in a given node.
     */
    int getClusterCellNeighbors(const SimplexId &clusterId) const;

    /**
     * Get the cell triangles for all cells in a given node.
     */
    int getClusterCellTriangles(const SimplexId &clusterId) const;
    /**
     * Get the edge links for all the edges in a given node.
     */
    int getClusterEdgeLinks(const SimplexId &clusterId) const;
    /**
     * Get the edge stars for all the edges in a given node.
     */
    int getClusterEdgeStars(const SimplexId &clusterId) const;

    /**
     * Get the edge triangles for all the edges in a given node.
     */
    int getClusterEdgeTriangles(const SimplexId &clusterId) const;
    /**
     * Get the triangle edges for all the triangles in a given node.
     */
    int getClusterTriangleEdges(const SimplexId &clusterId) const;

    /**
     * Get the triangle links for all the triangles in a given node.
     */
    int getClusterTriangleLinks(const SimplexId &clusterId) const;

    /**
     * Get the triangle stars for all the triangles in a given node.
     */
    int getClusterTriangleStars(const SimplexId &clusterId) const;

    /**
     * Get the vertex edges for all the vertices in a given node.
     */
    int getClusterVertexEdges(const SimplexId &clusterId) const;

    /**
     * Get the vertex links for all the vertices in a given node.
     */
    int getClusterVertexLinks(const SimplexId &clusterId) const;

    /**
     * Get the vertex neighbors for all the vertices in a given node.
     */
    int getClusterVertexNeighbors(const SimplexId &clusterId) const;

    /**
     * Get the vertex stars for all the vertices in a given node.
     * The function is similar as getVertexCells().
     */
    int getClusterVertexStars(const SimplexId &clusterId) const;
    int getClusterVertexStarsVector(const SimplexId &clusterId, 
                std::vector<std::vector<SimplexId>> &vertexStarVec) const;

    /**
     * Get the vertex triangles for all the vertices in a given node.
     */
    int getClusterVertexTriangles(const SimplexId &clusterId) const;

    /**
     * Get the boundary vertices in a given node.
     */
    int getClusterBoundaryVertices(const SimplexId &clusterId) const;

    /**
     * Get the boundary edges in a given node.
     */
    int getClusterBoundaryEdges(const SimplexId &clusterId) const;

    /**
     * Compute the required topological relation with the given cluster id and
     * relation type.
     */
    void computeRelation(const SimplexId &clusterId,
                         const RelationType &relation);

    /**
     * Wait on a topolgoical relation to be computed.
     */
    inline void waitForRelation(const SimplexId &clusterId, const RelationType &relationId) const {
      #ifdef TTK_ENABLE_OPENMP
      double start = omp_get_wtime();
      #else 
      Timer t;
      #endif
      std::unique_lock<std::mutex> llck(leaderMutex_);
      currCluster_ = clusterId;
      currRelation_ = relationId;
      waiting_ = true;
      llck.unlock();
      leaderCondVar_.notify_one();
      sem_wait(&semaphores_[0]);
      #ifdef TTK_ENABLE_OPENMP
      double end = omp_get_wtime();
      waitingTime_ += (end - start);
      #else 
      waitingTime_ += t.getElapsedTime();
      #endif
      missCnt_++;
    }

    /**
     * The procedure run by the head prodcuer, which needs to:
     * - Communicate with the consumer thread;
     * - Compute the relation for the current request; 
     * - Coordinate other workers;
     * - Clean the memory space.
     */
    void leaderProcedure();

    /**
     * The procedure run by the worker producer.
     */
    void preconditionFunc(const int &workerId);
    void workerProcedure(const int &workerId);

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
    std::mutex bufferMutex_;
    std::vector<std::vector<SimplexId>> connectivity_;
    std::list<SimplexId> sbuffer_;
    boost::unordered_set<SimplexId> bufferSet_;
    mutable size_t missCnt_;

    // Multithreading support
    int workMode_;
    int numProducers_;
    bool precondition_;
    mutable bool waiting_;
    mutable std::mutex leaderMutex_;
    mutable std::vector<std::thread> producers_;
    mutable std::vector<sem_t> semaphores_;
    mutable SimplexId currCluster_;
    mutable RelationType currRelation_;
    mutable std::condition_variable leaderCondVar_;

    // related to worker producers
    bool changed_;
    std::mutex workerMutex_;
    SimplexId sharedClusterId_;
    size_t workerClusterIdx_;
    size_t workerRelationIdx_;
    RelationType workerRelation_;
    std::vector<int> finished_;
    std::condition_variable workerCondVar_;
    std::vector<SimplexId> workerClusterVec_;
    mutable std::vector<RelationType> relationVec_;
  };
} // namespace ttk
