#include <AcTopo.h>
#include <set>
#include <algorithm>

#define findEdgeIdx(x, e) std::lower_bound(edgeList_.begin()+edgeIntervals_[x-1]+1, edgeList_.begin()+edgeIntervals_[x]+1, e) - edgeList_.begin() - edgeIntervals_[x-1] - 1
#define findTriangleIdx(x, f) std::lower_bound(triangleList_.begin()+triangleIntervals_[x-1]+1, triangleList_.begin()+triangleIntervals_[x]+1, f) - triangleList_.begin() - triangleIntervals_[x-1] - 1

using namespace ttk;

AcTopo::AcTopo() {
  // set output files
  // consumerFile_.open("consumer.txt", std::ofstream::trunc);
  // if(consumerFile_.fail()) {
  //   this->printErr("Cannot open file for consumer!");
  // }
  // producerFile_.open("producer.txt", std::ofstream::trunc);
  // if(producerFile_.fail()) {
  //   this->printErr("Cannot open file for producer!");
  // }
  // logFile_.open("log_topocluster.txt", std::ofstream::trunc);
  // if(logFile_.fail()) {
  //   this->printErr("Cannot open file for logging!");
  // }
  setDebugMsgPrefix("AcTopo");
  clear();
}

AcTopo::AcTopo(const AcTopo &rhs)
  : doublePrecision_(rhs.doublePrecision_), maxCellDim_(rhs.maxCellDim_),
    cellNumber_(rhs.cellNumber_), vertexNumber_(rhs.vertexNumber_),
    nodeNumber_(rhs.nodeNumber_), pointSet_(rhs.pointSet_),
    vertexIndices_(rhs.vertexIndices_), vertexIntervals_(rhs.vertexIntervals_),
    edgeIntervals_(rhs.edgeIntervals_),
    triangleIntervals_(rhs.triangleIntervals_),
    cellIntervals_(rhs.cellIntervals_), cellArray_(rhs.cellArray_),
    externalCells_(rhs.externalCells_) {
  // mutex and conditional variable are not copied
}

AcTopo::~AcTopo() {
  if(!producers_.empty()) {
    for(int i = 0; i < threadNumber_; i++) {
      std::unique_lock<std::mutex> llck(leaderMutexes_[i]);
      reqClusters_[i] = -1;
      waitings_[i] = 1;
      llck.unlock();
      leaderCondVars_[i].notify_one();
    }
    for(size_t j = 0; j < producers_.size(); j++) {
      if(producers_[j].joinable()) {
        producers_[j].join();
      }
    }
    #ifdef PRINT_WAITING_TIME
    for(int i = 0; i < threadNumber_; i++) {
      std::cout << "consumer " << i << ": " << waitingTimes_[i] << ", " << missCnts_[i] << std::endl;
    }
    #endif
  }
  // consumerFile_.close();
  // producerFile_.close();
}

AcTopo &AcTopo::operator=(const AcTopo &rhs) {
  if(this != &rhs) {
    doublePrecision_ = rhs.doublePrecision_;
    maxCellDim_ = rhs.maxCellDim_;
    cellNumber_ = rhs.cellNumber_;
    vertexNumber_ = rhs.vertexNumber_;
    nodeNumber_ = rhs.nodeNumber_;
    pointSet_ = rhs.pointSet_;
    vertexIndices_ = rhs.vertexIndices_;
    vertexIntervals_ = rhs.vertexIntervals_;
    edgeIntervals_ = rhs.edgeIntervals_;
    triangleIntervals_ = rhs.triangleIntervals_;
    cellIntervals_ = rhs.cellIntervals_;
    cellArray_ = rhs.cellArray_;
    externalCells_ = rhs.externalCells_;
    // mutex and conditional varaible are not assigned
  }
  return *this;
}

int AcTopo::clear() {
  vertexNumber_ = 0;
  cellNumber_ = 0;
  nodeNumber_ = 0;
  doublePrecision_ = false;
  bufferSize_ = 0.2;
  workMode_ = 1;
  prevClusterId_ = 0;
  numProducers_ = 1;
  producers_.clear();
  preconditions_ = std::vector<int>(threadNumber_, 0);
  sbuffers_ = std::vector<std::list<SimplexId>>(threadNumber_);
  bufferSets_ = std::vector<boost::unordered_set<SimplexId>>(threadNumber_);
  bufferMutexes_ = std::vector<std::mutex>(threadNumber_);
  leaderMutexes_ = std::vector<std::mutex>(threadNumber_);
  leaderCondVars_ = std::vector<std::condition_variable>(threadNumber_);
  workerCondVars_ = std::vector<std::condition_variable>(threadNumber_);
  changed_ = std::vector<int>(threadNumber_, 0);
  workerMutexes_ = std::vector<std::mutex>(threadNumber_);
  waitings_ = std::vector<int>(threadNumber_, 0);
  reqClusters_ = std::vector<SimplexId>(threadNumber_);
  reqRelations_ = std::vector<RelationType>(threadNumber_);
  return AbstractTriangulation::clear();
}

#ifdef TTK_CELL_ARRAY_NEW
// Layout with connectivity + offset array (new)
int AcTopo::setInputCells(const SimplexId &cellNumber,
                          const LongSimplexId *connectivity,
                          const LongSimplexId *offset) {

  // Cell Check
  {
    if(cellNumber > 0) {
      const auto &cellDimension = offset[1] - offset[0] - 1;

      if(cellDimension < 0 || cellDimension > 3) {
        this->printErr("Unable to create triangulation for cells of "
                        "dimension 4 or higher ("
                        + std::to_string(cellDimension) + ").");
        return -1;
      }

      bool error = false;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(SimplexId i = 0; i < cellNumber; i++) {
        if(offset[i + 1] - offset[i] - 1 != cellDimension) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
          error = true;
        }
      }

      if(error) {
        this->printErr("Unable to create triangulation for "
                        "inhomogeneous\ncell dimensions.");
        return -2;
      }
    }
  }

  if(cellNumber_)
    clear();

  cellNumber_ = cellNumber;
  std::vector<SimplexId> vertexMap(vertexNumber_);
  reorderVertices(vertexMap);
  reorderCells(vertexMap, cellNumber, connectivity, offset);
  cellArray_
    = std::make_shared<CellArray>(connectivity, offset, cellNumber);

  // ASSUME Regular Mesh Here to compute dimension!
  if(cellNumber) {
    if(cellArray_->getCellVertexNumber(0) == 3) {
      maxCellDim_ = 2;
    } else {
      maxCellDim_ = 3;
    }
  }

  // initialize all cluster pointers
  allClusters_ = std::vector<ImplicitCluster>(nodeNumber_ + 1);
  clusterMutexes_ = std::vector<std::mutex>(nodeNumber_ + 1);
  for(int i = 1; i <= nodeNumber_; i++) {
    allClusters_[i] = ImplicitCluster(i, threadNumber_);
  }

  // build the connectivity vector
  // add the self cluster in the adjacency list
  std::vector<std::set<SimplexId>> adjacency(nodeNumber_ + 1);

  // set up the adjacency list
  for(SimplexId i = 0; i < cellNumber_; i++) {
    std::vector<SimplexId> vertexClusters(maxCellDim_ + 1);
    for(SimplexId j = 0; j <= maxCellDim_; j++) {
      vertexClusters[j]
        = vertexIndices_[cellArray_->getCellVertex(i, j)];
    }
    for(SimplexId j = 0; j < maxCellDim_; j++) {
      for(SimplexId k = j + 1; k <= maxCellDim_; k++) {
        if(vertexClusters[j] != vertexClusters[k]) {
          adjacency[vertexClusters[j]].insert(vertexClusters[k]);
          adjacency[vertexClusters[k]].insert(vertexClusters[j]);
        }
      }
    }
  }

  connectivity_ = std::vector<std::vector<SimplexId>>(nodeNumber_ + 1);
  for(SimplexId i = 1; i <= nodeNumber_; i++) {
    connectivity_[i] = {adjacency[i].begin(), adjacency[i].end()};
  }

  return 0;
}
#else
// Flat layout with a single array (legacy & default one)
int AcTopo::setInputCells(const SimplexId &cellNumber,
                          const LongSimplexId *cellArray) {
  if(cellNumber_)
    clear();

  cellNumber_ = cellNumber;

  if(cellNumber) {
    // assume regular mesh here to compute dimension
    maxCellDim_ = cellArray[0] - 1;
    std::vector<SimplexId> vertexMap(vertexNumber_);
    reorderVertices(vertexMap);
    reorderCells(vertexMap, cellArray);
    cellArray_ = std::make_shared<CellArray>(
      cellArray, cellNumber, cellArray[0] - 1);
  }
  // ASSUME Regular Mesh Here to compute dimension!
  if(cellNumber) {
    if(cellArray_->getCellVertexNumber(0) == 3) {
      maxCellDim_ = 2;
    } else {
      maxCellDim_ = 3;
    }
  }

  // initialize all cluster pointers
  allClusters_ = std::vector<ImplicitCluster>(nodeNumber_ + 1);
  clusterMutexes_ = std::vector<std::mutex>(nodeNumber_ + 1);
  for(int i = 1; i <= nodeNumber_; i++) {
    allClusters_[i] = ImplicitCluster(i);
  }

  // build the connectivity vector
  // add the self cluster in the adjacency list
  std::vector<std::set<SimplexId>> adjacency(nodeNumber_ + 1);

  // set up the adjacency list
  for(SimplexId i = 0; i < cellNumber_; i++) {
    std::vector<SimplexId> vertexClusters(maxCellDim_ + 1);
    for(SimplexId j = 0; j <= maxCellDim_; j++) {
      vertexClusters[j]
        = vertexIndices_[cellArray_->getCellVertex(i, j)];
    }
    for(SimplexId j = 0; j < maxCellDim_; j++) {
      for(SimplexId k = j + 1; k <= maxCellDim_; k++) {
        if(vertexClusters[j] != vertexClusters[k]) {
          adjacency[vertexClusters[j]].insert(vertexClusters[k]);
          adjacency[vertexClusters[k]].insert(vertexClusters[j]);
        }
      }
    }
  }

  connectivity_ = std::vector<std::vector<SimplexId>>(nodeNumber_ + 1);
  for(SimplexId i = 1; i <= nodeNumber_; i++) {
    connectivity_[i] = {adjacency[i].begin(), adjacency[i].end()};
  }

  return 0;
}
#endif

int AcTopo::reorderVertices(std::vector<SimplexId> &vertexMap) {
  // get the number of nodes (the max value in the array)
  for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
    if(vertexIndices_[vid] > nodeNumber_) {
      nodeNumber_ = vertexIndices_[vid];
    }
  }
  nodeNumber_++; // since the index starts from 0
  std::vector<std::vector<SimplexId>> nodeVertices(nodeNumber_);
  for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
    nodeVertices[vertexIndices_[vid]].push_back(vid);
  }

  // update the vertex intervals
  vertexIntervals_.resize(nodeNumber_ + 1);
  vertexIntervals_[0] = -1;
  SimplexId vertexCount = 0;
  for(SimplexId nid = 0; nid < nodeNumber_; nid++) {
    for(SimplexId vid : nodeVertices[nid]) {
      vertexMap[vid] = vertexCount++;
    }
    vertexIntervals_[nid + 1] = vertexCount - 1;
  }

  // rearange the vertex coordinate values
  if(doublePrecision_) {
    double *newPointSet = new double[3 * vertexNumber_];
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        newPointSet[3 * vertexMap[vid] + j]
          = ((double *)pointSet_)[3 * vid + j];
      }
    }
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        ((double *)pointSet_)[3 * vid + j] = newPointSet[3 * vid + j];
      }
    }
    delete[] newPointSet;
  } else {
    float *newPointSet = new float[3 * vertexNumber_];
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        newPointSet[3 * vertexMap[vid] + j]
          = ((float *)pointSet_)[3 * vid + j];
      }
    }
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        ((float *)pointSet_)[3 * vid + j] = newPointSet[3 * vid + j];
      }
    }
    delete[] newPointSet;
  }

  // change the vertex indices
  for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
    for(SimplexId vid = vertexIntervals_[nid - 1] + 1;
        vid <= vertexIntervals_[nid]; vid++) {
      ((int *)vertexIndices_)[vid] = nid;
    }
  }

  return 0;
}

int AcTopo::reorderCells(const std::vector<SimplexId> &vertexMap,
                  const SimplexId &cellNumber,
                  const LongSimplexId *connectivity,
                  const LongSimplexId *offset) {
  // change the indices in cell array
  SimplexId cellCount = 0, verticesPerCell = offset[1] - offset[0];
  std::vector<std::vector<SimplexId>> nodeCells(nodeNumber_ + 1);
  LongSimplexId *cellArr = ((LongSimplexId *)connectivity);

  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    SimplexId cellId = offset[cid];
    for(int j = 0; j < verticesPerCell; j++) {
      cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
    }
    std::sort(cellArr + cellId, cellArr + cellId + verticesPerCell);
    nodeCells[vertexIndices_[cellArr[cellId]]].push_back(cid);
  }

  // rearange the cell std::array
  cellIntervals_.resize(nodeNumber_ + 1);
  externalCells_.resize(nodeNumber_ + 1);
  cellIntervals_[0] = -1;
  LongSimplexId *newCellArray
    = new LongSimplexId[verticesPerCell * cellNumber];
  for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
    for(SimplexId cid : nodeCells[nid]) {
      SimplexId cellId = verticesPerCell * cid;
      SimplexId newCellId = verticesPerCell * cellCount;
      for(int j = 0; j < verticesPerCell; j++) {
        newCellArray[newCellId + j] = connectivity[cellId + j];
        if(newCellArray[newCellId + j] > vertexIntervals_[nid]) {
          SimplexId nodeNum = vertexIndices_[newCellArray[newCellId + j]];
          if(externalCells_[nodeNum].empty()
              || externalCells_[nodeNum].back() != cid) {
            externalCells_[nodeNum].push_back(cid);
          }
        }
      }
      cellCount++;
    }
    cellIntervals_[nid] = cellCount - 1;
  }

  // copy the new cell array back to original one
  for(SimplexId i = 0; i < verticesPerCell * cellNumber; i++) {
    ((LongSimplexId *)connectivity)[i] = newCellArray[i];
  }
  delete[] newCellArray;

  return 0;
}

int AcTopo::reorderCells(const std::vector<SimplexId> &vertexMap,
                  const LongSimplexId *cellArray) {
  // change the indices in cell array
  SimplexId cellCount = 0, verticesPerCell = cellArray[0];
  std::vector<std::vector<SimplexId>> nodeCells(nodeNumber_ + 1);
  LongSimplexId *cellArr = ((LongSimplexId *)cellArray);

  for(SimplexId cid = 0; cid < cellNumber_; cid++) {
    SimplexId cellId = (verticesPerCell + 1) * cid;
    for(int j = 1; j <= verticesPerCell; j++) {
      cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
    }
    std::sort(cellArr + cellId + 1, cellArr + cellId + 1 + verticesPerCell);
    nodeCells[vertexIndices_[cellArr[cellId + 1]]].push_back(cid);
  }

  // rearange the cell array
  cellIntervals_.resize(nodeNumber_ + 1);
  externalCells_.resize(nodeNumber_ + 1);
  cellIntervals_[0] = -1;
  LongSimplexId *newCellArray
    = new LongSimplexId[(verticesPerCell + 1) * cellNumber_];
  for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
    for(SimplexId cid : nodeCells[nid]) {
      SimplexId cellId = (verticesPerCell + 1) * cid;
      SimplexId newCellId = (verticesPerCell + 1) * cellCount;
      newCellArray[newCellId] = verticesPerCell;
      for(int j = 1; j <= verticesPerCell; j++) {
        newCellArray[newCellId + j] = cellArray[cellId + j];
        if(newCellArray[newCellId + j] > vertexIntervals_[nid]) {
          SimplexId nodeNum = vertexIndices_[newCellArray[newCellId + j]];
          if(externalCells_[nodeNum].empty()
              || externalCells_[nodeNum].back() != cid) {
            externalCells_[nodeNum].push_back(cid);
          }
        }
      }
      cellCount++;
    }
    cellIntervals_[nid] = cellCount - 1;
  }

  // copy the new cell array back to original one
  for(SimplexId i = 0; i < (verticesPerCell + 1) * cellNumber_; i++) {
    ((LongSimplexId *)cellArray)[i] = newCellArray[i];
  }
  delete[] newCellArray;

  return 0;
}

int AcTopo::buildInternalEdgeList(const SimplexId &clusterId, bool buildBoundary) {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId edgeCount = 0;
  // if we need to build the boundary edge list
  if(buildBoundary) {
    std::map<std::array<SimplexId, 2>, int> localEdgeMap;

    // loop through the internal cell list
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 2> edgeIds;

      // loop through each edge of the cell
      for(SimplexId j = 0; j < maxCellDim_; j++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        // the edge does not belong to the current node
        if(edgeIds[0] > vertexIntervals_[clusterId]) {
          break;
        }
        for(SimplexId k = j + 1; k < maxCellDim_ + 1; k++) {
          edgeIds[1] = cellArray_->getCellVertex(cid, k);
          localEdgeMap[edgeIds]++;
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 2> edgeIds;

      // loop through each edge of the cell from cell[1]
      for(SimplexId j = 1; j < maxCellDim_; j++) {
        for(SimplexId k = j + 1; k < maxCellDim_ + 1; k++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          edgeIds[1] = cellArray_->getCellVertex(cid, k);

          // the edge is in the current node
          if(edgeIds[0] > vertexIntervals_[clusterId - 1]
              && edgeIds[0] <= vertexIntervals_[clusterId]) {
            localEdgeMap[edgeIds]++;
          }
        }
      }
    }

    // compute boundary edges
    edgeCount = localEdgeMap.size();
    std::vector<bool> localBoundEdges(edgeCount, false);
    int k = 0;
    for(auto &m : localEdgeMap) {
      clusterEdges_[clusterId].push_back(std::move(m.first));
      if(m.second == 1) {
        localBoundEdges[k] = true;
      }
      k++;
    }
    allClusters_[clusterId].boundaryEdges_[0] = std::move(localBoundEdges);
  } 
  // otherwise, we just build the edge list
  else {
    std::set<std::array<SimplexId, 2>> localEdgeSet;

    // loop through the internal cell list
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 2> edgeIds;

      // loop through each edge of the cell
      for(SimplexId j = 0; j < maxCellDim_; j++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        // the edge does not belong to the current node
        if(edgeIds[0] > vertexIntervals_[clusterId]) {
          break;
        }
        for(SimplexId k = j + 1; k < maxCellDim_ + 1; k++) {
          edgeIds[1] = cellArray_->getCellVertex(cid, k);

          // not found in the edge map - assign new edge id
          if(localEdgeSet.find(edgeIds) == localEdgeSet.end()) {
            edgeCount++;
            localEdgeSet.insert(edgeIds);
          }
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 2> edgeIds;

      // loop through each edge of the cell from cell[1]
      for(SimplexId j = 1; j < maxCellDim_; j++) {
        for(SimplexId k = j + 1; k < maxCellDim_ + 1; k++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          edgeIds[1] = cellArray_->getCellVertex(cid, k);

          // the edge is in the current node
          if(edgeIds[0] > vertexIntervals_[clusterId - 1]
              && edgeIds[0] <= vertexIntervals_[clusterId]) {
            if(localEdgeSet.find(edgeIds) == localEdgeSet.end()) {
              edgeCount++;
              localEdgeSet.insert(edgeIds);
            }
          }
        }
      }
    }

    clusterEdges_[clusterId] = {localEdgeSet.begin(), localEdgeSet.end()};
  }
  return edgeCount;
}

int AcTopo::buildInternalTriangleList(const SimplexId &clusterId, bool buildBoundary) {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId triangleCount = 0;

  // if we need to build the boundary triangle list
  if(buildBoundary) {
    std::map<std::array<SimplexId, 3>, int> localTriangleMap;

    // loop through the internal cell list
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 3> triangleIds;

      // loop through each triangle of the cell
      for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
        triangleIds[0] = cellArray_->getCellVertex(cid, j);
        // the triangle does not belong to the current node
        if(triangleIds[0] > vertexIntervals_[clusterId]) {
          break;
        }
        for(SimplexId k = j + 1; k < maxCellDim_; k++) {
          for(SimplexId l = k + 1; l < maxCellDim_ + 1; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);
            localTriangleMap[triangleIds]++;
          }
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      // only check the cell[1,2,3] since cell[0] belongs to the previous cluster?
      std::array<SimplexId, 3> triangleIds;
        triangleIds[0] = cellArray_->getCellVertex(cid, 1);
      if(triangleIds[0] > vertexIntervals_[clusterId - 1]
            && triangleIds[0] <= vertexIntervals_[clusterId]) {
        triangleIds[1] = cellArray_->getCellVertex(cid, 2);
        triangleIds[2] = cellArray_->getCellVertex(cid, 3);
        localTriangleMap[triangleIds]++;
      }
    }

    // compute boundary triangles
    triangleCount = localTriangleMap.size();
    std::vector<bool> localBoundTriangles(triangleCount, false);
    int k = 0;
    for(auto &m : localTriangleMap) {
      clusterTriangles_[clusterId].push_back(std::move(m.first));
      if(m.second == 1) {
        localBoundTriangles[k] = true;
      }
      k++;
    }
    allClusters_[clusterId].boundaryTriangles_[0] = std::move(localBoundTriangles);
  }
  // otherwise, we just build the triangle list
  else {
    std::set<std::array<SimplexId, 3>> localTriangleSet;

    // loop through the internal cell list
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 3> triangleIds;

      // loop through each triangle of the cell
      for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
        triangleIds[0] = cellArray_->getCellVertex(cid, j);
        // the triangle does not belong to the current node
        if(triangleIds[0] > vertexIntervals_[clusterId]) {
          break;
        }
        for(SimplexId k = j + 1; k < maxCellDim_; k++) {
          for(SimplexId l = k + 1; l < maxCellDim_ + 1; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);

            if(localTriangleSet.find(triangleIds) == localTriangleSet.end()) {
              triangleCount++;
              localTriangleSet.insert(triangleIds);
            }
          }
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      // only check the cell[1,2,3] since cell[0] belongs to the previous cluster?
      std::array<SimplexId, 3> triangleIds;
        triangleIds[0] = cellArray_->getCellVertex(cid, 1);
      if(triangleIds[0] > vertexIntervals_[clusterId - 1]
            && triangleIds[0] <= vertexIntervals_[clusterId]) {
        triangleIds[1] = cellArray_->getCellVertex(cid, 2);
        triangleIds[2] = cellArray_->getCellVertex(cid, 3);
        if(localTriangleSet.find(triangleIds) == localTriangleSet.end()) {
          triangleCount++;
          localTriangleSet.insert(triangleIds);
        }
      }
    }
    clusterTriangles_[clusterId] = {localTriangleSet.begin(), localTriangleSet.end()};
  }
  return triangleCount;
}

int AcTopo::buildBoundaryTopCellList(const SimplexId &clusterId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  
  if(maxCellDim_ == 2) {
    SimplexId localEdgeNum = edgeIntervals_[clusterId] - edgeIntervals_[clusterId-1];
    SimplexId localVertexNum 
      = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];

    std::vector<bool> localBoundEdges(localEdgeNum, false);
    // [localVertexId][globalVertexId] = number of incident cells 
    std::vector<boost::unordered_map<SimplexId, int>> localEdgeStarMaps(localVertexNum);

    // loop through the internal cell list
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 2> edgeIds;
      for(SimplexId j = 0; j < maxCellDim_; j++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        // the edge does not belong to the current node
        if(edgeIds[0] > vertexIntervals_[clusterId]) {
          break;
        }
        SimplexId localVertexId = edgeIds[0] - vertexIntervals_[clusterId - 1] - 1;
        for(SimplexId k = j + 1; k < maxCellDim_+1; k++) {
          edgeIds[1] = cellArray_->getCellVertex(cid, k);
          localEdgeStarMaps[localVertexId][edgeIds[1]]++;
        }
      }
    }
    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 2> edgeIds;
      for(SimplexId j = 0; j < maxCellDim_; j++) {
        for(SimplexId k = j + 1; k < maxCellDim_ + 1; k++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          edgeIds[1] = cellArray_->getCellVertex(cid, k);
          if(edgeIds[0] > vertexIntervals_[clusterId - 1]
              && edgeIds[0] <= vertexIntervals_[clusterId]) {
            SimplexId localVertexId = edgeIds[0] - vertexIntervals_[clusterId - 1] - 1;
            localEdgeStarMaps[localVertexId][edgeIds[1]]++;
          }
        }
      }
    }

    for(SimplexId i = 0; i < localVertexNum; i++) {
      for(auto &ele : localEdgeStarMaps[i]) {
        if(ele.second == 1) {
          std::array<SimplexId, 2> edgeIds{i + vertexIntervals_[clusterId-1] + 1, ele.first};
          SimplexId idx = findEdgeIdx(clusterId, edgeIds);
          localBoundEdges[idx] = true;
        }
      }
    }
    cluster.boundaryEdges_[0] = std::move(localBoundEdges);
  }
  else if(maxCellDim_ == 3) {
    // get the boundary triangles first
    SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];
    SimplexId localVertexNum = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];

    std::vector<bool> localBoundTriangles(localTriangleNum, false);
    std::vector<boost::unordered_map<SimplexId, boost::unordered_map<SimplexId, int>>> localTriangleStarMaps(localVertexNum);


    // loop through the internal cell list
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 3> triangleIds;
      for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
        triangleIds[0] = cellArray_->getCellVertex(cid, j);
        // the triangle does not belong to the current node
        if(triangleIds[0] > vertexIntervals_[clusterId]) {
          break;
        }
        SimplexId localVertexId = triangleIds[0] - vertexIntervals_[clusterId-1] - 1;
        for(SimplexId k = j + 1; k < maxCellDim_ ; k++) {
          for(SimplexId l = k + 1; l < maxCellDim_ + 1; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);
            localTriangleStarMaps[localVertexId][triangleIds[1]][triangleIds[2]]++;
          }
        }
      }
    }
    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 3> triangleIds;
      triangleIds[0] = cellArray_->getCellVertex(cid, 1);
      if(triangleIds[0] > vertexIntervals_[clusterId - 1]
          && triangleIds[0] <= vertexIntervals_[clusterId]) {
        SimplexId localVertexId = triangleIds[0] - vertexIntervals_[clusterId-1] - 1;
        triangleIds[1] = cellArray_->getCellVertex(cid, 2);
        triangleIds[2] = cellArray_->getCellVertex(cid, 3);
        localTriangleStarMaps[localVertexId][triangleIds[1]][triangleIds[2]]++;
      }
    }

    for(SimplexId i = 0; i < localVertexNum; i++) {
      for(auto &m : localTriangleStarMaps[i]) {
        for(auto &ele : m.second) {
          if(ele.second == 1) {
            std::array<SimplexId, 3> triangleIds{i + vertexIntervals_[clusterId-1] + 1, m.first, ele.first};
            SimplexId idx = findTriangleIdx(clusterId, triangleIds);
            localBoundTriangles[idx] = true;
          }
        }
      }
    }
    cluster.boundaryTriangles_[0] = std::move(localBoundTriangles);

  } else {
    return -1;
  }

  return 0;
}

int AcTopo::getClusterCellEdges(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = 4;
  std::vector<std::array<SimplexId, 6>> localTetraEdges_(cellIntervals_[clusterId] - cellIntervals_[clusterId - 1]);
  ImplicitCluster &cluster = allClusters_.at(clusterId);

  for(SimplexId i = cellIntervals_[clusterId - 1] + 1;
      i <= cellIntervals_[clusterId]; i++) {
    int cnt = 0;
    // get the internal edge id from the map
    for(SimplexId k = 1; k < verticesPerCell; k++) {
      std::array<SimplexId, 2> edgePair
        = {(SimplexId)cellArray_->getCellVertex(i, 0),
            (SimplexId)cellArray_->getCellVertex(i, k)};
      SimplexId idx = findEdgeIdx(clusterId, edgePair);
      localTetraEdges_[i - cellIntervals_[clusterId - 1] - 1][cnt++] = idx + 1 + edgeIntervals_[clusterId - 1];
    }
    for(SimplexId j = 1; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        std::array<SimplexId, 2> edgePair
          = {(SimplexId)cellArray_->getCellVertex(i, j),
              (SimplexId)cellArray_->getCellVertex(i, k)};
        if(edgePair[0] <= vertexIntervals_[clusterId]) {
          SimplexId idx = findEdgeIdx(clusterId, edgePair);
          localTetraEdges_[i - cellIntervals_[clusterId - 1] - 1][cnt++]
            = idx + 1 + edgeIntervals_[clusterId - 1];
        } else {
          SimplexId nodeNum = vertexIndices_[edgePair[0]];
          SimplexId idx = findEdgeIdx(nodeNum, edgePair);
          localTetraEdges_[i - cellIntervals_[clusterId - 1] - 1][cnt++]
            = idx + 1 + edgeIntervals_[nodeNum - 1];
        }
      }
    }
  }

  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.tetraEdges_[threadId].empty()) {
      cluster.tetraEdges_[threadId] = std::move(localTetraEdges_);
    }
  }

  return 0;
}

int AcTopo::getClusterCellNeighbors(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  std::vector<std::vector<SimplexId>> localCellNeighbors(
    cellIntervals_[clusterId] - cellIntervals_[clusterId - 1]);
  ImplicitCluster &cluster = allClusters_.at(clusterId);

  // get the vertex stars into vector of vector
  SimplexId verticesPerCell = maxCellDim_ + 1;
  std::vector<std::vector<SimplexId>> localVertexStars;
  getClusterVertexStarsVector(clusterId, localVertexStars);

  boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
    nodeMaps;
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      if(cellArray_->getCellVertex(cid, j) > vertexIntervals_[clusterId]) {
        SimplexId nodeId
          = vertexIndices_[cellArray_->getCellVertex(cid, j)];
        if(nodeMaps.find(nodeId) == nodeMaps.end()) {
          std::vector<std::vector<SimplexId>> tmpVec;
          getClusterVertexStarsVector(nodeId, tmpVec);
          nodeMaps[nodeId] = std::move(tmpVec);
        }
      }
    }
  }

  if(maxCellDim_ == 2) {
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      for(SimplexId j = 0; j < 3; j++) {

        SimplexId v0 = cellArray_->getCellVertex(cid, j);
        SimplexId v1 = cellArray_->getCellVertex(cid, (j + 1) % 3);
        SimplexId localV0 = v0 - vertexIntervals_[clusterId - 1] - 1;
        SimplexId localV1 = v1 - vertexIntervals_[clusterId - 1] - 1;

        std::vector<SimplexId> stars0, stars1;
        if(v0 <= vertexIntervals_[clusterId]) {
          stars0 = localVertexStars[localV0];
        } else {
          localV0 = v0 - vertexIntervals_[vertexIndices_[v0] - 1] - 1;
          stars0 = nodeMaps[vertexIndices_[v0]][localV0];
        }
        if(v1 <= vertexIntervals_[clusterId]) {
          stars1 = localVertexStars[localV1];
        } else {
          localV1 = v1 - vertexIntervals_[vertexIndices_[v1] - 1] - 1;
          stars1 = nodeMaps[vertexIndices_[v1]][localV1];
        }

        // perform an intersection of the 2 sorted star lists
        SimplexId pos0 = 0, pos1 = 0;
        SimplexId intersection = -1;

        while((pos0 < (SimplexId)stars0.size())
              && (pos1 < (SimplexId)stars1.size())) {
          SimplexId biggest = stars0[pos0];
          if(stars1[pos1] > biggest) {
            biggest = stars1[pos1];
          }

          for(SimplexId l = pos0; l < (SimplexId)stars0.size(); l++) {
            if(stars0[l] < biggest) {
              pos0++;
            } else {
              break;
            }
          }
          for(SimplexId l = pos1; l < (SimplexId)stars1.size(); l++) {
            if(stars1[l] < biggest) {
              pos1++;
            } else {
              break;
            }
          }

          if(pos0 >= (SimplexId)stars0.size()
              || pos1 >= (SimplexId)stars1.size())
            break;

          if(stars0[pos0] == stars1[pos1]) {
            if(stars0[pos0] != cid) {
              intersection = stars0[pos0];
              break;
            }
            pos0++;
            pos1++;
          }
        }

        if(intersection != -1) {
          localCellNeighbors[cid - cellIntervals_[clusterId - 1] - 1]
            .push_back(intersection);
        }
      }
    }
  }

  else if(maxCellDim_ == 3) {
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      // go triangle by triangle
      for(SimplexId j = 0; j < 4; j++) {

        SimplexId v0 = cellArray_->getCellVertex(cid, j % 4);
        SimplexId v1 = cellArray_->getCellVertex(cid, (j + 1) % 4);
        SimplexId v2 = cellArray_->getCellVertex(cid, (j + 2) % 4);

        SimplexId localV0 = v0 - vertexIntervals_[clusterId - 1] - 1;
        SimplexId localV1 = v1 - vertexIntervals_[clusterId - 1] - 1;
        SimplexId localV2 = v2 - vertexIntervals_[clusterId - 1] - 1;

        std::vector<SimplexId> stars0, stars1, stars2;
        if(v0 <= vertexIntervals_[clusterId]) {
          stars0 = localVertexStars[localV0];
        } else {
          localV0 = v0 - vertexIntervals_[vertexIndices_[v0] - 1] - 1;
          stars0 = nodeMaps[vertexIndices_[v0]][localV0];
        }
        if(v1 <= vertexIntervals_[clusterId]) {
          stars1 = localVertexStars[localV1];
        } else {
          localV1 = v1 - vertexIntervals_[vertexIndices_[v1] - 1] - 1;
          stars1 = nodeMaps[vertexIndices_[v1]][localV1];
        }
        if(v2 <= vertexIntervals_[clusterId]) {
          stars2 = localVertexStars[localV2];
        } else {
          localV2 = v2 - vertexIntervals_[vertexIndices_[v2] - 1] - 1;
          stars2 = nodeMaps[vertexIndices_[v2]][localV2];
        }

        // perform an intersection of the 3 (sorted) star lists
        SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
        SimplexId intersection = -1;

        while((pos0 < (SimplexId)stars0.size())
              && (pos1 < (SimplexId)stars1.size())
              && (pos2 < (SimplexId)stars2.size())) {

          SimplexId biggest = stars0[pos0];
          if(stars1[pos1] > biggest) {
            biggest = stars1[pos1];
          }
          if(stars2[pos2] > biggest) {
            biggest = stars2[pos2];
          }

          for(SimplexId l = pos0; l < (SimplexId)stars0.size(); l++) {
            if(stars0[l] < biggest) {
              pos0++;
            } else {
              break;
            }
          }
          for(SimplexId l = pos1; l < (SimplexId)stars1.size(); l++) {
            if(stars1[l] < biggest) {
              pos1++;
            } else {
              break;
            }
          }
          for(SimplexId l = pos2; l < (SimplexId)stars2.size(); l++) {
            if(stars2[l] < biggest) {
              pos2++;
            } else {
              break;
            }
          }

          if(pos0 >= (SimplexId)stars0.size()
              || pos1 >= (SimplexId)stars1.size()
              || pos2 >= (SimplexId)stars2.size())
            break;

          if((stars0[pos0] == stars1[pos1])
              && (stars0[pos0] == stars2[pos2])) {
            if(stars0[pos0] != cid) {
              intersection = stars0[pos0];
              break;
            }
            pos0++;
            pos1++;
            pos2++;
          }
        }

        if(intersection != -1) {
          localCellNeighbors[cid - cellIntervals_[clusterId - 1] - 1]
            .push_back(intersection);
        }
      }
    }
  }

  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.cellNeighbors_[threadId].empty()) {
      cluster.cellNeighbors_[threadId] = FlatJaggedArray(localCellNeighbors);
    }
  }

  return 0;
}

int AcTopo::getClusterCellTriangles(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = 4;
  std::vector<std::array<SimplexId, 4>> localTetraTriangles_(cellIntervals_[clusterId] - cellIntervals_[clusterId - 1]);
  ImplicitCluster &cluster = allClusters_.at(clusterId);

  for(SimplexId i = cellIntervals_[clusterId - 1] + 1;
      i <= cellIntervals_[clusterId]; i++) {
    std::array<SimplexId, 3> triangleVec;
    // get the internal triangle from the map
    triangleVec[0] = cellArray_->getCellVertex(i, 0);
    for(SimplexId k = 1; k < verticesPerCell - 1; k++) {
      triangleVec[1] = cellArray_->getCellVertex(i, k);
      for(SimplexId l = k + 1; l < verticesPerCell; l++) {
        triangleVec[2] = cellArray_->getCellVertex(i, l);
        SimplexId idx = findTriangleIdx(clusterId, triangleVec);
        localTetraTriangles_[i - cellIntervals_[clusterId - 1] - 1][k + l - 3]
          = idx + 1 + triangleIntervals_[clusterId - 1];
      }
    }
    // group the external triangles by node id
    triangleVec[0] = cellArray_->getCellVertex(i, 1);
    triangleVec[1] = cellArray_->getCellVertex(i, 2);
    triangleVec[2] = cellArray_->getCellVertex(i, 3);
    if(triangleVec[0] <= vertexIntervals_[clusterId]) {
      SimplexId idx = findTriangleIdx(clusterId, triangleVec);
      localTetraTriangles_[i - cellIntervals_[clusterId - 1] - 1].back()
        = idx + 1 + triangleIntervals_[clusterId - 1];
    } else {
      SimplexId nodeNum = vertexIndices_[triangleVec[0]];
      SimplexId idx = findTriangleIdx(nodeNum, triangleVec);
      localTetraTriangles_[i - cellIntervals_[clusterId - 1] - 1].back()
        = idx + 1 + triangleIntervals_[nodeNum - 1];
    }
  }

  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.tetraTriangles_[threadId].empty()) {
      cluster.tetraTriangles_[threadId] = std::move(localTetraTriangles_);
    }
  }

  return 0;
}

int AcTopo::getClusterEdgeLinks(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  SimplexId localEdgeNum = edgeIntervals_[clusterId] - edgeIntervals_[clusterId-1];
  std::vector<SimplexId> offsets(localEdgeNum + 1, 0),
    linksCount(localEdgeNum, 0);

  if(maxCellDim_ == 2) {
    if(cluster.edgeStars_[threadId].empty()) {
      getClusterEdgeStars(clusterId, threadId);
    }
    // set the offsets std::vector
    for(SimplexId eid = 0; eid < localEdgeNum; eid++) {
      SimplexId globalEid = eid + edgeIntervals_[clusterId-1] + 1;
      for(SimplexId j = 0; j < cluster.edgeStars_[threadId].size(eid); j++) {
        SimplexId cellId = cluster.edgeStars_[threadId].get(eid, j);
        for(int k = 0; k < 3; k++) {
          SimplexId vertexId = cellArray_->getCellVertex(cellId, k);
          if((vertexId != edgeList_[globalEid][0]) 
              && (vertexId != edgeList_[globalEid][1])) {
            offsets[eid + 1]++;
            break;
          }
        }
      }
    }

    // compute partial sum of number of links per edge
    for(SimplexId i = 1; i <= localEdgeNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for edge link data
    std::vector<SimplexId> edgeLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId eid = 0; eid < localEdgeNum; eid++) {
      SimplexId globalEid = eid + edgeIntervals_[clusterId-1] + 1;
      for(SimplexId j = 0; j < cluster.edgeStars_[threadId].size(eid); j++) {
        SimplexId cellId = cluster.edgeStars_[threadId].get(eid, j);
        for(int k = 0; k < 3; k++) {
          SimplexId vertexId = cellArray_->getCellVertex(cellId, k);
          if((vertexId != edgeList_[globalEid][0]) 
              && (vertexId != edgeList_[globalEid][1])) {
            SimplexId localVertexId
              = vertexId - vertexIntervals_[clusterId - 1] - 1;
            edgeLinkData[offsets[localVertexId] + linksCount[localVertexId]]
              = vertexId;
            break;
          }
        }
      }
    }

    // fill FlatJaggedArray struct
    if(cluster.edgeLinks_[threadId].empty()) {
      cluster.edgeLinks_[threadId] = FlatJaggedArray(std::move(edgeLinkData), std::move(offsets));
    }

  } else if(maxCellDim_ == 3) {
    if(cluster.tetraEdges_[threadId].empty()) {
      getClusterCellEdges(clusterId, threadId);
    }

    // set the offsets vector
    SimplexId localCellNum = cellIntervals_[clusterId] - cellIntervals_[clusterId - 1];
    for(SimplexId cid = 0; cid < localCellNum; cid++) {
      SimplexId cellId = cid + cellIntervals_[clusterId - 1] + 1;
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cellId, 0),
            (SimplexId)cellArray_->getCellVertex(cellId, 1),
            (SimplexId)cellArray_->getCellVertex(cellId, 2),
            (SimplexId)cellArray_->getCellVertex(cellId, 3)};
      std::array<SimplexId, 2> edgePair;
      edgePair[0] = vertexIds[0];
      for(SimplexId j = 1; j < 4; j++) {
        edgePair[1] = vertexIds[j];
        SimplexId idx = findEdgeIdx(clusterId, edgePair);
        offsets[idx+1]++;
      }
      if(vertexIds[1] <= vertexIntervals_[clusterId]) {
        edgePair[0] = vertexIds[1];
        for(int j = 2; j < 4; j++) {
          edgePair[1] = vertexIds[j];
          SimplexId idx = findEdgeIdx(clusterId, edgePair);
          offsets[idx+1]++;
        }
        if(vertexIds[2] <= vertexIntervals_[clusterId]) {
          edgePair = {vertexIds[2], vertexIds[3]};
          SimplexId idx = findEdgeIdx(clusterId, edgePair);
          offsets[idx+1]++;
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 2> edgePair;
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
            (SimplexId)cellArray_->getCellVertex(cid, 1),
            (SimplexId)cellArray_->getCellVertex(cid, 2),
            (SimplexId)cellArray_->getCellVertex(cid, 3)};

      // loop through each edge of the cell
      for(SimplexId j = 0; j < 3; j++) {
        for(SimplexId k = j + 1; k < 4; k++) {
          edgePair[0] = vertexIds[j];
          edgePair[1] = vertexIds[k];

          // the edge is in the current node
          if(edgePair[0] > vertexIntervals_[clusterId - 1]
              && edgePair[0] <= vertexIntervals_[clusterId]) {
            SimplexId idx = findEdgeIdx(clusterId, edgePair);
            offsets[idx+1]++;
          }
        }
      }
    }

    // compute partial sum of number of neighbors per vertex
    for(SimplexId i = 1; i <= localEdgeNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for edge link data
    std::vector<SimplexId> edgeLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId cid = 0; cid < localCellNum; cid++) {
      SimplexId cellId = cid + cellIntervals_[clusterId - 1] + 1;
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cellId, 0),
            (SimplexId)cellArray_->getCellVertex(cellId, 1),
            (SimplexId)cellArray_->getCellVertex(cellId, 2),
            (SimplexId)cellArray_->getCellVertex(cellId, 3)};
      std::array<SimplexId, 2> edgePair;
      edgePair[0] = vertexIds[0];
      for(SimplexId j = 1; j < 4; j++) {
        edgePair[1] = vertexIds[j];
        SimplexId localEdgeId = findEdgeIdx(clusterId, edgePair);
        edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
          = cluster.tetraEdges_[threadId].at(cid)[6 - j];
        linksCount[localEdgeId]++;
      }
      if(vertexIds[1] <= vertexIntervals_[clusterId]) {
        edgePair[0] = vertexIds[1];
        for(int j = 2; j < 4; j++) {
          edgePair[1] = vertexIds[j];
          SimplexId localEdgeId = findEdgeIdx(clusterId, edgePair);
          edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
            = cluster.tetraEdges_[threadId].at(cid)[4 - j];
          linksCount[localEdgeId]++;
        }
        if(vertexIds[2] <= vertexIntervals_[clusterId]) {
          edgePair = {vertexIds[2], vertexIds[3]};
          SimplexId localEdgeId = findEdgeIdx(clusterId, edgePair);
          edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
            = cluster.tetraEdges_[threadId].at(cid)[0];
          linksCount[localEdgeId]++;
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
            (SimplexId)cellArray_->getCellVertex(cid, 1),
            (SimplexId)cellArray_->getCellVertex(cid, 2),
            (SimplexId)cellArray_->getCellVertex(cid, 3)};

      std::array<SimplexId, 2> edgeIds;
      // loop through each edge of the cell
      for(SimplexId j = 0; j < 3; j++) {
        for(SimplexId k = j + 1; k < 4; k++) {
          edgeIds[0] = vertexIds[j];
          edgeIds[1] = vertexIds[k];

          // the edge is in the current node
          if(edgeIds[0] > vertexIntervals_[clusterId - 1]
              && edgeIds[0] <= vertexIntervals_[clusterId]) {
            std::array<SimplexId, 2> otherEdge = {-1, -1};
            for(int i = 0; i < 4; i++) {
              if(vertexIds[i] != edgeIds[0] && vertexIds[i] != edgeIds[1]) {
                if(otherEdge[0] == -1) {
                  otherEdge[0] = vertexIds[i];
                } else if(otherEdge[1] == -1) {
                  otherEdge[1] = vertexIds[i];
                } else {
                  printErr("[AcTopo] More than two other vertices are "
                            "found in the edge!\n");
                }
              }
            }
            SimplexId nodeId = vertexIndices_[otherEdge[0]];
            SimplexId localEdgeId = findEdgeIdx(clusterId, edgeIds);
            edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
              = findEdgeIdx(nodeId, otherEdge);
            linksCount[localEdgeId]++;
          }
        }
      }
    }

    // fill FlatJaggedArray struct
    {
      std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
      if(cluster.edgeLinks_[threadId].empty()) {
        cluster.edgeLinks_[threadId]
          = FlatJaggedArray(std::move(edgeLinkData), std::move(offsets));
      }
    }
  }

  return 0;
}

int AcTopo::getClusterEdgeStars(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  SimplexId verticesPerCell = maxCellDim_ + 1;
  SimplexId localEdgeNum = edgeIntervals_[clusterId] - edgeIntervals_[clusterId-1];
  std::vector<std::vector<SimplexId>> localEdgeStarVec(localEdgeNum);

  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    std::array<SimplexId, 2> edgeIds;
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      edgeIds[0] = cellArray_->getCellVertex(cid, j);
      // the edge does not belong to the current node
      if(edgeIds[0] > vertexIntervals_[clusterId]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[1] = cellArray_->getCellVertex(cid, k);
        SimplexId idx = findEdgeIdx(clusterId, edgeIds);
        localEdgeStarVec[idx].push_back(cid);
      }
    }
  }
  // loop through the external cell list
  for(const SimplexId &cid : externalCells_[clusterId]) {
    std::array<SimplexId, 2> edgeIds;
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        edgeIds[1] = cellArray_->getCellVertex(cid, k);
        if(edgeIds[0] > vertexIntervals_[clusterId - 1]
            && edgeIds[0] <= vertexIntervals_[clusterId]) {
          SimplexId idx = findEdgeIdx(clusterId, edgeIds);
          localEdgeStarVec[idx].push_back(cid);
        }
      }
    }
  }

  std::vector<SimplexId> edgeStarData;
  std::vector<SimplexId> offsets;
  offsets.resize(localEdgeStarVec.size() + 1);
  for(size_t i = 0; i < localEdgeStarVec.size(); ++i) {
    offsets[i + 1] = offsets[i] + localEdgeStarVec[i].size();
  }
  edgeStarData.resize(offsets.back());
  for(size_t i = 0; i < localEdgeStarVec.size(); ++i) {
    for(size_t j = 0; j < localEdgeStarVec[i].size(); ++j) {
      edgeStarData[offsets[i] + j] = localEdgeStarVec[i][j];
    }
  }

  // fill FlatJaggedArray struct
  
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.edgeStars_[threadId].empty()) {
      cluster.edgeStars_[threadId] = FlatJaggedArray{std::move(edgeStarData), std::move(offsets)};
    }
  }

  return 0;
}

int AcTopo::getClusterEdgeTriangles(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  SimplexId localEdgeNum = edgeIntervals_[clusterId] - edgeIntervals_[clusterId-1];
  SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];
  std::vector<std::vector<SimplexId>> localEdgeTriangleVec(localEdgeNum);
  boost::unordered_set<std::array<SimplexId, 3>> exTriangleSet;

  // internal triangles
  for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
    SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
    std::array<SimplexId, 2> edge1 = {triangleList_[globalFid][0], triangleList_[globalFid][1]};
    std::array<SimplexId, 2> edge2 = {triangleList_[globalFid][0], triangleList_[globalFid][2]};
    SimplexId localEdgeId = findEdgeIdx(clusterId, edge1), triangleId = fid + 1 + triangleIntervals_[clusterId - 1];
    localEdgeTriangleVec[localEdgeId].push_back(triangleId);
    localEdgeId = findEdgeIdx(clusterId, edge2);
    localEdgeTriangleVec[localEdgeId].push_back(triangleId);

    if(triangleList_[globalFid][1] <= vertexIntervals_[clusterId]) {
      edge1 = {triangleList_[globalFid][1], triangleList_[globalFid][2]};
      localEdgeId = findEdgeIdx(clusterId, edge1);
      localEdgeTriangleVec[localEdgeId].push_back(triangleId);
    }
  }

  // loop through the external cell list
  for(const SimplexId &cid : externalCells_[clusterId]) {
    std::array<SimplexId, 3> triangleIds;
    // suppose cell (1,6,11,16) and 6 is in the current cluster 
    // triangle like (1,6,11) needs to be added for edge (6,11)
    for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] <= vertexIntervals_[clusterId - 1]) {
        for(SimplexId k = j + 1; k < maxCellDim_; k++) {
          for(SimplexId l = k + 1; l < maxCellDim_+1; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);
            if (exTriangleSet.count(triangleIds)) continue;
            if(triangleIds[1] > vertexIntervals_[clusterId-1]
                && triangleIds[1] <= vertexIntervals_[clusterId]) {
              std::array<SimplexId, 2> edgePair = {triangleIds[1], triangleIds[2]};
              SimplexId localEdgeId = findEdgeIdx(clusterId, edgePair);
              SimplexId nodeNum = vertexIndices_[triangleIds[0]];
              SimplexId localTriangleId = findTriangleIdx(nodeNum, triangleIds);
              localEdgeTriangleVec[localEdgeId].push_back(localTriangleId + 1 + triangleIntervals_[nodeNum-1]);
              exTriangleSet.insert(triangleIds);
            }
          }
        }
      }
    }
  }

  std::vector<SimplexId> edgeTraingleData;
  std::vector<SimplexId> offsets;
  offsets.resize(localEdgeTriangleVec.size() + 1);
  for(size_t i = 0; i < localEdgeTriangleVec.size(); ++i) {
    offsets[i + 1] = offsets[i] + localEdgeTriangleVec[i].size();
  }
  edgeTraingleData.resize(offsets.back());
  for(size_t i = 0; i < localEdgeTriangleVec.size(); ++i) {
    for(size_t j = 0; j < localEdgeTriangleVec[i].size(); ++j) {
      edgeTraingleData[offsets[i] + j] = localEdgeTriangleVec[i][j];
    }
  }

  // fill FlatJaggedArray struct
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.edgeTriangles_[threadId].empty()) {
      cluster.edgeTriangles_[threadId] = FlatJaggedArray(std::move(edgeTraingleData), std::move(offsets));
    }
  }

  return 0;
}

int AcTopo::getClusterTriangleEdges(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  std::vector<std::array<SimplexId, 3>> localTriangleEdges(triangleIntervals_[clusterId] - triangleIntervals_[clusterId - 1]);
  ImplicitCluster &cluster = allClusters_.at(clusterId);
  SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];

  for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
    // since the first vertex of the triangle is in the node ...
    SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
    std::array<SimplexId, 2> edgePair = {triangleList_[globalFid][0], triangleList_[globalFid][1]};
    localTriangleEdges[fid][0] = findEdgeIdx(clusterId, edgePair) + 1 + edgeIntervals_[clusterId - 1];
    edgePair[1] = triangleList_[globalFid][2];
    localTriangleEdges[fid][1] = findEdgeIdx(clusterId, edgePair) + 1 + edgeIntervals_[clusterId - 1];
    edgePair[0] = triangleList_[globalFid][1];
    if(edgePair[0] > vertexIntervals_[clusterId - 1]
        && edgePair[0] <= vertexIntervals_[clusterId]) {
      localTriangleEdges[fid][2] = findEdgeIdx(clusterId, edgePair) + 1 + edgeIntervals_[clusterId - 1];
    } else {
      SimplexId nodeNum = vertexIndices_[edgePair[0]];
      localTriangleEdges[fid][2] = findEdgeIdx(nodeNum, edgePair) + 1 + edgeIntervals_[nodeNum - 1];
    }
  }


  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.triangleEdges_[threadId].empty()) {
      cluster.triangleEdges_[threadId] = std::move(localTriangleEdges);
    }
  }

  return 0;
}

int AcTopo::getClusterTriangleLinks(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];
  std::vector<SimplexId> offsets(localTriangleNum + 1, 0),
    linksCount(localTriangleNum, 0);

  if(cluster.triangleStars_[threadId].empty()) {
    getClusterTriangleStars(clusterId, threadId);
  }

  // set the offsets vector
  for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
    SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
    for(SimplexId i = 0; i < cluster.triangleStars_[threadId].size(fid); i++) {
      SimplexId cellId = cluster.triangleStars_[threadId].get(fid, i);
      for(int j = 0; j < 4; j++) {
        SimplexId vertexId = cellArray_->getCellVertex(cellId, j);
        if((vertexId != triangleList_[globalFid][0]) 
            && (vertexId != triangleList_[globalFid][1])
            && (vertexId != triangleList_[globalFid][2])) {
          offsets[fid+1]++;
          break;
        }
      }
    }
  }

  // compute partial sum of number of links per triangle
  for(SimplexId i = 1; i <= localTriangleNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for triangle link data
  std::vector<SimplexId> triangleLinkData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
    SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
    for(SimplexId i = 0; i < cluster.triangleStars_[threadId].size(fid); i++) {
      SimplexId cellId = cluster.triangleStars_[threadId].get(fid, i);
      for(int j = 0; j < 4; j++) {
        SimplexId vertexId = cellArray_->getCellVertex(cellId, j);
        if((vertexId != triangleList_[globalFid][0]) 
            && (vertexId != triangleList_[globalFid][1])
            && (vertexId != triangleList_[globalFid][2])) {
          triangleLinkData[offsets[fid] + linksCount[fid]]
            = vertexId;
          linksCount[fid]++;
          break;
        }
      }
    }
  }

  // fill FlatJaggedArray struct
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.triangleLinks_[threadId].empty()) {
      cluster.triangleLinks_[threadId] = FlatJaggedArray(
        std::move(triangleLinkData), std::move(offsets));
    }
  }

  return 0;
}

int AcTopo::getClusterTriangleStars(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];
  FlatJaggedArray localTriangleStars;
  // use vector to avoid duplicate for-loops
  std::vector<std::vector<SimplexId>> localTriangleStarVec(localTriangleNum);

  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    std::array<SimplexId, 3> triangleIds;
    for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      // the triangle does not belong to the current node
      if(triangleIds[0] > vertexIntervals_[clusterId]) {
        break;
      }
      for(SimplexId k = j + 1; k < maxCellDim_; k++) {
        for(SimplexId l = k + 1; l < maxCellDim_ + 1; l++) {
          triangleIds[1] = cellArray_->getCellVertex(cid, k);
          triangleIds[2] = cellArray_->getCellVertex(cid, l);
          SimplexId localTriangleId = findTriangleIdx(clusterId, triangleIds);
          localTriangleStarVec[localTriangleId].push_back(cid);
        }
      }
    }
  }
  // loop through the external cell list
  for(const SimplexId &cid : externalCells_[clusterId]) {
    std::array<SimplexId, 3> triangleIds;
    triangleIds[0] = cellArray_->getCellVertex(cid, 1);
    if(triangleIds[0] > vertexIntervals_[clusterId - 1]
        && triangleIds[0] <= vertexIntervals_[clusterId]) {
      triangleIds[1] = cellArray_->getCellVertex(cid, 2);
      triangleIds[2] = cellArray_->getCellVertex(cid, 3);
      SimplexId localTriangleId = findTriangleIdx(clusterId, triangleIds);
      localTriangleStarVec[localTriangleId].push_back(cid);
    }
  }

  localTriangleStars.fillFrom(localTriangleStarVec);

  // fill FlatJaggedArray struct
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.triangleStars_[threadId].empty()) {
      cluster.triangleStars_[threadId] = std::move(localTriangleStars);
    }
  }

  return 0;
}

int AcTopo::getClusterVertexEdges(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId localVertexNum
    = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
  SimplexId localEdgeNum 
    = edgeIntervals_[clusterId] - edgeIntervals_[clusterId - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    edgesCount(localVertexNum, 0);
  ImplicitCluster &cluster = allClusters_.at(clusterId);
  boost::unordered_set<std::array<SimplexId, 2>> exEdgeSet;

  // set the offsets vector
  for(SimplexId eid = 0; eid < localEdgeNum; eid++) {
    SimplexId globalEid = eid + edgeIntervals_[clusterId-1] + 1;
    offsets[edgeList_[globalEid][0] - vertexIntervals_[clusterId - 1]]++;
    if(edgeList_[globalEid][1] <= vertexIntervals_[clusterId]) {
      offsets[edgeList_[globalEid][1] - vertexIntervals_[clusterId - 1]]++;
    }
  }
  // loop through the external cell list
  for(const SimplexId &cid : externalCells_[clusterId]) {
    std::array<SimplexId, 2> edgeIds;

    // loop through each edge of the cell
    for(SimplexId j = 0; j < maxCellDim_; j++) {
      edgeIds[0] = cellArray_->getCellVertex(cid, j);
      if(edgeIds[0] <= vertexIntervals_[clusterId-1]) {
        for(SimplexId k = j + 1; k <= maxCellDim_; k++) {
          edgeIds[1] = cellArray_->getCellVertex(cid, k);
          if(exEdgeSet.count(edgeIds)) continue;
          // check if the edge is an external edge
          if(edgeIds[1] > vertexIntervals_[clusterId - 1]
              && edgeIds[1] <= vertexIntervals_[clusterId]) {
            exEdgeSet.insert(edgeIds);
            offsets[edgeIds[1] - vertexIntervals_[clusterId-1]]++;
          }
        }
      }
    }
  }

  // compute partial sum of number of edges per vertex
  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for vertex edge data
  std::vector<SimplexId> vertexEdgeData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(SimplexId eid = 0; eid < localEdgeNum; eid++) {
    SimplexId globalEid = eid + edgeIntervals_[clusterId-1] + 1;
    SimplexId localVertexId = edgeList_[globalEid][0] - vertexIntervals_[clusterId - 1] - 1;
    vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
      = eid + 1 + edgeIntervals_[clusterId - 1];
    edgesCount[localVertexId]++;
    if(edgeList_[globalEid][1] <= vertexIntervals_[clusterId]) {
      localVertexId = edgeList_[globalEid][1] - vertexIntervals_[clusterId - 1] - 1;
      vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
        = eid + 1 + edgeIntervals_[clusterId - 1];
      edgesCount[localVertexId]++;
    }
  }

  for (auto &edgeIds : exEdgeSet) {
    SimplexId nodeNum = vertexIndices_[edgeIds[0]];
    SimplexId localVertexId = edgeIds[1] - vertexIntervals_[clusterId - 1] - 1;
    vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
      = findEdgeIdx(nodeNum, edgeIds) + 1 + edgeIntervals_[nodeNum-1];
    edgesCount[localVertexId]++;
  }

  // fill FlatJaggedArray struct
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.vertexEdges_[threadId].empty()) {
      cluster.vertexEdges_[threadId] = FlatJaggedArray(std::move(vertexEdgeData), std::move(offsets));
    }
  }

  return 0;
}

int AcTopo::getClusterVertexLinks(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId localVertexNum
    = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    linksCount(localVertexNum, 0);
  ImplicitCluster &cluster = allClusters_.at(clusterId);
  // triangle mesh
  if(maxCellDim_ == 2) {
    // set the offsets vector
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      offsets[cellArray_->getCellVertex(cid, 0)
              - vertexIntervals_[clusterId - 1]]++;
      for(SimplexId j = 1; j < 3; j++) {
        SimplexId vid = cellArray_->getCellVertex(cid, j);
        if(vid <= vertexIntervals_[clusterId]) {
          offsets[vid - vertexIntervals_[clusterId - 1]]++;
        }
      }
    }
    for(const SimplexId &cid : externalCells_[clusterId]) {
      for(SimplexId j = 1; j < 3; j++) {
        SimplexId vid = cellArray_->getCellVertex(cid, j);
        if(vid > vertexIntervals_[clusterId - 1]
            && vid <= vertexIntervals_[clusterId]) {
          offsets[vid - vertexIntervals_[clusterId - 1]]++;
        }
      }
    }

    // compute partial sum of number of links per vertex
    for(SimplexId i = 1; i <= localVertexNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for vertex link data
    std::vector<SimplexId> vertexLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 3> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
            (SimplexId)cellArray_->getCellVertex(cid, 1),
            (SimplexId)cellArray_->getCellVertex(cid, 2)};
      // the first vertex of the cell must be in the cluster
      std::array<SimplexId, 2> edgePair = {vertexIds[1], vertexIds[2]};
      SimplexId nodeId = vertexIndices_[vertexIds[1]];
      SimplexId localVertexId
        = vertexIds[0] - vertexIntervals_[clusterId - 1] - 1;
      vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
        = findEdgeIdx(nodeId, edgePair) + 1 + edgeIntervals_[nodeId - 1];
      linksCount[localVertexId]++;

      if(vertexIds[1] <= vertexIntervals_[clusterId]) {
        edgePair[0] = vertexIds[0];
        localVertexId = vertexIds[1] - vertexIntervals_[clusterId - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findEdgeIdx(clusterId, edgePair) + 1 + edgeIntervals_[clusterId - 1];
        linksCount[localVertexId]++;
        if(vertexIds[2] <= vertexIntervals_[clusterId]) {
          localVertexId
            = vertexIds[2] - vertexIntervals_[clusterId - 1] - 1;
          edgePair[1] = vertexIds[1];
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = findEdgeIdx(clusterId, edgePair) + 1 + edgeIntervals_[clusterId - 1];
          linksCount[localVertexId]++;
        }
      }
    }
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 3> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
            (SimplexId)cellArray_->getCellVertex(cid, 1),
            (SimplexId)cellArray_->getCellVertex(cid, 2)};
      std::array<SimplexId, 2> edgePair = {vertexIds[0], vertexIds[2]};
      SimplexId nodeId = vertexIndices_[edgePair[0]];
      SimplexId localVertexId
        = vertexIds[1] - vertexIntervals_[clusterId - 1] - 1;
      if(vertexIds[1] > vertexIntervals_[clusterId - 1]
          && vertexIds[1] <= vertexIntervals_[clusterId]) {
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findEdgeIdx(nodeId, edgePair) + 1 + edgeIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
      if(vertexIds[2] > vertexIntervals_[clusterId - 1]
          && vertexIds[2] <= vertexIntervals_[clusterId]) {
        edgePair[1] = vertexIds[1];
        localVertexId = vertexIds[2] - -vertexIntervals_[clusterId - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findEdgeIdx(nodeId, edgePair) + 1 + edgeIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
    }

    // fill FlatJaggedArray struct
    {
      std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
      if(cluster.vertexLinks_[threadId].empty()) {
        cluster.vertexLinks_[threadId] = FlatJaggedArray(
          std::move(vertexLinkData), std::move(offsets));
      }
    }

    // tetrahedral mesh
  } else if(maxCellDim_ == 3) {

    // set the offsets vector
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      offsets[cellArray_->getCellVertex(cid, 0)
              - vertexIntervals_[clusterId - 1]]++;
      for(SimplexId j = 1; j < 4; j++) {
        SimplexId vid = cellArray_->getCellVertex(cid, j);
        if(vid <= vertexIntervals_[clusterId]) {
          offsets[vid - vertexIntervals_[clusterId - 1]]++;
        }
      }
    }
    for(const SimplexId &cid : externalCells_[clusterId]) {
      for(SimplexId j = 1; j < 4; j++) {
        SimplexId vid = cellArray_->getCellVertex(cid, j);
        if(vid > vertexIntervals_[clusterId - 1]
            && vid <= vertexIntervals_[clusterId]) {
          offsets[vid - vertexIntervals_[clusterId - 1]]++;
        }
      }
    }

    // compute partial sum of number of links per vertex
    for(SimplexId i = 1; i <= localVertexNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for vertex link data
    std::vector<SimplexId> vertexLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
        cid <= cellIntervals_[clusterId]; cid++) {
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
            (SimplexId)cellArray_->getCellVertex(cid, 1),
            (SimplexId)cellArray_->getCellVertex(cid, 2),
            (SimplexId)cellArray_->getCellVertex(cid, 3)};

      // v1: (v2, v3, v4)
      std::array<SimplexId, 3> triangleVec
        = {vertexIds[1], vertexIds[2], vertexIds[3]};
      SimplexId nodeId = vertexIndices_[vertexIds[1]];
      SimplexId localVertexId
        = vertexIds[0] - vertexIntervals_[clusterId - 1] - 1;
      vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
        = findTriangleIdx(nodeId, triangleVec) + 1 + triangleIntervals_[nodeId - 1];
      linksCount[localVertexId]++;
      // v2: (v1, v3, v4)
      if(vertexIds[1] <= vertexIntervals_[clusterId]) {
        triangleVec[0] = vertexIds[0];
        localVertexId = vertexIds[1] - vertexIntervals_[clusterId - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findTriangleIdx(clusterId, triangleVec) + 1 + triangleIntervals_[clusterId - 1];
        linksCount[localVertexId]++;
        // v3: (v1, v2, v4)
        if(vertexIds[2] <= vertexIntervals_[clusterId]) {
          triangleVec[1] = vertexIds[1];
          localVertexId = vertexIds[2] - vertexIntervals_[clusterId - 1] - 1;
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = findTriangleIdx(clusterId, triangleVec) + 1 
              + triangleIntervals_[clusterId - 1];
          linksCount[localVertexId]++;
        }
        // v4: (v1, v2, v3)
        if(vertexIds[3] <= vertexIntervals_[clusterId]) {
          triangleVec[2] = vertexIds[2];
          localVertexId = vertexIds[3] - vertexIntervals_[clusterId - 1] - 1;
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = findTriangleIdx(clusterId, triangleVec) + 1
              + triangleIntervals_[clusterId - 1];
          linksCount[localVertexId]++;
        }
      }
    }

    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
            (SimplexId)cellArray_->getCellVertex(cid, 1),
            (SimplexId)cellArray_->getCellVertex(cid, 2),
            (SimplexId)cellArray_->getCellVertex(cid, 3)};
      // start from v2
      std::array<SimplexId, 3> triangleVec
        = {vertexIds[0], vertexIds[2], vertexIds[3]};
      SimplexId nodeId = vertexIndices_[vertexIds[0]];
      SimplexId localVertexId
        = vertexIds[1] - vertexIntervals_[clusterId - 1] - 1;
      if(vertexIds[1] > vertexIntervals_[clusterId - 1]
          && vertexIds[1] <= vertexIntervals_[clusterId]) {
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findTriangleIdx(nodeId, triangleVec) + 1
            + triangleIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
      if(vertexIds[2] > vertexIntervals_[clusterId - 1]
          && vertexIds[2] <= vertexIntervals_[clusterId]) {
        triangleVec[1] = vertexIds[1];
        localVertexId = vertexIds[2] - vertexIntervals_[clusterId - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findTriangleIdx(nodeId, triangleVec) + 1
            + triangleIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
      if(vertexIds[3] > vertexIntervals_[clusterId - 1]
          && vertexIds[3] <= vertexIntervals_[clusterId]) {
        triangleVec[1] = vertexIds[1];
        triangleVec[2] = vertexIds[2];
        localVertexId = vertexIds[3] - vertexIntervals_[clusterId - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = findTriangleIdx(nodeId, triangleVec) + 1
            + triangleIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
    }

    // fill FlatJaggedArray struct
    {
      std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
      if(cluster.vertexLinks_[threadId].empty()) {
        cluster.vertexLinks_[threadId] = FlatJaggedArray(
          std::move(vertexLinkData), std::move(offsets));
      }
    }
  }

  return 0;
}

int AcTopo::getClusterVertexNeighbors(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = maxCellDim_ + 1;
  SimplexId localVertexNum
    = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
  std::vector<SimplexId> vertexNeighborData, offsets(localVertexNum + 1, 0);
  ImplicitCluster &cluster = allClusters_.at(clusterId);

  SimplexId v1, v2;
  std::vector<boost::unordered_set<SimplexId>> vertexNeighborSet(
    localVertexNum);

  // loop through the internal cells
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      v1 = cellArray_->getCellVertex(cid, j);
      if(v1 <= vertexIntervals_[clusterId]) {
        for(SimplexId k = j + 1; k < verticesPerCell; k++) {
          v2 = cellArray_->getCellVertex(cid, k);
          vertexNeighborSet[v1 - vertexIntervals_[clusterId - 1] - 1]
            .insert(v2);
          if(v2 <= vertexIntervals_[clusterId]) {
            vertexNeighborSet[v2 - vertexIntervals_[clusterId - 1] - 1]
              .insert(v1);
          }
        }
      }
    }
  }

  // loop through external cells
  for(const SimplexId &cid : externalCells_[clusterId]) {
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        v1 = cellArray_->getCellVertex(cid, j);
        v2 = cellArray_->getCellVertex(cid, k);
        if(v1 > vertexIntervals_[clusterId - 1]
            && v1 <= vertexIntervals_[clusterId])
          vertexNeighborSet[v1 - vertexIntervals_[clusterId - 1] - 1]
            .insert(v2);
        if(v2 > vertexIntervals_[clusterId - 1]
            && v2 <= vertexIntervals_[clusterId])
          vertexNeighborSet[v2 - vertexIntervals_[clusterId - 1] - 1]
            .insert(v1);
      }
    }
  }

  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] = offsets[i - 1] + vertexNeighborSet[i - 1].size();
    vertexNeighborData.insert(vertexNeighborData.end(),
                              vertexNeighborSet[i-1].begin(),
                              vertexNeighborSet[i-1].end());
    // std::move(vertexNeighborSet[i-1].begin(), vertexNeighborSet[i-1].end(), vertexNeighborData.end());
  }

  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.vertexNeighbors_[threadId].empty()) {
      cluster.vertexNeighbors_[threadId] = FlatJaggedArray(
        std::move(vertexNeighborData), std::move(offsets));
    }
  }

  return 0;
}

int AcTopo::getClusterVertexStars(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = maxCellDim_ + 1;
  SimplexId localVertexNum
    = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    starsCount(localVertexNum, 0);
  ImplicitCluster &cluster = allClusters_.at(clusterId);

  // set the offsets vector
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    SimplexId vid = cellArray_->getCellVertex(cid, 0);
    offsets[vid - vertexIntervals_[clusterId - 1]]++;
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      vid = cellArray_->getCellVertex(cid, j);
      if(vid > vertexIntervals_[clusterId - 1]
          && vid <= vertexIntervals_[clusterId]) {
        offsets[vid - vertexIntervals_[clusterId - 1]]++;
      }
    }
  }
  for(const SimplexId &cid : externalCells_[clusterId]) {
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      SimplexId vid = cellArray_->getCellVertex(cid, j);
      if(vid > vertexIntervals_[clusterId - 1]
          && vid <= vertexIntervals_[clusterId]) {
        offsets[vid - vertexIntervals_[clusterId - 1]]++;
      }
    }
  }

  // compute partial sum of number of stars per vertex
  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for vertex star data
  std::vector<SimplexId> vertexStarData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    SimplexId localVertexId = cellArray_->getCellVertex(cid, 0)
                              - vertexIntervals_[clusterId - 1] - 1;
    vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
      = cid;
    starsCount[localVertexId]++;
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      // see if it is in the current node
      if(vertexId > vertexIntervals_[clusterId - 1]
          && vertexId <= vertexIntervals_[clusterId]) {
        localVertexId = vertexId - vertexIntervals_[clusterId - 1] - 1;
        vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
          = cid;
        starsCount[localVertexId]++;
      }
    }
  }
  for(const SimplexId &cid : externalCells_[clusterId]) {
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      // see if it is in the current node
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      if(vertexId > vertexIntervals_[clusterId - 1]
          && vertexId <= vertexIntervals_[clusterId]) {
        SimplexId localVertexId
          = vertexId - vertexIntervals_[clusterId - 1] - 1;
        vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
          = cid;
        starsCount[localVertexId]++;
      }
    }
  }

  // fill FlatJaggedArray struct
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.vertexStars_[threadId].empty()) {
      cluster.vertexStars_[threadId]
        = FlatJaggedArray(std::move(vertexStarData), std::move(offsets));
    }
  }

  return 0;
}

int AcTopo::getClusterVertexStarsVector(const SimplexId &clusterId, 
                std::vector<std::vector<SimplexId>> &vertexStarVec) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  // get the vertex stars into a vector of vector
  SimplexId verticesPerCell = maxCellDim_ + 1;
  SimplexId localVertexNum
    = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
  vertexStarVec = std::vector<std::vector<SimplexId>>(localVertexNum);

  // internal cells
  for(SimplexId cid = cellIntervals_[clusterId - 1] + 1;
      cid <= cellIntervals_[clusterId]; cid++) {
    SimplexId localVertexId = cellArray_->getCellVertex(cid, 0)
                              - vertexIntervals_[clusterId - 1] - 1;
    
    vertexStarVec[localVertexId].push_back(cid);
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      // see if it is in the current node
      if(vertexId > vertexIntervals_[clusterId - 1]
          && vertexId <= vertexIntervals_[clusterId]) {
        localVertexId = vertexId - vertexIntervals_[clusterId - 1] - 1;
        vertexStarVec[localVertexId].push_back(cid);
      }
    }
  }
  for(const SimplexId &cid : externalCells_[clusterId]) {
    for(SimplexId j = 0; j < verticesPerCell; j++) {
      // see if it is in the current node
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      if(vertexId > vertexIntervals_[clusterId - 1]
          && vertexId <= vertexIntervals_[clusterId]) {
        SimplexId localVertexId
          = vertexId - vertexIntervals_[clusterId - 1] - 1;
        vertexStarVec[localVertexId].push_back(cid);
      }
    }
  }

  for(size_t i = 0; i < vertexStarVec.size(); i++) {
    std::sort(vertexStarVec[i].begin(), vertexStarVec[i].end());
  }

  return 0;
}

int AcTopo::getClusterVertexTriangles(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  SimplexId localVertexNum
    = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
  SimplexId localTriangleNum
    = triangleIntervals_[clusterId] - triangleIntervals_[clusterId - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    trianglesCount(localVertexNum, 0);
  boost::unordered_set<std::array<SimplexId, 3>> exTriangleSet{};
  ImplicitCluster &cluster = allClusters_.at(clusterId);

  // set the offsets vector
  for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
    SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
    for(SimplexId j = 0; j < 3; j++) {
      SimplexId vid = triangleList_[globalFid][j];
      if(vid > vertexIntervals_[clusterId - 1]
          && vid <= vertexIntervals_[clusterId])
        offsets[vid - vertexIntervals_[clusterId - 1]]++;
    }
  }
  // loop through the external cell list
  for(const SimplexId &cid : externalCells_[clusterId]) {
    std::array<SimplexId, 3> triangleIds;
    // suppose cell (1,6,11,16) and 6 is in the current cluster 
    for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] <= vertexIntervals_[clusterId - 1]) {
        for(SimplexId k = j + 1; k < maxCellDim_; k++) {
          for(SimplexId l = k + 1; l < maxCellDim_+1; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);
            if (exTriangleSet.count(triangleIds)) continue;
            bool added = false;
            if(triangleIds[1] > vertexIntervals_[clusterId-1]
                && triangleIds[1] <= vertexIntervals_[clusterId]) {
              exTriangleSet.insert(triangleIds);
              added = true;
              offsets[triangleIds[1] - vertexIntervals_[clusterId-1]]++;
            }
            if(triangleIds[2] > vertexIntervals_[clusterId-1]
                && triangleIds[2] <= vertexIntervals_[clusterId]) {
              if(!added) exTriangleSet.insert(triangleIds);
              offsets[triangleIds[2] - vertexIntervals_[clusterId-1]]++;
            }
          }
        }
      }
    }
  }

  // compute partial sum of number of triangles per vertex
  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for vertex triangle data
  std::vector<SimplexId> vertexTriangleData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
    SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
    for(SimplexId j = 0; j < 3; j++) {
      if(triangleList_[globalFid][j] > vertexIntervals_[clusterId - 1]
          && triangleList_[globalFid][j] <= vertexIntervals_[clusterId]) {
        SimplexId localVertexId
          = triangleList_[globalFid][j] - vertexIntervals_[clusterId - 1] - 1;
        vertexTriangleData[offsets[localVertexId] + trianglesCount[localVertexId]]
          = fid + 1 + triangleIntervals_[clusterId - 1];
        trianglesCount[localVertexId]++;
      }
    }
  }
  for(const auto &triangleIds : exTriangleSet) {
    SimplexId nodeNum = vertexIndices_[triangleIds[0]];
    for(SimplexId j = 1; j < 3; j++) {
      if(triangleIds[j] > vertexIntervals_[clusterId - 1]
          && triangleIds[j] <= vertexIntervals_[clusterId]) {
        SimplexId localVertexId = triangleIds[j] - vertexIntervals_[clusterId - 1] - 1;
        vertexTriangleData[offsets[localVertexId]
                            + trianglesCount[localVertexId]]
          = findTriangleIdx(nodeNum, triangleIds) + 1 + triangleIntervals_[nodeNum - 1];
        trianglesCount[localVertexId]++;
      }
    }
  }

  // fill FlatJaggedArray struct
  {
    std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
    if(cluster.vertexTriangles_[threadId].empty()) {
      cluster.vertexTriangles_[threadId] = FlatJaggedArray(
        std::move(vertexTriangleData), std::move(offsets));
    }
  }

  return 0;
}

int AcTopo::getClusterBoundaryVertices(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  if(cluster.boundaryVertices_[threadId].empty()) {
    SimplexId localVertexNum = vertexIntervals_[clusterId] - vertexIntervals_[clusterId - 1];
    SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];
    std::vector<bool> localBoundVertices(localVertexNum, false);
    if(maxCellDim_ == 2) {
      SimplexId localEdgeNum = edgeIntervals_[clusterId] - edgeIntervals_[clusterId-1];
      // internal edges
      for(SimplexId eid = 0; eid < localEdgeNum; eid++) {
        SimplexId globalEid = eid + edgeIntervals_[clusterId-1] + 1;
        if(cluster.boundaryEdges_[0][eid]) {
          localBoundVertices[edgeList_[globalEid][0] - vertexIntervals_[clusterId - 1] - 1]
            = true;
          if(edgeList_[globalEid][1] <= vertexIntervals_[clusterId]) {
            localBoundVertices
              [edgeList_[globalEid][1] - vertexIntervals_[clusterId - 1] - 1]
              = true;
          }
        }
      }
      // loop through the external cell list
      for(const SimplexId &cid : externalCells_[clusterId]) {
        std::array<SimplexId, 2> edgeIds;
        // loop through each edge of the cell
        for(SimplexId j = 0; j < maxCellDim_; j++) {
          edgeIds[0] = cellArray_->getCellVertex(cid, j);
          if(edgeIds[0] <= vertexIntervals_[clusterId-1]) {
            for(SimplexId k = j + 1; k < maxCellDim_+1; k++) {
              edgeIds[1] = cellArray_->getCellVertex(cid, k);
              // check if the edge is an external edge
              if(edgeIds[1] > vertexIntervals_[clusterId - 1]
                  && edgeIds[1] <= vertexIntervals_[clusterId]) {
                SimplexId nodeNum = vertexIndices_[edgeIds[0]];
                SimplexId idx = findEdgeIdx(nodeNum, edgeIds);
                if(allClusters_[nodeNum].boundaryEdges_[0][idx]) {
                  localBoundVertices[edgeIds[1] - vertexIntervals_[clusterId - 1] - 1] = true;
                }
              }
            }
          }
        }
      }
    } else if (maxCellDim_ == 3) {
      // internal triangles
      for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
        SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
        if(cluster.boundaryTriangles_[0][fid]) {
          for(int j = 0; j < 3; j++) {
            SimplexId vid = triangleList_[globalFid][j];
            if(vid <= vertexIntervals_[clusterId]) {
              localBoundVertices[vid - vertexIntervals_[clusterId - 1] - 1]
                = true;
            }
          }
        }
      }
      // loop through the external cell list
      for(const SimplexId &cid : externalCells_[clusterId]) {
        std::array<SimplexId, 3> triangleIds;
        // suppose cell (1,6,11,16) and 6 is in the current cluster 
        // triangle like (1,6,11) needs to be added for vertex 6
        for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
          triangleIds[0] = cellArray_->getCellVertex(cid, j);
          if(triangleIds[0] <= vertexIntervals_[clusterId - 1]) {
            for(SimplexId k = j + 1; k < maxCellDim_; k++) {
              for(SimplexId l = k + 1; l < maxCellDim_+1; l++) {
                triangleIds[1] = cellArray_->getCellVertex(cid, k);
                triangleIds[2] = cellArray_->getCellVertex(cid, l);
                SimplexId nodeNum = vertexIndices_[triangleIds[0]];
                if(triangleIds[1] > vertexIntervals_[clusterId-1]
                    && triangleIds[1] <= vertexIntervals_[clusterId]) {
                  SimplexId localTriangleId = findTriangleIdx(nodeNum, triangleIds);
                  if(allClusters_[nodeNum].boundaryTriangles_[0][localTriangleId])
                    localBoundVertices[triangleIds[1] - vertexIntervals_[clusterId-1] - 1] = true;
                }
                if(triangleIds[2] > vertexIntervals_[clusterId-1]
                    && triangleIds[2] <= vertexIntervals_[clusterId]) {
                    SimplexId localTriangleId = findTriangleIdx(nodeNum, triangleIds);
                    if(allClusters_[nodeNum].boundaryTriangles_[0][localTriangleId])
                      localBoundVertices[triangleIds[2] - vertexIntervals_[clusterId-1] - 1] = true;
                }
              }
            }
          }
        }
      }
    }

    {
      std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
      if(cluster.boundaryVertices_[threadId].empty())
        cluster.boundaryVertices_[threadId] = std::move(localBoundVertices);
    }
  }

  return 0;
}

/**
 * Get the boundary edges in a given node.
 */
int AcTopo::getClusterBoundaryEdges(const SimplexId &clusterId, const ThreadId &threadId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(clusterId <= 0 || clusterId > nodeNumber_)
    return -1;
#endif

  ImplicitCluster &cluster = allClusters_.at(clusterId);
  if(cluster.boundaryEdges_[threadId].empty()) {
    SimplexId localEdgeNum = edgeIntervals_[clusterId] - edgeIntervals_[clusterId-1];
    SimplexId localTriangleNum = triangleIntervals_[clusterId] - triangleIntervals_[clusterId-1];

    std::vector<bool> localBoundEdges(localEdgeNum, false);
    for(SimplexId fid = 0; fid < localTriangleNum; fid++) {
      SimplexId globalFid = fid + triangleIntervals_[clusterId-1] + 1;
      if(cluster.boundaryTriangles_[0][fid]) {
        std::array<SimplexId, 2> edgePair
          = {triangleList_[globalFid][0], 
              triangleList_[globalFid][1]};
        SimplexId localEdgeId = findEdgeIdx(clusterId, edgePair);
        localBoundEdges[localEdgeId] = true;
        edgePair[1] = triangleList_[globalFid][2];
        localEdgeId = findEdgeIdx(clusterId, edgePair);
        localBoundEdges[localEdgeId] = true;
        if(triangleList_[globalFid][1] <= vertexIntervals_[clusterId]) {
          edgePair[0] = triangleList_[globalFid][1];
          localEdgeId = findEdgeIdx(clusterId, edgePair);
          localBoundEdges[localEdgeId] = true;
        }
      }
    }
    // loop through the external cell list
    for(const SimplexId &cid : externalCells_[clusterId]) {
      std::array<SimplexId, 3> triangleIds;
      for(SimplexId j = 0; j < maxCellDim_ - 1; j++) {
        triangleIds[0] = cellArray_->getCellVertex(cid, j);
        if(triangleIds[0] <= vertexIntervals_[clusterId - 1]) {
          for(SimplexId k = j + 1; k < maxCellDim_; k++) {
            for(SimplexId l = k + 1; l < maxCellDim_+1; l++) {
              triangleIds[1] = cellArray_->getCellVertex(cid, k);
              triangleIds[2] = cellArray_->getCellVertex(cid, l);
              if(triangleIds[1] > vertexIntervals_[clusterId-1]
                  && triangleIds[1] <= vertexIntervals_[clusterId]) {
                SimplexId nodeNum = vertexIndices_[triangleIds[0]];
                SimplexId idx = findTriangleIdx(nodeNum, triangleIds);
                if(allClusters_[nodeNum].boundaryTriangles_[0][idx]) {
                  std::array<SimplexId, 2> edgePair = {triangleIds[1], triangleIds[2]};
                  localBoundEdges[findEdgeIdx(clusterId, edgePair)] = true;
                }
              }
            }
          }
        }
      }
    }

    {
      std::lock_guard<std::mutex> clck(clusterMutexes_[clusterId]);
      if(cluster.boundaryEdges_[threadId].empty())
        cluster.boundaryEdges_[threadId] = std::move(localBoundEdges);
    }
  }
  
  return 0;
}

void AcTopo::computeRelation(const SimplexId &clusterId,
                      const RelationType &relation, 
                      const ThreadId &threadId) {
  ImplicitCluster &cluster = allClusters_.at(clusterId);
  switch(relation) {

    /* internal lists */
    case RelationType::IEList:
      if(edgeList_.empty()) {
        edgeIntervals_[clusterId] = buildInternalEdgeList(clusterId);
      }
      break;

    case RelationType::IFList:
      if(triangleList_.empty()) {
        triangleIntervals_[clusterId] = buildInternalTriangleList(clusterId);
      }
      break;
    
    case RelationType::IBEList:
      if(edgeList_.empty()) {
        edgeIntervals_[clusterId] = buildInternalEdgeList(clusterId, true);
      }
      break;

    case RelationType::IBFList:
      if(triangleList_.empty()) {
        triangleIntervals_[clusterId] = buildInternalTriangleList(clusterId, true);
      }
      break;

    case RelationType::IBCList:
      if(maxCellDim_ == 3 && cluster.boundaryTriangles_[0].empty()) {
        buildBoundaryTopCellList(clusterId);
      }
      else if(maxCellDim_ == 2 && cluster.boundaryEdges_[0].empty()) {
        buildBoundaryTopCellList(clusterId);
      }
      break;

    /* vertex-related relations */
    case RelationType::VERelation:
      if(cluster.vertexEdges_[threadId].empty()) {
        getClusterVertexEdges(clusterId, threadId);
      }
      break;

    case RelationType::VFRelation:
      if(cluster.vertexTriangles_[threadId].empty()) {
        getClusterVertexTriangles(clusterId, threadId);
      }
      break;
    
    case RelationType::VLRelation:
      if(cluster.vertexLinks_[threadId].empty()) {
        getClusterVertexLinks(clusterId, threadId);
      }
      break;

    case RelationType::VTRelation:
      if(cluster.vertexStars_[threadId].empty()) {
        getClusterVertexStars(clusterId, threadId);
      }
      break;

    case RelationType::VVRelation:
      if(cluster.vertexNeighbors_[threadId].empty()) {
        getClusterVertexNeighbors(clusterId, threadId);
      }
      break;

    /* edge-related relations */
    case RelationType::ELRelation:
      if(cluster.edgeLinks_[threadId].empty()) {
        getClusterEdgeLinks(clusterId, threadId);
      }
      break;

    case RelationType::EFRelation:
      if(cluster.edgeTriangles_[threadId].empty()) {
        getClusterEdgeTriangles(clusterId, threadId);
      }
      break;

    case RelationType::ETRelation:
      if(cluster.edgeStars_[threadId].empty()) {
        getClusterEdgeStars(clusterId, threadId);
      }
      break;

    /* triangle-related relations */
    case RelationType::FERelation:
      if(cluster.triangleEdges_[threadId].empty()) {
        getClusterTriangleEdges(clusterId, threadId);
      }
      break;

    case RelationType::FLRelation:
      if(cluster.triangleLinks_[threadId].empty()) {
        getClusterTriangleLinks(clusterId, threadId);
      }
      break;

    case RelationType::FTRelation:
      if(cluster.triangleStars_[threadId].empty()) {
        getClusterTriangleStars(clusterId, threadId);
      }
      break;

    /* tetrahedron-related relations */
    // reserved for TV relation
    case RelationType::TVRelation:
      break;

    case RelationType::TERelation:
      if(cluster.tetraEdges_[threadId].empty()) {
        getClusterCellEdges(clusterId, threadId);
      }
      break;

    case RelationType::TFRelation:
      if(cluster.tetraTriangles_[threadId].empty()) {
        getClusterCellTriangles(clusterId, threadId);
      }
      break;

    case RelationType::TTRelation:
      if(cluster.cellNeighbors_[threadId].empty()) {
        getClusterCellNeighbors(clusterId, threadId);
      }
      break;

    /* boundary relations */
    case RelationType::BVRelation:
      if(cluster.boundaryVertices_[threadId].empty()) {
        getClusterBoundaryVertices(clusterId, threadId);
      }
      break;

    case RelationType::BERelation:
      if(cluster.boundaryEdges_[threadId].empty()) {
        getClusterBoundaryEdges(clusterId, threadId);
      }
      break;

    // case RelationType::BFRelation:
    //   if(cluster.boundaryTriangles_.empty()) {
    //     getBoundaryCells(clusterId, 2);
    //   }
    //   break;

    default:
      break;
  }
}

void AcTopo::leaderProcedure(const int &consumerId) {
  // create the worker threads
  const int startIdx = consumerId * numProducers_;
  SimplexId leaderCluster = 0;
  RelationType leaderFunc = RelationType::BlankRelation;

  if(numProducers_ > 1) {
    // Initialize the worker semaphores
    for(int i = 1; i < numProducers_; i++) {
      if(sem_init(&semaphores_[startIdx + i], 0, 0)) {
        this->printErr("Cannot initialize the worker semaphore vector for " + std::to_string(consumerId) + "-th leader producer!");
      }
    }

    // Start the worker threads 
    changed_[consumerId] = false;
    sharedClusterIds_[consumerId] = 0;
    for(int i = 1; i < numProducers_; i++) {
      producers_[startIdx + i] = std::thread(&AcTopo::workerProcedure, this, i, consumerId);
    }

    while(true) {
      std::unique_lock<std::mutex> llck(leaderMutexes_[consumerId]);
      while(!preconditions_[consumerId] && !waitings_[consumerId]) {
        leaderCondVars_[consumerId].wait(llck);
      }

      // preconditioning mode
      if(preconditions_[consumerId]) {
        if(leaderFunc != reqRelations_[consumerId]) {
          // set the preconditioning for workers
          leaderFunc = reqRelations_[consumerId];
          std::unique_lock<std::mutex> wlck(workerMutexes_[consumerId]);
          workerRelations_[consumerId] = leaderFunc;
          for(int i = 1; i < numProducers_; i++) 
            finished_[startIdx + i] = 0;
          changed_[consumerId] = 1;
          wlck.unlock();
          workerCondVars_[consumerId].notify_all();

          SimplexId startCluster = clustersPerThread_ * startIdx + 1;
          SimplexId endCluster = startCluster + clustersPerThread_;
          for(SimplexId nid = startCluster; nid < endCluster; nid++) {
            computeRelation(nid, leaderFunc, consumerId);
          }
          
          // wait for the worker producers
          for(int i = 1; i < numProducers_; i++) {
            sem_wait(&semaphores_[startIdx + i]);
          }

          sem_post(&semaphores_[startIdx]);
          changed_[consumerId] = 0;
          wlck.lock();
          preconditions_[consumerId] = 0;
          wlck.unlock();
        }
      }

      // computing mode
      else if(waitings_[consumerId]) {
        if(leaderCluster != reqClusters_[consumerId] || leaderFunc != reqRelations_[consumerId]) {
          changed_[consumerId] = 0;      // make sure the worker will stop for the leader
          leaderCluster = reqClusters_[consumerId], leaderFunc = reqRelations_[consumerId];
          if(leaderCluster == -1) {
            std::unique_lock<std::mutex> wlck(workerMutexes_[consumerId]);
            sharedClusterIds_[consumerId] = -1;
            if(workMode_ > 2) {
              workerClusterVecs_[consumerId].push_back(-1), workerClusterIds_[consumerId] = 0;
            }
            changed_[consumerId] = 1;
            wlck.unlock();
            workerCondVars_[consumerId].notify_all();
            return;
          }
          computeRelation(leaderCluster, leaderFunc, consumerId);
        
          // update the worker cluster and function
          {
            std::lock_guard<std::mutex> wlck(workerMutexes_[consumerId]);
            sharedClusterIds_[consumerId] = leaderCluster;
            if(workMode_ == 1) {
              workerRelations_[consumerId] = leaderFunc;
            }
            else if(workMode_ == 2) {
              workerRelationIds_[consumerId] = 0;
            }
            else if(workMode_ == 3) {
              workerClusterVecs_[consumerId] = connectivity_[sharedClusterIds_[consumerId]];
              workerClusterIds_[consumerId] = 0, workerRelationIds_[consumerId] = 0;
            }
            else if(workMode_ == 4) {
              workerClusterVecs_[consumerId] = connectivity_[sharedClusterIds_[consumerId]];
              workerClusterIds_[consumerId] = 0, workerRelations_[consumerId] = leaderFunc;
            }
            changed_[consumerId] = 1;
          }
          workerCondVars_[consumerId].notify_all();
          sem_post(&semaphores_[startIdx]);

          // update the buffer
          std::unique_lock<std::mutex> bufferLock(bufferMutexes_[consumerId]);
          if(bufferSets_[consumerId].find(leaderCluster) == bufferSets_[consumerId].end()) {
            // clean the buffer first
            if(sbuffers_[consumerId].size() >= bufferSize_) {
              for(size_t i = 0; i < bufferSize_/2; i++) {
              allClusters_[sbuffers_[consumerId].front()].clear(consumerId);
              bufferSets_[consumerId].erase(sbuffers_[consumerId].front());
              sbuffers_[consumerId].pop_front();
              }
            }
            // add the current cluster to the buffer
            sbuffers_[consumerId].push_back(leaderCluster);
            bufferSets_[consumerId].insert(leaderCluster);
          }
          bufferLock.unlock();
        }
      }
      llck.unlock();
    }
  }

  else {
    while(true) {
      std::unique_lock<std::mutex> llck(leaderMutexes_[consumerId]);
      while(!preconditions_[consumerId] && !waitings_[consumerId]) {
        leaderCondVars_[consumerId].wait(llck);
      }

      // preconditioning mode
      if(preconditions_[consumerId]) {
        if(leaderFunc != reqRelations_[consumerId]) {
          // set the preconditioning for workers
          leaderFunc = reqRelations_[consumerId];
          SimplexId startCluster = clustersPerThread_ * startIdx + 1;
          SimplexId endCluster = startCluster + clustersPerThread_;
          for(SimplexId nid = startCluster; nid < endCluster; nid++) {
            computeRelation(nid, leaderFunc, consumerId);
          }
          sem_post(&semaphores_[startIdx]);
          preconditions_[consumerId] = 0;
        }
      }

      // computing mode
      else if(waitings_[consumerId]) {
        if(leaderCluster != reqClusters_[consumerId] || leaderFunc != reqRelations_[consumerId]) {
          leaderCluster = reqClusters_[consumerId], leaderFunc = reqRelations_[consumerId];
          if(leaderCluster == -1) {
            return;
          }
          computeRelation(leaderCluster, leaderFunc, consumerId);
          sem_post(&semaphores_[startIdx]);

          // update the buffer 
          std::unique_lock<std::mutex> bufferLock(bufferMutexes_[consumerId]);
          if(bufferSets_[consumerId].find(leaderCluster) == bufferSets_[consumerId].end()) {
            // clean the buffer first
            if(sbuffers_[consumerId].size() >= bufferSize_) {
              for(size_t i = 0; i < bufferSize_/2; i++) {
              allClusters_[sbuffers_[consumerId].front()].clear(consumerId);
              bufferSets_[consumerId].erase(sbuffers_[consumerId].front());
              sbuffers_[consumerId].pop_front();
              }
            }
            // add the current cluster to the buffer
            sbuffers_[consumerId].push_back(leaderCluster);
            bufferSets_[consumerId].insert(leaderCluster);
          }
          bufferLock.unlock();
        }
      }
      llck.unlock();
    }
  }

}

void AcTopo::preconditionFunc(const int &workerId, const int &consumerId) {
  const int startIdx = consumerId * numProducers_;
  int isFinished = finished_[startIdx + workerId];
  if(isFinished == 0) {
    SimplexId start = clustersPerThread_* (startIdx + workerId) + 1;
    SimplexId end = start + clustersPerThread_;
    for(SimplexId nid = start; nid < end; nid++) {
      computeRelation(nid, workerRelations_[consumerId], consumerId);
    }
    sem_post(&semaphores_[startIdx + workerId]);
    // make sure all workers stop after finishing the preconditioning task
    finished_[startIdx + workerId] = 1;
  }
}

void AcTopo::workerProcedure(const int &workerId, const int &consumerId) {
  
  SimplexId clusterId = 0;
  RelationType functionId = RelationType::BlankRelation;

  if(workMode_ == 1) {
    while(true) {
      std::unique_lock<std::mutex> wlock(workerMutexes_[consumerId]);
      while(!changed_[consumerId]) {
        workerCondVars_[consumerId].wait(wlock);
      }
      if(preconditions_[consumerId] ) {
        wlock.unlock();
        preconditionFunc(workerId, consumerId);
      }
      else if(sharedClusterIds_[consumerId] < 0) {
        wlock.unlock();
        return;
      }
      else {
        if(sbuffers_[consumerId].size() >= bufferSize_) {
          wlock.unlock();
          continue;
        }
        else if(sharedClusterIds_[consumerId] <= nodeNumber_) {
          clusterId = sharedClusterIds_[consumerId]++, functionId = workerRelations_[consumerId];
          wlock.unlock();
          computeRelation(clusterId, functionId, consumerId);
          // update the buffer
          {
            std::lock_guard<std::mutex> blck(bufferMutexes_[consumerId]);
            if(!bufferSets_[consumerId].count(clusterId)) {
              sbuffers_[consumerId].push_back(clusterId);
              bufferSets_[consumerId].insert(clusterId);
            }
          }
        }
        else {
          wlock.unlock();
        }
      }
    }
  }
  else if(workMode_ == 2) {
    while(true) {
      std::unique_lock<std::mutex> wlock(workerMutexes_[consumerId]);
      while(!changed_[consumerId]) {
        workerCondVars_[consumerId].wait(wlock);
      }
      if(preconditions_[consumerId] ) {
        wlock.unlock();
        preconditionFunc(workerId, consumerId);
      }
      else if(sharedClusterIds_[consumerId] < 0) {
        wlock.unlock();
        return;
      }
      else {
        if(sbuffers_[consumerId].size() >= bufferSize_) {
          wlock.unlock();
          continue;
        }
        else if(sharedClusterIds_[consumerId] <= nodeNumber_) {
          clusterId = sharedClusterIds_[consumerId] , functionId = relationVec_[workerRelationIds_[consumerId]++];
          if(workerRelationIds_[consumerId] >= relationVec_.size()) {
            workerRelationIds_[consumerId] = 0, sharedClusterIds_[consumerId] ++;
          }
          wlock.unlock();
          computeRelation(clusterId, functionId, consumerId);
          // update the buffer
          {
            std::lock_guard<std::mutex> blck(bufferMutexes_[consumerId]);
            if(!bufferSets_[consumerId].count(clusterId)) {
              sbuffers_[consumerId].push_back(clusterId);
              bufferSets_[consumerId].insert(clusterId);
            }
          }
        }
        else {
          wlock.unlock();
        }
      }
    }
  }
  else if(workMode_ == 3) {
    while(true) {
      std::unique_lock<std::mutex> wlock(workerMutexes_[consumerId]);
      while(!changed_[consumerId]) {
        workerCondVars_[consumerId].wait(wlock);
      }
      if(preconditions_[consumerId] ) {
        wlock.unlock();
        preconditionFunc(workerId, consumerId);
      }
      else if(sharedClusterIds_[consumerId] < 0) {
        wlock.unlock();
        return;
      }
      else {
        if(sbuffers_[consumerId].size() >= bufferSize_) {
          wlock.unlock();
          continue;
        }
        else if(!workerClusterVecs_[consumerId].empty()) {
          if(workerClusterIds_[consumerId] < workerClusterVecs_[consumerId].size()) {
            clusterId = workerClusterVecs_[consumerId][workerClusterIds_[consumerId]], functionId = relationVec_[workerRelationIds_[consumerId]++];
            if(workerRelationIds_[consumerId] >= relationVec_.size()) {
              workerRelationIds_[consumerId] = 0, workerClusterIds_[consumerId]++;
            }
            wlock.unlock();
            computeRelation(clusterId, functionId, consumerId);
            // update the buffer
            {
              std::lock_guard<std::mutex> blck(bufferMutexes_[consumerId]);
              if(!bufferSets_[consumerId].count(clusterId)) {
                sbuffers_[consumerId].push_back(clusterId);
                bufferSets_[consumerId].insert(clusterId);
              }
            }
          }
        }
        else {
          wlock.unlock();
        }
      }
    }
  }
  else if(workMode_ == 4) {
    while(true) {
      std::unique_lock<std::mutex> wlock(workerMutexes_[consumerId]);
      while(!changed_[consumerId]) {
        workerCondVars_[consumerId].wait(wlock);
      }
      if(preconditions_[consumerId] ) {
        wlock.unlock();
        preconditionFunc(workerId, consumerId);
      }
      else if(sharedClusterIds_[consumerId] < 0) {
        wlock.unlock();
        return;
      }
      else {
        if(sbuffers_[consumerId].size() >= bufferSize_) {
          wlock.unlock();
          continue;
        }
        else if(!workerClusterVecs_[consumerId].empty()) {
          if(workerClusterIds_[consumerId] < workerClusterVecs_[consumerId].size()) {
            clusterId = workerClusterVecs_[consumerId][workerClusterIds_[consumerId]++], functionId = workerRelations_[consumerId];
            wlock.unlock();
            computeRelation(clusterId, functionId, consumerId);
            // update the buffer
            {
              std::lock_guard<std::mutex> blck(bufferMutexes_[consumerId]);
              if(!bufferSets_[consumerId].count(clusterId)) {
                sbuffers_[consumerId].push_back(clusterId);
                bufferSets_[consumerId].insert(clusterId);
              }
            }
          }
        }
        else {
          wlock.unlock();
        }
      }
    }
  }
}
