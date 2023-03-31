
#pragma once

// base code includes
#include <Debug.h>
#include <MemoryUsage.h>
#include <Triangulation.h>

namespace ttk {

  class TestTopoRelations : public virtual Debug {

  public:
    TestTopoRelations() {
      this->setDebugMsgPrefix("TestTopoRelations");
    }

    ~TestTopoRelations(){};

    template <class triangulationType = AbstractTriangulation>
    int execute(const triangulationType *triangulation) const;

    inline int preconditionTriangulation(Triangulation *triangulation) {
      if(triangulation) {
        Timer t;
        // build edges and triangles
        triangulation->preconditionEdges();
        triangulation->preconditionTriangles();
        // // boundary relations
        // triangulation->preconditionBoundaryVertices();
        // triangulation->preconditionBoundaryEdges();
        // triangulation->preconditionBoundaryTriangles();
        // vertex related relations
        triangulation->preconditionVertexEdges();
        triangulation->preconditionVertexStars();
        triangulation->preconditionVertexNeighbors();
        triangulation->preconditionVertexTriangles();
        // edge related relations
        triangulation->preconditionEdgeStars();
        triangulation->preconditionEdgeTriangles();
        // triangle related relations
        triangulation->preconditionTriangleEdges();
        triangulation->preconditionTriangleStars();
        // cell related relations
        triangulation->preconditionCellEdges();
        triangulation->preconditionCellNeighbors();
        triangulation->preconditionCellTriangles();
        // // links
        // triangulation->preconditionVertexLinks();
        // triangulation->preconditionEdgeLinks();
        // triangulation->preconditionTriangleLinks();

        std::cout << "[TestTopoRelations] Time usage for preconditioning: "
                  << t.getElapsedTime() << " s." << std::endl;
      }
      return 0;
    }
  };
} // namespace ttk

// template functions
template <class triangulationType>
int ttk::TestTopoRelations::execute(
  const triangulationType *triangulation) const {

  Timer t, tot;

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation)
    return -1;
#endif

  std::cout << "[TestTopoRelations] Dimension: "
            << triangulation->getDimensionality() << std::endl;

  SimplexId vertexNumber = triangulation->getNumberOfVertices();
  SimplexId edgeNumber = triangulation->getNumberOfEdges();
  SimplexId triangleNumber = triangulation->getNumberOfTriangles();
  SimplexId cellNumber = triangulation->getNumberOfCells();

  std::cout << "[TestTopoRelations] vertex num: " << vertexNumber
            << ", edge num: " << edgeNumber
            << ", triangle num: " << triangleNumber
            << ", cell num: " << cellNumber << std::endl;

//   // VL relation
//   if(triangulation->getDimensionality() == 2) {
//     for(SimplexId vid = 0; vid < vertexNumber; vid++) {
//       int linkNum = triangulation->getVertexLinkNumber(vid);
//       for(int i = 0; i < linkNum; i++) {
//         SimplexId linkId;
//         int res = triangulation->getVertexLink(vid, i, linkId);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get vertex link for vertex id "
//                     << vid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         for(int j = 0; j < 2; j++) {
//           SimplexId vertexId;
//           triangulation->getEdgeVertex(linkId, j, vertexId);
//           if(vertexId == vid) {
//             std::cout << "[TestTopoRelations] Find vertex " << vid << " in link "
//                       << linkId << std::endl;
//             return -4;
//           }
//         }
//       }
//     }
//   } else if(triangulation->getDimensionality() == 3) {
//     for(SimplexId vid = 0; vid < vertexNumber; vid++) {
//       int linkNum = triangulation->getVertexLinkNumber(vid);
//       for(int i = 0; i < linkNum; i++) {
//         SimplexId linkId;
//         int res = triangulation->getVertexLink(vid, i, linkId);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get vertex link for vertex id "
//                     << vid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         for(int j = 0; j < 3; j++) {
//           SimplexId vertexId;
//           triangulation->getTriangleVertex(linkId, j, vertexId);
//           if(vertexId == vid) {
//             std::cout << "[TestTopoRelations] Find vertex " << vid << " in link "
//                       << linkId << std::endl;
//             return -4;
//           }
//         }
//       }
//     }
//   }
//   printMsg("Finish testing VL relation!");

//   // EL relation
//   if(triangulation->getDimensionality() == 2) {
//     for(SimplexId eid = 0; eid < edgeNumber; eid++) {
//       int linkNum = triangulation->getEdgeLinkNumber(eid);
//       for(int i = 0; i < linkNum; i++) {
//         SimplexId linkId;
//         int res = triangulation->getEdgeLink(eid, i, linkId);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get edge link for edge id "
//                     << eid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         for(int j = 0; j < 2; j++) {
//           SimplexId vertexId;
//           triangulation->getEdgeVertex(eid, j, vertexId);
//           if(vertexId == linkId) {
//             std::cout << "[TestTopoRelations] Find vertex " << vertexId
//                       << " as link " << linkId << std::endl;
//             return -4;
//           }
//         }
//       }
//     }
//   } else if(triangulation->getDimensionality() == 3) {
//     for(SimplexId eid = 0; eid < edgeNumber; eid++) {
//       int linkNum = triangulation->getEdgeLinkNumber(eid);
//       for(int i = 0; i < linkNum; i++) {
//         SimplexId linkId;
//         int res = triangulation->getEdgeLink(eid, i, linkId);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get edge link for edge id "
//                     << eid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         if(eid == linkId) {
//           std::cout << "[TestTopoRelations] Find edge " << eid << " as link "
//                     << linkId << std::endl;
//           break;
//         }
//       }
//     }
//   }
//   printMsg("Finish testing EL relation!");

//   // FL relation
//   if(triangulation->getDimensionality() == 3) {
//     for(SimplexId tid = 0; tid < triangleNumber; tid++) {
//       int linkNum = triangulation->getTriangleLinkNumber(tid);
//       for(int i = 0; i < linkNum; i++) {
//         SimplexId linkId;
//         int res = triangulation->getTriangleLink(tid, i, linkId);
//         if(res) {
//           std::cout
//             << "[TestTopoRelations] Cannot get triangle link for triangle id "
//             << tid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         for(int j = 0; j < 3; j++) {
//           SimplexId vertexId;
//           triangulation->getTriangleVertex(tid, j, vertexId);
//           if(vertexId == linkId) {
//             std::cout << "[TestTopoRelations] Find vertex " << vertexId
//                       << " as link " << linkId << std::endl;
//             break;
//           }
//         }
//       }
//     }
//   }
//   printMsg("Finish testing FL relation!");

  // // Boundary relations
  // t.reStart();
  // int boundaryVertexNum = 0;
  // for(SimplexId vid = 0; vid < vertexNumber; vid++) {
  //   if(triangulation->isVertexOnBoundary(vid)) {
  //     boundaryVertexNum++;
  //   }
  // }
  // std::cout << "[TestTopoRelations] Boundary vertex number: " << boundaryVertexNum
  //           << std::endl;

  // int boundaryEdgeNum = 0;
  // for(SimplexId eid = 0; eid < edgeNumber; eid++) {
  //   if(triangulation->isEdgeOnBoundary(eid)) {
  //     boundaryEdgeNum++;
  //   }
  // }
  // std::cout << "[TestTopoRelations] Boundary edge number: " << boundaryEdgeNum
  //           << std::endl;

  // int boundaryTriangleNum = 0;
  // for(SimplexId tid = 0; tid < triangleNumber; tid++) {
  //   if(triangulation->isTriangleOnBoundary(tid)) {
  //     boundaryTriangleNum++;
  //   }
  // }
  // std::cout << "[TestTopoRelations] Boundary triangle number: "
  //           << boundaryTriangleNum << std::endl;
  // std::cout << "[TestTopoRelations] Time usage for getting boundaries: "
  //           << t.getElapsedTime() << " s." << std::endl;

//   // VE relation
//   for(SimplexId vid = 0; vid < vertexNumber; vid++) {
//     int edgeNum = triangulation->getVertexEdgeNumber(vid);
//     for(int i = 0; i < edgeNum; i++) {
//       SimplexId edgeId;
//       triangulation->getVertexEdge(vid, i, edgeId);
//       if(edgeId < 0) {
//         std::cout << "[TestTopoRelations] Cannot get vertex edge for vertex id "
//                   << vid << ", error code: " << edgeId << std::endl;
//         return -4;
//       }
//       bool hasFound = false;
//       for(int j = 0; j < 2; j++) {
//         SimplexId vertexId;
//         triangulation->getEdgeVertex(edgeId, j, vertexId);
//         if(vertexId == vid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find vertex " << vid
//                   << " in edge " << edgeId << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing VE relation!");

//   // VF relation
//   for(SimplexId vid = 0; vid < vertexNumber; vid++) {
//     int triangleNum = triangulation->getVertexTriangleNumber(vid);
//     for(int i = 0; i < triangleNum; i++) {
//       SimplexId triangleId;
//       int res = triangulation->getVertexTriangle(vid, i, triangleId);
//       if(res) {
//         std::cout
//           << "[TestTopoRelations] Cannot get vertex triangle for vertex id "
//           << vid << ", error code: " << res << std::endl;
//         return -4;
//       }
//       bool hasFound = false;
//       for(int j = 0; j < 3; j++) {
//         SimplexId vertexId;
//         triangulation->getTriangleVertex(triangleId, j, vertexId);
//         if(vertexId == vid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find vertex " << vid
//                   << " in triangle " << triangleId << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing VF relation!");

  // // VT relation
  // for(SimplexId vid = 0; vid < vertexNumber; vid++) {
  //   int starNum = triangulation->getVertexStarNumber(vid);
  //   for(int i = 0; i < starNum; i++) {
  //     SimplexId cellId;
  //     int res = triangulation->getVertexStar(vid, i, cellId);
  //     if(res) {
  //       std::cout << "[TestTopoRelations] Cannot get vertex edge for vertex id "
  //                 << vid << ", error code: " << res << std::endl;
  //       return -4;
  //     }
  //     bool hasFound = false;
  //     for(int j = 0; j < 4; j++) {
  //       SimplexId vertexId;
  //       triangulation->getCellVertex(cellId, j, vertexId);
  //       if(vertexId == vid) {
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound) {
  //       std::cout << "[TestTopoRelations] Cannot find vertex " << vid
  //                 << " in cell " << cellId << std::endl;
  //       break;
  //     }
  //   }
  // }
  // printMsg("Finish testing VT relation!");

//   // EV relation
//   for(SimplexId eid = 0; eid < edgeNumber; eid++) {
//     for(int i = 0; i < 2; i++) {
//       SimplexId vid;
//       triangulation->getEdgeVertex(eid, i, vid);
//       SimplexId edgeNum = triangulation->getVertexEdgeNumber(vid);
//       bool hasFound = false;
//       for(int j = 0; j < edgeNum; j++) {
//         SimplexId edgeid;
//         int res = triangulation->getVertexEdge(vid, j, edgeid);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get vertex edge for vertex id "
//                     << vid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         if(edgeid == eid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find edge " << eid
//                   << " in vertex " << vid << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing EV relation!");

//   // EF relation
//   for(SimplexId eid = 0; eid < edgeNumber; eid++) {
//     SimplexId triangleNum = triangulation->getEdgeTriangleNumber(eid);
//     for(int i = 0; i < triangleNum; i++) {
//       SimplexId triangleId;
//       triangulation->getEdgeTriangle(eid, i, triangleId);
//       if(triangleId < 0) {
//         std::cout << "[TestTopoRelations] Cannot get edge triangle for edge id "
//                   << eid << ", error code: " << triangleId << std::endl;
//         return -4;
//       }
//       bool findEdge = false;
//       SimplexId edgeId;
//       for(int j = 0; j < 3; j++) {
//         triangulation->getTriangleEdge(triangleId, j, edgeId);
//         if(edgeId == eid) {
//           findEdge = true;
//           break;
//         }
//       }
//       if(!findEdge) {
//         std::cout << "[TestTopoRelations] Cannot find edge " << eid
//                   << " in triangle " << triangleId << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing EF relation!");

//   // ET relation
//   for(SimplexId eid = 0; eid < edgeNumber; eid++) {
//     SimplexId starNum = triangulation->getEdgeStarNumber(eid);
//     for(int i = 0; i < starNum; i++) {
//       SimplexId cellId;
//       int res = triangulation->getEdgeStar(eid, i, cellId);
//       if(res) {
//         std::cout << "[TestTopoRelations] Cannot get edge star for edge id "
//                   << eid << ", error code: " << res << std::endl;
//         return -4;
//       }
//       bool findEdge = false;
//       SimplexId edgeId;
//       for(int j = 0; j < 6; j++) {
//         triangulation->getCellEdge(cellId, j, edgeId);
//         if(edgeId == eid) {
//           findEdge = true;
//           break;
//         }
//       }
//       if(!findEdge) {
//         std::cout << "[TestTopoRelations] Cannot find edge " << eid << " in cell "
//                   << cellId << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing ET relation!");

//   // FV relation
//   for(SimplexId tid = 0; tid < triangleNumber; tid++) {
//     for(int i = 0; i < 3; i++) {
//       SimplexId vid;
//       triangulation->getTriangleVertex(tid, i, vid);
//       SimplexId triangleNum = triangulation->getVertexTriangleNumber(vid);
//       bool hasFound = false;
//       for(int j = 0; j < triangleNum; j++) {
//         SimplexId triangleid;
//         int res = triangulation->getVertexTriangle(vid, j, triangleid);
//         if(res) {
//           std::cout
//             << "[TestTopoRelations] Cannot get vertex triangle for vertex id "
//             << vid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         if(triangleid == tid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find triangle " << tid
//                   << " in vertex " << vid << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing FV relation!");

//   // FE relation
//   for(SimplexId tid = 0; tid < triangleNumber; tid++) {
//     for(int i = 0; i < 3; i++) {
//       SimplexId eid;
//       triangulation->getTriangleEdge(tid, i, eid);
//       SimplexId triangleNum = triangulation->getEdgeTriangleNumber(eid);
//       bool hasFound = false;
//       for(int j = 0; j < triangleNum; j++) {
//         SimplexId triangleid;
//         int res = triangulation->getEdgeTriangle(eid, j, triangleid);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get edge triangle for edge id "
//                     << eid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         if(triangleid == tid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find triangle " << tid
//                   << " in edge " << eid << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing FE relation!");

//   // FT relation
//   for(SimplexId tid = 0; tid < triangleNumber; tid++) {
//     SimplexId starNum = triangulation->getTriangleStarNumber(tid);
//     for(int i = 0; i < starNum; i++) {
//       SimplexId cellId;
//       int res = triangulation->getTriangleStar(tid, i, cellId);
//       if(res) {
//         std::cout
//           << "[TestTopoRelations] Cannot get traingle star for traingle id "
//           << tid << ", error code: " << res << std::endl;
//         return -4;
//       }
//       bool findTriangle = false;
//       SimplexId triangleid;
//       for(int j = 0; j < 4; j++) {
//         triangulation->getCellTriangle(cellId, j, triangleid);
//         if(triangleid == tid) {
//           findTriangle = true;
//           break;
//         }
//       }
//       if(!findTriangle) {
//         std::cout << "[TestTopoRelations] Cannot find triangle " << tid
//                   << " in cell " << cellId << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing FT relation!");

//   // TV relation
//   for(SimplexId cid = 0; cid < cellNumber; cid++) {
//     for(int i = 0; i < 4; i++) {
//       SimplexId vid;
//       triangulation->getCellVertex(cid, i, vid);
//       SimplexId starNum = triangulation->getVertexStarNumber(vid);
//       bool hasFound = false;
//       for(int j = 0; j < starNum; j++) {
//         SimplexId starid;
//         int res = triangulation->getVertexStar(vid, j, starid);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get vertex star for vertex id "
//                     << vid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         if(starid == cid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find cell " << cid
//                   << " in vertex " << vid << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing TV relation!");

//   // TE relation
//   for(SimplexId cid = 0; cid < cellNumber; cid++) {
//     for(int i = 0; i < 6; i++) {
//       SimplexId eid;
//       triangulation->getCellEdge(cid, i, eid);
//       SimplexId starNum = triangulation->getEdgeStarNumber(eid);
//       bool hasFound = false;
//       for(int j = 0; j < starNum; j++) {
//         SimplexId starid;
//         int res = triangulation->getEdgeStar(eid, j, starid);
//         if(res) {
//           std::cout << "[TestTopoRelations] Cannot get edge star for edge id "
//                     << eid << ", error code: " << res << std::endl;
//           return -4;
//         }
//         if(starid == cid) {
//           hasFound = true;
//           break;
//         }
//       }
//       if(!hasFound) {
//         std::cout << "[TestTopoRelations] Cannot find cell " << cid << " in edge "
//                   << eid << std::endl;
//         break;
//       }
//     }
//   }
//   printMsg("Finish testing TE relation!");

//   // TF relation
//   if(triangulation->getDimensionality() == 3) {
//     for(SimplexId cid = 0; cid < cellNumber; cid++) {
//       for(int i = 0; i < 4; i++) {
//         SimplexId tid;
//         triangulation->getCellTriangle(cid, i, tid);
//         SimplexId starNum = triangulation->getTriangleStarNumber(tid);
//         bool hasFound = false;
//         for(int j = 0; j < starNum; j++) {
//           SimplexId starid;
//           int res = triangulation->getTriangleStar(tid, j, starid);
//           if(res) {
//             std::cout
//               << "[TestTopoRelations] Cannot get traingle star for traingle id "
//               << tid << ", error code: " << res << std::endl;
//             return -4;
//           }
//           if(starid == cid) {
//             hasFound = true;
//             break;
//           }
//         }
//         if(!hasFound) {
//           std::cout << "[TestTopoRelations] Cannot find cell " << cid
//                     << " in traingle " << tid << std::endl;
//           break;
//         }
//       }
//     }
//   }
//   printMsg("Finish testing TF relation!");
  
//   // cell neighbor relation
//   t.reStart();
// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(threadNumber_)
// // #endif
//   for(SimplexId cellId = 0; cellId < cellNumber; cellId++) {
//     SimplexId neighborId;
//     SimplexId neighborNum = triangulation->getCellNeighborNumber(cellId);
//     for(SimplexId j = 0; j < neighborNum; j++) {
//       int result1 = triangulation->getCellNeighbor(cellId, j, neighborId);
//       if(result1) {
//         std::cout << "[TestTopoRelations] cellId " << cellId
//                   << " Something wrong in getCellNeighbor()! Error code: "
//                   << result1 << std::endl;
//       }
//     }
//   }
//   std::cout << "[TestTopoRelations] Time usage for CN: " << t.getElapsedTime()
//             << " s." << std::endl;

//   // Vertex neighbor relation
//   t.reStart();
// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(threadNumber_)
// // #endif
//   for(SimplexId vid = 0; vid < vertexNumber; vid++) {
//     SimplexId neighborId;
//     SimplexId neighborNum = triangulation->getVertexNeighborNumber(vid);
//     for(int j = 0; j < neighborNum; j++) {
//       triangulation->getVertexNeighbor(vid, j, neighborId);
//     }
//   }
//   std::cout << "[TestTopoRelations] Time usage for VN: " << t.getElapsedTime()
//             << " s." << std::endl;
  
  /* test vertex tetrahedron relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++) {
    SimplexId cellNum = triangulation->getVertexStarNumber(vertexId);
    SimplexId cellId;
    for(SimplexId j = 0; j < cellNum; j++) {
      int result1 = triangulation->getVertexStar(vertexId, j, cellId);
      if(result1) {
        std::cout << "[TestTopoRelations] vertexId " << vertexId
                  << " Something wrong in getVertexStar()! Error code: "
                  << result1 << std::endl;
        return 1;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for VT: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test edge vertex relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++) {
    SimplexId vertexId;
    for(SimplexId j = 0; j < 2; j++) {
      int result1 = triangulation->getEdgeVertex(edgeId, j, vertexId);
      if(result1) {
        std::cout << "[TestTopoRelations] edgeId " << edgeId
                  << " Something wrong in getEdgeVertex()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for EV: " << t.getElapsedTime()
            << " s." << std::endl;
  
  /* test edge tetrahedron relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++) {
    SimplexId cellNum = triangulation->getEdgeStarNumber(edgeId);
    SimplexId cellId;
    for(SimplexId j = 0; j < cellNum; j++) {
      int result1 = triangulation->getEdgeStar(edgeId, j, cellId);
      if(result1) {
        std::cout << "[TestTopoRelations] edgeId " << edgeId
                  << " Something wrong in getEdgeStar()! Error code: "
                  << result1 << std::endl;
        return 1;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for ET: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test triangle vertex relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++) {
    SimplexId vertexId;
    for(SimplexId j = 0; j < 3; j++) {
      int result1 = triangulation->getTriangleVertex(triangleId, j, vertexId);
      if(result1) {
        std::cout << "[TestTopoRelations] triangleId " << triangleId
                  << " Something wrong in getTriangleVertex()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for FV: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test triangle tetrahedron relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++) {
    SimplexId cellNum = triangulation->getTriangleStarNumber(triangleId);
    SimplexId cellId;
    for(SimplexId j = 0; j < cellNum; j++) {
      int result1 = triangulation->getTriangleStar(triangleId, j, cellId);
      if(result1) {
        std::cout << "[TestTopoRelations] triangleId " << triangleId
                  << " Something wrong in getTriangleStar()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for FT: " << t.getElapsedTime()
            << " s." << std::endl;

  SimplexId verticesPerCell = triangulation->getCellVertexNumber(0);

  /* test tetrahedron vertex relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId cellId = 0; cellId < cellNumber; cellId++) {
    SimplexId vertexId;
    for(SimplexId j = 0; j < verticesPerCell; j++) {
      int result1 = triangulation->getCellVertex(cellId, j, vertexId);
      if(result1) {
        std::cout << "[TestTopoRelations] cellId " << cellId
                  << " Something wrong in getCellVertex()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for TV: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test vertex edge relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++) {
    SimplexId edgeNum = triangulation->getVertexEdgeNumber(vertexId);
    SimplexId edgeId;
    // if(std::is_same<triangulationType, ttk::AcTopo>::value) {
    //   AcTopo *topo = (AcTopo *)triangulation;
    //   topo->printBuffer();
    // }
    
    for(SimplexId j = 0; j < edgeNum; j++) {
      int result1 = triangulation->getVertexEdge(vertexId, j, edgeId);
      if(result1) {
        std::cerr << "[TestTopoRelations] vertexId " << vertexId
                  << " Something wrong in getVertexEdge()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for VE: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test vertex triangle relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++) {
    SimplexId triangleNum = triangulation->getVertexTriangleNumber(vertexId);
    SimplexId triangleId;
    for(SimplexId j = 0; j < triangleNum; j++) {
      int result1 = triangulation->getVertexTriangle(vertexId, j, triangleId);
      if(result1) {
        std::cerr << "[TestTopoRelations] vertexId " << vertexId
                  << " Something wrong in getVertexTriangle()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for VF: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test edge triangle relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++) {
    SimplexId triangleNum = triangulation->getEdgeTriangleNumber(edgeId);
    SimplexId triangleId;
    for(SimplexId j = 0; j < triangleNum; j++) {
      int result1 = triangulation->getEdgeTriangle(edgeId, j, triangleId);
      if(result1) {
        std::cout << "[TestTopoRelations] edgeId " << edgeId
                  << " Something wrong in getEdgeTriangle()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for EF: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test triangle edge relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++) {
    SimplexId edgeId;
    for(SimplexId j = 0; j < 3; j++) {
      int result1 = triangulation->getTriangleEdge(triangleId, j, edgeId);
      if(result1) {
        std::cout << "[TestTopoRelations] triangleId " << triangleId
                  << " Something wrong in getTriangleEdge()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for FE: " << t.getElapsedTime()
            << " s." << std::endl;

  SimplexId edgesPerCell = triangulation->getCellEdgeNumber(0);
  SimplexId trianglesPerCell = triangulation->getCellTriangleNumber(0);

  /* test tetrahedron edge relationship */
  t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
  for(SimplexId cellId = 0; cellId < cellNumber; cellId++) {
    SimplexId edgeId;
    for(SimplexId j = 0; j < edgesPerCell; j++) {
      int result1 = triangulation->getCellEdge(cellId, j, edgeId);
      if(result1) {
        std::cout << "[TestTopoRelations] cellId " << cellId
                  << " Something wrong in getCellEdge()! Error code: "
                  << result1 << std::endl;
      }
    }
  }
  std::cout << "[TestTopoRelations] Time usage for TE: " << t.getElapsedTime()
            << " s." << std::endl;

  /* test tetrahedron triangle relationship */
  if(triangulation->getDimensionality() == 3) {
    t.reStart();
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
    for(SimplexId cellId = 0; cellId < cellNumber; cellId++) {
      SimplexId triangleId;
      for(SimplexId j = 0; j < trianglesPerCell; j++) {
        int result1 = triangulation->getCellTriangle(cellId, j, triangleId);
        if(result1) {
          std::cout << "[TestTopoRelations] cellId " << cellId
                    << " Something wrong in getCellTriangle()! Error code: "
                    << result1 << std::endl;
        }
      }
    }
    std::cout << "[TestTopoRelations] Time usage for TF: " << t.getElapsedTime()
              << " s." << std::endl;
  }

    std::cout << "[TestTopoRelations] Data-set (" << vertexNumber
        << " points) processed in " << tot.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
  
  return 0;
}