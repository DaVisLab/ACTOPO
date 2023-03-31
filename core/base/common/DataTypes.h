/// \ingroup base
/// \class ttk::DataTypes
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date May 2018.
///
///\brief TTK base package defining the standard types.

#ifndef _DATATYPES_H
#define _DATATYPES_H

namespace ttk {
  /// \brief Identifier type for simplices of any dimension.
  using LongSimplexId = long long int;

  /// \brief Identifier type for simplices of any dimension.
#ifdef TTK_ENABLE_64BIT_IDS
  using SimplexId = long long int;
#else
  using SimplexId = int;
#endif

  /// \brief Identifier type for threads (i.e. with OpenMP).
  using ThreadId = int;

  /// \brief Identifier type for tasks (i.e. with OpenMP).
  using TaskId = int;

  /// default name for mask scalar field
  const char MaskScalarFieldName[] = "ttkMaskScalarField";

  /// default name for vertex scalar field
  const char VertexScalarFieldName[] = "ttkVertexScalarField";

  /// default name for offset scalar field
  const char OffsetScalarFieldName[] = "ttkOffsetScalarField";

  /// default name for bivariate offset fields
  const char OffsetFieldUName[] = "ttkOffsetFieldU";
  const char OffsetFieldVName[] = "ttkOffsetFieldV";

  // default names for the Morse-Smale complex
  const char MorseSmaleCellDimensionName[] = "CellDimension";
  const char MorseSmaleCellIdName[] = "CellId";
  const char MorseSmaleBoundaryName[] = "IsOnBoundary";
  const char MorseSmaleManifoldSizeName[] = "ManifoldSize";
  const char MorseSmaleSourceIdName[] = "SourceId";
  const char MorseSmaleDestinationIdName[] = "DestinationId";
  const char MorseSmaleSeparatrixIdName[] = "SeparatrixId";
  const char MorseSmaleSeparatrixTypeName[] = "SeparatrixType";
  const char MorseSmaleSeparatrixMaximumName[] = "SeparatrixFunctionMaximum";
  const char MorseSmaleSeparatrixMinimumName[] = "SeparatrixFunctionMinimum";
  const char MorseSmaleSeparatrixDifferenceName[]
    = "SeparatrixFunctionDifference";
  const char MorseSmaleCriticalPointsOnBoundaryName[]
    = "NumberOfCriticalPointsOnBoundary";
  const char MorseSmaleAscendingName[] = "AscendingManifold";
  const char MorseSmaleDescendingName[] = "DescendingManifold";
  const char MorseSmaleManifoldName[] = "MorseSmaleManifold";

  // default names for persistence diagram meta data
  const char PersistenceCriticalTypeName[] = "CriticalType";
  const char PersistenceBirthName[] = "Birth";
  const char PersistenceDeathName[] = "Death";
  const char PersistenceCoordinatesName[] = "Coordinates";
  const char PersistencePairIdentifierName[] = "PairIdentifier";
  const char PersistenceName[] = "Persistence";
  const char PersistencePairTypeName[] = "PairType";

  // default name for compact triangulation index
  const char compactTriangulationIndex[] = "_index";

  /// default value for critical index
  enum class CriticalType {
    Local_minimum = 0,
    Saddle1,
    Saddle2,
    Local_maximum,
    Degenerate,
    Regular
  };
  /// number of different critical types
  const int CriticalTypeNumber = 6;

  /// relation names for TopoClsuter
  enum class RelationType {
    BlankRelation, // Empty -- 0
    VVRelation, // vertexNeighbor -- 1
    VTRelation, // vertexStar -- 2
    EVRelation, // edgeVertex -- 3
    ETRelation, // edgeStar -- 4
    FVRelation, // triangleVertex -- 5
    FTRelation, // triangleStar -- 6
    TVRelation, // cellVertex -- 7
    TTRelation, // cellNeighbor -- 8
    // edges and triangles computation
    IEList, // internalEdgeList -- 9
    IFList, // internalTriangleList -- 10
    IEFList, // intenralEdgeList and internalTriangleList -- 11
    IBCList, // internalBoundaryCellList -- 12
    // require external maps
    VERelation, // vertexEdge -- 13
    VFRelation, // vertexTriangle -- 14
    EFRelation, // edgeTriangle -- 15
    FERelation, // triangleEdge -- 16
    TERelation, // cellEdge -- 17
    TFRelation, // cellTriangle -- 18
    VLRelation, // vertexLink -- 19
    ELRelation, // edgeLink -- 20
    FLRelation, // triangleLink -- 21
    BVRelation, // boundaryVertex -- 22
    BERelation, // boundaryEdge -- 23
    BFRelation, // boundaryTriangle -- 24
    RelationNum // total number of relations
  };

} // namespace ttk

#endif // _DATATYPES_H
