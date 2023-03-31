#include <ttkTestTopoRelations.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTestTopoRelations);

ttkTestTopoRelations::ttkTestTopoRelations() {
  // init
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  debugLevel_ = 4;
}

ttkTestTopoRelations::~ttkTestTopoRelations() {
}

int ttkTestTopoRelations::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTestTopoRelations::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTestTopoRelations::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  MemoryUsage m;

  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return -1;

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");

  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return -1;

  this->preconditionTriangulation(triangulation);

  int status = -1;
  ttkTemplateMacro(
    triangulation->getType(),
    (status = this->execute((TTK_TT *)triangulation->getData())));

  if(status != 0)
    return -1;

  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  {
    printMsg("Memory usage: " + std::to_string(m.getValue_in_MB(false)) + " MB.");
  }

  return 1;
}
