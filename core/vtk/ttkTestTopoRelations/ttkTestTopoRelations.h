#pragma once

// VTK Module
#include <ttkTestTopoRelationsModule.h>

// VTK includes -- to adapt
#include <ttkAlgorithm.h>

// ttk code includes
#include <TestTopoRelations.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKTESTTOPORELATIONS_EXPORT ttkTestTopoRelations
  : public ttkAlgorithm,
    protected ttk::TestTopoRelations {
private:
  std::string ScalarField;

public:
  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  static ttkTestTopoRelations *New();
  vtkTypeMacro(ttkTestTopoRelations, ttkAlgorithm);

protected:
  ttkTestTopoRelations();
  ~ttkTestTopoRelations() override;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};