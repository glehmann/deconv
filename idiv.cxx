#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkIDivergenceCalculator.h"
#include "itkTestingMacros.h"

int main(int argc, char * argv[])
{

  if( argc < 3 )
    {
    std::cerr << "usage: " << argv[0] << " intput new" << std::endl;
    std::cerr << " input: the input image" << std::endl;
    std::cerr << " output: the output image" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  typedef itk::IDivergenceCalculator< IType > CalculatorType;
  CalculatorType::Pointer idiv = CalculatorType::New();
  idiv->SetImage( reader->GetOutput() );
  idiv->SetNewImage( reader2->GetOutput() );
  idiv->Compute();
  std::cout << "I-Divergence: " << idiv->GetOutput() << std::endl;  

  return 0;
}

