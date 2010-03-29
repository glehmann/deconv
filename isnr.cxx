#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImprovementInSignalToNoiseRatioCalculator.h"
#include "itkTestingMacros.h"

int main(int argc, char * argv[])
{

  if( argc < 3 )
    {
    std::cerr << "usage: " << argv[0] << " original blurred restored" << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer original = ReaderType::New();
  original->SetFileName( argv[1] );
  original->Update();

  ReaderType::Pointer blurred = ReaderType::New();
  blurred->SetFileName( argv[2] );
  blurred->Update();

  ReaderType::Pointer restored = ReaderType::New();
  restored->SetFileName( argv[3] );
  restored->Update();

  typedef itk::ImprovementInSignalToNoiseRatioCalculator< IType > CalculatorType;
  CalculatorType::Pointer isnr = CalculatorType::New();
  isnr->SetOriginalImage( original->GetOutput() );
  isnr->SetBlurredImage( blurred->GetOutput() );
  isnr->SetRestoredImage( restored->GetOutput() );
  isnr->Compute();
  std::cout << "ISNR: " << isnr->GetOutput() << " dB" << std::endl;  

  return 0;
}

