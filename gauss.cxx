#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTestingMacros.h"


int main(int argc, char * argv[])
{

  if( argc < 3 )
    {
    std::cerr << "usage: " << argv[0] << " output sigma" << std::endl;
    std::cerr << " output: the output image" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef float RPType;
  typedef itk::Image< RPType, dim > RIType;
  
  RIType::Pointer img = RIType::New();
  RIType::SizeType s;
  s.Fill( atoi(argv[2]) * 3 * 2 + 1 );
  img->SetRegions( s );
  img->Allocate();
  img->FillBuffer( 0 );
  RIType::IndexType idx;
  idx.Fill( atoi(argv[2]) * 3 );
  img->SetPixel( idx, 1.0 );

  typedef itk::SmoothingRecursiveGaussianImageFilter< RIType, RIType > GaussType;
  GaussType::Pointer gauss = GaussType::New();
  gauss->SetInput( img );
  gauss->SetSigma( atoi(argv[2]) );

  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::RescaleIntensityImageFilter< RIType, IType > RescaleType;
  RescaleType::Pointer rescale = RescaleType::New();
  rescale->SetInput( gauss->GetOutput() );
  
  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rescale->GetOutput() );
  writer->SetFileName( argv[1] );
  
  writer->Update();

  return 0;
}
