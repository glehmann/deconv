#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

// #include "itkFFTWRealToComplexConjugateImageFilter.h"
// #include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkMaximumEntropyRichardsonLucyDeconvolutionImageFilter.h"
#include "itkTestingMacros.h"


int main(int argc, char * argv[])
{

  if( argc < 4 )
    {
    std::cerr << "usage: " << argv[0] << " intput kernel output [lambda [nbIt [th [method [gpf]]]]]" << std::endl;
    std::cerr << " input: the input image" << std::endl;
    std::cerr << " output: the output image" << std::endl;
    std::cerr << " method: padding method - 1: zero flux, 2: zero, 3: mirror, 4: wrap" << std::endl;
    std::cerr << " gpf: greatest prime factor of the size of the padded image" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  // itk::SimpleFilterWatcher watcher_reader(reader, "reader");

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  // itk::SimpleFilterWatcher watcher_reader2(reader2, "reader2");

  typedef itk::MaximumEntropyRichardsonLucyDeconvolutionImageFilter< IType > FFTConvolutionType;
  FFTConvolutionType::Pointer conv = FFTConvolutionType::New();
  conv->SetInput( reader->GetOutput() );
  conv->SetPointSpreadFunction( reader2->GetOutput() );
  // test default value
  TEST_SET_GET_VALUE( 13, conv->GetGreatestPrimeFactor() );
  TEST_SET_GET_VALUE( FFTConvolutionType::ZERO_FLUX_NEUMANN, conv->GetPadMethod() );
  TEST_SET_GET_VALUE( 10, conv->GetNumberOfIterations() );
  TEST_SET_GET_VALUE( 0.0, conv->GetRelativeChangeThreshold() );
  // TEST_SET_GET_VALUE( 0.005, conv->GetLambda() );

  if( argc >= 5 )
    {
    conv->SetLambda( atof(argv[4]) );
    // TEST_SET_GET_VALUE( atof(argv[4]), conv->GetLambda() );
    }

  if( argc >= 6 )
    {
    conv->SetNumberOfIterations( atoi(argv[5]) );
    TEST_SET_GET_VALUE( atoi(argv[5]), conv->GetNumberOfIterations() );
    }

  if( argc >= 7 )
    {
    conv->SetRelativeChangeThreshold( atof(argv[6]) );
    TEST_SET_GET_VALUE( atof(argv[6]), conv->GetRelativeChangeThreshold() );
    }

  if( argc >= 8 )
    {
    conv->SetPadMethod( atoi(argv[7]) );
    TEST_SET_GET_VALUE( atoi(argv[7]), conv->GetPadMethod() );
    }

  if( argc >= 9 )
    {
    conv->SetGreatestPrimeFactor( atoi(argv[8]) );
    TEST_SET_GET_VALUE( atoi(argv[8]), conv->GetGreatestPrimeFactor() );
    }
  itk::SimpleFilterWatcher watcher_norm(conv, "conv");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( conv->GetOutput() );
  writer->SetFileName( argv[3] );
  // itk::SimpleFilterWatcher watcher_writer(writer, "writer");
  
  writer->Update();

  return 0;
}