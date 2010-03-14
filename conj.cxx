#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkComplexConjugateImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkTestingMacros.h"


int main(int argc, char * argv[])
{

  if( argc != 3 )
    {
    std::cerr << "usage: " << argv[0] << " intput output" << std::endl;
    std::cerr << " input: the input image" << std::endl;
    std::cerr << " output: the output image" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef float PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  // itk::SimpleFilterWatcher watcher_reader(reader, "reader");

  typedef itk::FFTRealToComplexConjugateImageFilter< PType, 2 > FFTType;
  FFTType::Pointer fft = FFTType::New();
  fft->SetInput( reader->GetOutput() );
  
  typedef itk::ComplexConjugateImageFilter< FFTType::OutputImageType > ConjType;
  ConjType::Pointer conj = ConjType::New();
  conj->SetInput( fft->GetOutput() );
  itk::SimpleFilterWatcher watcher_norm(conj, "conj");

  typedef itk::FFTComplexConjugateToRealImageFilter< PType, 2 > IFFTType;
  IFFTType::Pointer ifft = IFFTType::New();
  ifft->SetInput( conj->GetOutput() );
  
  typedef itk::Image< unsigned char, 2 > OutputImageType;
  typedef itk::CastImageFilter< IType, OutputImageType > CastType;
  CastType::Pointer cast = CastType::New();
  cast->SetInput( ifft->GetOutput() );
  
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( cast->GetOutput() );
  writer->SetFileName( argv[2] );
  // itk::SimpleFilterWatcher watcher_writer(writer, "writer");
  
  writer->Update();

  return 0;
}
