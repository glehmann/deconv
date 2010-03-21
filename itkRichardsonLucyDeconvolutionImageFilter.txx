/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRichardsonLucyDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRichardsonLucyDeconvolutionImageFilter_txx
#define __itkRichardsonLucyDeconvolutionImageFilter_txx

#include "itkRichardsonLucyDeconvolutionImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkFlipImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkRelativeChangeCalculator.h"

namespace itk {

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::RichardsonLucyDeconvolutionImageFilter()
{
  m_Normalize = true;
  m_GreatestPrimeFactor = 13;
  m_PadMethod = ZERO_FLUX_NEUMANN;
  m_NumberOfIterations = 10;
  m_RelativeChangeThreshold = 0;
  m_SmoothingFilter = NULL;
  m_Iteration = 0;
  m_RelativeChange = 0;
  m_SmoothingPeriod = 1;
  this->SetNumberOfRequiredInputs(2);
}

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void 
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput());
  InputImageType * input1 = const_cast<PointSpreadFunctionType *>(this->GetPointSpreadFunction());
  if ( !input0 || !input1 )
    { 
    return;
    }
  input0->SetRequestedRegion( input0->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( input1->GetLargestPossibleRegion() );
}


template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  // members used to monitor the iterations
  m_Iteration = 0;
  m_RelativeChange = 0;
  
  this->AllocateOutputs();
  const InputImageType * input = this->GetInput();
  const PointSpreadFunctionType * kernel = this->GetPointSpreadFunction();
  OutputImageType * output = this->GetOutput();

  // Create a process accumulator for tracking the progress of this minipipeline
  // ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  // progress->SetMiniPipelineFilter(this);

  typedef itk::FlipImageFilter< PointSpreadFunctionType > FlipType;
  typename FlipType::Pointer flip = FlipType::New();
  flip->SetInput( kernel );
  typename FlipType::FlipAxesArrayType axes;
  axes.Fill( true ); // we must flip all the axes
  flip->SetFlipAxes( axes );
  flip->SetNumberOfThreads( this->GetNumberOfThreads() );
  flip->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( flip, 0.005f );

  typedef itk::ImageToImageFilter< PointSpreadFunctionType, InternalImageType > NormType;
  typedef itk::NormalizeToConstantImageFilter< PointSpreadFunctionType, InternalImageType > NormConstType;
  typedef itk::CastImageFilter< PointSpreadFunctionType, InternalImageType > CastType;
  typename NormType::Pointer norm;
  if( m_Normalize )
    {
    norm = NormConstType::New();
    }
  else
    {
    norm = CastType::New();
    }
  norm->SetInput( flip->GetOutput() );
  norm->SetNumberOfThreads( this->GetNumberOfThreads() );
  norm->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( norm, 0.005f );

  typedef itk::FFTPadImageFilter< InputImageType, InternalImageType, InternalImageType, InternalImageType > PadType;
  typename PadType::Pointer pad = PadType::New();
  pad->SetInput( input );
  pad->SetInputKernel( norm->GetOutput() );
  pad->SetNumberOfThreads( this->GetNumberOfThreads() );
  // the padded image is reused several time, so it shouldn't be released
  // pad->SetReleaseDataFlag( true );
  pad->SetPadMethod( m_PadMethod );
  // progress->RegisterInternalFilter( pad, 0.05f );

  typedef itk::FFTShiftImageFilter< InternalImageType, InternalImageType > ShiftType;
  typename ShiftType::Pointer shift = ShiftType::New();
  shift->SetInput( pad->GetOutputKernel() );
  shift->SetInverse( true );
  shift->SetNumberOfThreads( this->GetNumberOfThreads() );
  shift->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( shift, 0.04f );
  
  typedef itk::FFTRealToComplexConjugateImageFilter< InternalPrecisionType, ImageDimension > FFTType;
  typename FFTType::Pointer fftk = FFTType::New();
  fftk->SetInput( shift->GetOutput() );
  fftk->SetNumberOfThreads( this->GetNumberOfThreads() );
  // don't release the data here: we must keep it for the iterations
  // fftk->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( fftk, 0.25f );
  
  // vnl filters need a size which is a power of 2
  if( std::string(fftk->GetNameOfClass()).find("Vnl") == 0 )
    {
    pad->SetPadToPowerOfTwo( true );
    // fake the cyclic behavior
    if( m_PadMethod == NO_PADDING )
      {
      pad->SetPadMethod( WRAP );
      }
    }
  else
    {
    pad->SetGreatestPrimeFactor( m_GreatestPrimeFactor );
    }


  // iterated code from here
  
  // first convolve the input image by the psf
  
  typename FFTType::Pointer fft = FFTType::New();
  fft->SetInput( pad->GetOutput() );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( fft, 0.25f );

  typedef typename FFTType::OutputImagePixelType ComplexType;
  
  typedef itk::MultiplyImageFilter< typename FFTType::OutputImageType,
                typename FFTType::OutputImageType,
                typename FFTType::OutputImageType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, fftk->GetOutput() );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  // progress->RegisterInternalFilter( mult, 0.1f );
  
  typedef itk::FFTComplexConjugateToRealImageFilter< InternalPrecisionType, ImageDimension > IFFTType;
  typename IFFTType::Pointer ifft = IFFTType::New();
  ifft->SetInput( mult->GetOutput() );
  // we can't run a single update here: we have to set the
  // ActualXDimensionIsOdd attribute, which can be done only after the update of
  // output information of the pad filter
  pad->UpdateOutputInformation();
  ifft->SetActualXDimensionIsOdd( pad->GetOutput()->GetLargestPossibleRegion().GetSize()[0] % 2 );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( ifft, 0.25f );

  // input convolution completed
  // divide the input by (the convolved image + epsilon)

  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::EpsilonDivide< TInternalPrecision > >
                  EpsilonDivideType;
  typename EpsilonDivideType::Pointer ediv = EpsilonDivideType::New();
  ediv->SetInput( 1, pad->GetOutput() );
  ediv->SetInput( 0, ifft->GetOutput() );
  ediv->SetNumberOfThreads( this->GetNumberOfThreads() );
  ediv->SetReleaseDataFlag( true );
  ediv->SetInPlace( true );
  // progress->RegisterInternalFilter( ediv, 0.1f );
  
  // convolve the divided image by the psf
  
  typename FFTType::Pointer fft2 = FFTType::New();
  fft2->SetInput( ediv->GetOutput() );
  fft2->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft2->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( fft, 0.25f );

  typedef itk::BinaryFunctorImageFilter< typename FFTType::OutputImageType,
                typename FFTType::OutputImageType,
                typename FFTType::OutputImageType,
                typename Functor::MultiplyByComplexConjugate< ComplexType > >
                  MultiplyByComplexConjugateType;
  typename MultiplyByComplexConjugateType::Pointer cmult = MultiplyByComplexConjugateType::New();
  cmult->SetInput( 0, fft2->GetOutput() );
  cmult->SetInput( 1, fftk->GetOutput() );
  cmult->SetNumberOfThreads( this->GetNumberOfThreads() );
  cmult->SetReleaseDataFlag( true );
  cmult->SetInPlace( true );
  // progress->RegisterInternalFilter( cmult, 0.1f );
  
  typename IFFTType::Pointer ifft2 = IFFTType::New();
  ifft2->SetInput( cmult->GetOutput() );
  // we can't run a single update here: we have to set the
  // ActualXDimensionIsOdd attribute, which can be done only after the update of
  // output information of the pad filter
  pad->UpdateOutputInformation();
  ifft2->SetActualXDimensionIsOdd( pad->GetOutput()->GetLargestPossibleRegion().GetSize()[0] % 2 );
  ifft2->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft2->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( ifft, 0.25f );
  
  // divided image convolution completed
  // multiply the result with the input
  
  typedef itk::MultiplyImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > RealMultType;
  typename RealMultType::Pointer rmult = RealMultType::New();
  rmult->SetInput( 0, ifft2->GetOutput() );
  rmult->SetInput( 1, pad->GetOutput() );
  rmult->SetNumberOfThreads( this->GetNumberOfThreads() );
  // can't be released, it is required by two filters on next iteration
  // rmult->SetReleaseDataFlag( true );
  rmult->SetInPlace( true );
  // progress->RegisterInternalFilter( rmult, 0.1f );
  
  // begin the iterations
  typedef typename itk::RelativeChangeCalculator< InternalImageType > ChangeType;
  typename ChangeType::Pointer change = ChangeType::New();
  typename InternalImageType::Pointer img = pad->GetOutput();
  for( m_Iteration=1; m_Iteration<=m_NumberOfIterations; m_Iteration++ )
    {
    // should we use regularisation filter? -- tested in the iteration on purpose, to be able to
    // change the filter by looking at the iteration event
    typename SmoothingFilterType::Pointer last = rmult.GetPointer();
    if( m_SmoothingFilter.IsNotNull() && m_Iteration % m_SmoothingPeriod == 0 )
      {
      m_SmoothingFilter->SetInput( rmult->GetOutput() );
      last = m_SmoothingFilter;
      }
    last->Update();
    
    // do we have to stop the iterations based on the relative change?
    change->SetImage( img );
    change->SetNewImage( last->GetOutput() );
    change->Compute();
    // std::cout << change->GetOutput() << std::endl;
    m_RelativeChange = change->GetOutput();
    if( m_RelativeChangeThreshold > 0 && m_RelativeChange < m_RelativeChangeThreshold )
      {
      break;
      }
    else
      {
      // ok, lets go for another round
      img = last->GetOutput();
      img->DisconnectPipeline();
      fft->SetInput( img );
      rmult->SetInput( 1, img );
      this->UpdateProgress( m_Iteration/(float)m_NumberOfIterations );
      this->InvokeEvent( IterationEvent() );
      }
    }
  img->SetReleaseDataFlag( true );
  this->UpdateProgress( 1.0 );
  m_Iteration--; // to have the real number of iterations
   
  typedef itk::IntensityWindowingImageFilter< InternalImageType, OutputImageType > WindowType;
  typename WindowType::Pointer window = WindowType::New();
  window->SetInput( img );
  window->SetWindowMinimum( NumericTraits< OutputImagePixelType >::Zero );
  window->SetWindowMaximum( NumericTraits< OutputImagePixelType >::max() );
  window->SetNumberOfThreads( this->GetNumberOfThreads() );
  window->SetReleaseDataFlag( true );
  window->SetInPlace( true );
  // progress->RegisterInternalFilter( window, 0.025f );
  
  typedef itk::RegionFromReferenceImageFilter< OutputImageType, OutputImageType > CropType;
  typename CropType::Pointer crop = CropType::New();
  crop->SetInput( window->GetOutput() );
  crop->SetReferenceImage( input );
  crop->SetNumberOfThreads( this->GetNumberOfThreads() );
  // progress->RegisterInternalFilter( crop, 0.025f );
  
  crop->GraftOutput( output );
  crop->Update();
  this->GraftOutput( crop->GetOutput() );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Iteration: "  << m_Iteration << std::endl;
  os << indent << "NumberOfIterations: "  << m_NumberOfIterations << std::endl;
  os << indent << "RelativeChange: "  << m_RelativeChange << std::endl;
  os << indent << "RelativeChangeThreshold: "  << m_RelativeChangeThreshold << std::endl;
  os << indent << "SmoothingPeriod: "  << m_SmoothingPeriod << std::endl;
//  os << indent << "SmoothingFilter: "  << m_SmoothingFilter << std::endl;
  os << indent << "SmoothingFilter: ";
  if( m_SmoothingFilter.IsNull() )
    {
    std::cout << "NULL" << std::endl;
    }
  else
    {
    m_SmoothingFilter->Print( os, indent.GetNextIndent() );
    }
  os << indent << "Normalize: "  << m_Normalize << std::endl;
  os << indent << "GreatestPrimeFactor: "  << m_GreatestPrimeFactor << std::endl;
  os << indent << "PadMethod: "  << m_PadMethod << std::endl;
}

}// end namespace itk
#endif
