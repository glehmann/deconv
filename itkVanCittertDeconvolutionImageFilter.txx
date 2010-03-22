/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVanCittertDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVanCittertDeconvolutionImageFilter_txx
#define __itkVanCittertDeconvolutionImageFilter_txx

#include "itkVanCittertDeconvolutionImageFilter.h"
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
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRelativeChangeCalculator.h"

namespace itk {

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::VanCittertDeconvolutionImageFilter()
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
  m_RegularizationFilter = NULL;
  m_Alpha = 1;
  m_NonNegativity = true;
  this->SetNumberOfRequiredInputs(2);
}

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void 
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
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
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
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
  
  // compute the residual
  typedef itk::SubtractImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > SubtractType;
  typename SubtractType::Pointer sub = SubtractType::New();
  sub->SetInput( 0, pad->GetOutput() );
  sub->SetInput( 1, ifft->GetOutput() );
  sub->SetNumberOfThreads( this->GetNumberOfThreads() );
  sub->SetReleaseDataFlag( true );
  // don't run in place - we need to keep the input image
  // sub->SetInPlace( true );

  if( m_RegularizationFilter.IsNotNull() )
    {
    // connect the regularization filter
    m_RegularizationFilter->SetInput( sub->GetOutput() );
    }
  
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::VanCittert< TInternalPrecision > >
                  VanCittertType;
  typename VanCittertType::Pointer add = VanCittertType::New();
  add->SetInput( 1, pad->GetOutput() );
  if( m_RegularizationFilter.IsNull() )
    {
    add->SetInput( 0, sub->GetOutput() );
    }
  else
    {
    add->SetInput( 0, m_RegularizationFilter->GetOutput() );
    }
  add->GetFunctor().m_Alpha = m_Alpha;
  add->GetFunctor().m_NonNegativity = m_NonNegativity;
  add->SetNumberOfThreads( this->GetNumberOfThreads() );
  add->SetReleaseDataFlag( true );
  add->SetInPlace( true );
  // progress->RegisterInternalFilter( add, 0.1f );
  
  // begin the iterations
  typedef typename itk::RelativeChangeCalculator< InternalImageType > ChangeType;
  typename ChangeType::Pointer change = ChangeType::New();
  typename InternalImageType::Pointer img = pad->GetOutput();
  for( m_Iteration=1; m_Iteration<=m_NumberOfIterations; m_Iteration++ )
    {
    // should we use smoothing filter? -- tested in the iteration on purpose, to be able to
    // change the filter by looking at the iteration event
    typename InternalFilterType::Pointer last = add.GetPointer();
    if( m_SmoothingFilter.IsNotNull() && m_Iteration % m_SmoothingPeriod == 0 )
      {
      m_SmoothingFilter->SetInput( add->GetOutput() );
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
      add->SetInput( 1, img );
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
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Iteration: "  << m_Iteration << std::endl;
  os << indent << "NumberOfIterations: "  << m_NumberOfIterations << std::endl;
  os << indent << "Alpha: "  << m_Alpha << std::endl;
  os << indent << "NonNegativity: "  << m_NonNegativity << std::endl;
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
  os << indent << "RegularizationFilter: ";
  if( m_RegularizationFilter.IsNull() )
    {
    std::cout << "NULL" << std::endl;
    }
  else
    {
    m_RegularizationFilter->Print( os, indent.GetNextIndent() );
    }
  os << indent << "Normalize: "  << m_Normalize << std::endl;
  os << indent << "GreatestPrimeFactor: "  << m_GreatestPrimeFactor << std::endl;
  os << indent << "PadMethod: "  << m_PadMethod << std::endl;
}

}// end namespace itk
#endif
