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

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision, class TFunctor>
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision, TFunctor>
::VanCittertDeconvolutionImageFilter()
{
  m_Alpha = 1;
  m_NonNegativity = true;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision, class TFunctor>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision, TFunctor>
::GenerateData()
{
  // members used to monitor the iterations
  this->SetIteration( 0 );
  this->SetRelativeChange( 0.0 );
  
  InternalImagePointerType input;
  ComplexImagePointerType psf;
  bool xIsOdd;
  
  this->Init( input, psf, xIsOdd, 0 );

  // iterated code from here
  
  // first convolve the input image by the psf
  
  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  fft->SetInput( input );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( fft, 0.25f );

  typedef typename FFTFilterType::OutputImagePixelType ComplexType;
  
  typedef itk::MultiplyImageFilter< typename FFTFilterType::OutputImageType,
                typename FFTFilterType::OutputImageType,
                typename FFTFilterType::OutputImageType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, psf );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  // progress->RegisterInternalFilter( mult, 0.1f );
  
  typedef itk::FFTComplexConjugateToRealImageFilter< InternalPrecisionType, ImageDimension > IFFTFilterType;
  typename IFFTFilterType::Pointer ifft = IFFTFilterType::New();
  ifft->SetInput( mult->GetOutput() );
  ifft->SetActualXDimensionIsOdd( xIsOdd );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( ifft, 0.25f );

  // input convolution completed
  
  // compute the residual
  typedef itk::SubtractImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > SubtractType;
  typename SubtractType::Pointer sub = SubtractType::New();
  sub->SetInput( 0, input );
  sub->SetInput( 1, ifft->GetOutput() );
  sub->SetNumberOfThreads( this->GetNumberOfThreads() );
  sub->SetReleaseDataFlag( true );
  // don't run in place - we need to keep the input image
  // sub->SetInPlace( true );

  if( this->GetRegularizationFilter() != NULL )
    {
    // connect the regularization filter
    this->GetRegularizationFilter()->SetInput( sub->GetOutput() );
    }
  
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                FunctorType >
                  VanCittertType;
  typename VanCittertType::Pointer add = VanCittertType::New();
  add->SetInput( 1, input );
  if( this->GetRegularizationFilter() == NULL )
    {
    add->SetInput( 0, sub->GetOutput() );
    }
  else
    {
    add->SetInput( 0, this->GetRegularizationFilter()->GetOutput() );
    }
  this->InitFunctor( add->GetFunctor() );
  add->SetNumberOfThreads( this->GetNumberOfThreads() );
  add->SetReleaseDataFlag( true );
  add->SetInPlace( true );
  // progress->RegisterInternalFilter( add, 0.1f );
  
  // begin the iterations
  typedef typename itk::RelativeChangeCalculator< InternalImageType > ChangeType;
  typename ChangeType::Pointer change = ChangeType::New();
  typename InternalImageType::Pointer img = input;
  for( int i=1; i<=this->GetNumberOfIterations(); i++ )
    {
    this->SetIteration( i );
    // should we use smoothing filter? -- tested in the iteration on purpose, to be able to
    // change the filter by looking at the iteration event
    typename InternalFilterType::Pointer last = add.GetPointer();
    if( this->GetSmoothingFilter() != NULL && i % this->GetSmoothingPeriod() == 0 )
      {
      this->GetSmoothingFilter()->SetInput( add->GetOutput() );
      last = this->GetSmoothingFilter();
      }
    last->Update();
    
    // do we have to stop the iterations based on the relative change?
    change->SetImage( img );
    change->SetNewImage( last->GetOutput() );
    change->Compute();
    // std::cout << change->GetOutput() << std::endl;
    this->SetRelativeChange( change->GetOutput() );
    if( this->GetRelativeChangeThreshold() > 0 && this->GetRelativeChange() < this->GetRelativeChangeThreshold() )
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
      this->UpdateProgress( i/(float)this->GetNumberOfIterations() );
      this->InvokeEvent( IterationEvent() );
      }
    }
  img->SetReleaseDataFlag( true );
  this->UpdateProgress( 1.0 );
   
  this->End( img, 0 );
}


template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision, class TFunctor>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision, TFunctor>
::InitFunctor( FunctorType & functor )
{
  functor.m_Alpha = m_Alpha;
  functor.m_NonNegativity = m_NonNegativity;
  // std::cout << "functor.m_Alpha: " << functor.m_Alpha << std::endl;
  // std::cout << "functor.m_NonNegativity: " << functor.m_NonNegativity << std::endl;
}


template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision, class TFunctor>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision, TFunctor>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: "  << m_Alpha << std::endl;
  os << indent << "NonNegativity: "  << m_NonNegativity << std::endl;
}

}// end namespace itk
#endif
