/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPoissonMaximumAPosterioriDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPoissonMaximumAPosterioriDeconvolutionImageFilter_txx
#define __itkPoissonMaximumAPosterioriDeconvolutionImageFilter_txx

#include "itkPoissonMaximumAPosterioriDeconvolutionImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkFFTConvolveByOpticalTransferFunctionImageFilter.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
PoissonMaximumAPosterioriDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Init()
{
  Superclass::Init();

  InternalImagePointerType input;
  InternalImagePointerType firstEstimate;
  ComplexImagePointerType psf;
  
  this->PrepareInputs( input, firstEstimate, psf );

  // first convolve the input image by the psf
  typedef itk::FFTConvolveByOpticalTransferFunctionImageFilter< InternalPrecisionType, ImageDimension > ConvolutionType;
  typename ConvolutionType::Pointer conv = ConvolutionType::New();
  // conv->SetInput( input );
  conv->SetOpticalTransferFunction( psf );
  conv->SetNumberOfThreads( this->GetNumberOfThreads() );
  conv->SetReleaseDataFlag( true );
  m_Convolution = conv;
  
  // divide the input by (the convolved image + epsilon)
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::PoissonMaximumAPosteriori< TInternalPrecision > >
                  PoissonMaximumAPosterioriType;
  typename PoissonMaximumAPosterioriType::Pointer div = PoissonMaximumAPosterioriType::New();
  div->SetInput( 1, input );
  div->SetInput( 0, conv->GetOutput() );
  div->SetNumberOfThreads( this->GetNumberOfThreads() );
  div->SetReleaseDataFlag( true );
  div->SetInPlace( true );
  m_Division = div;
  
  // convolve the divided image by the psf
  typename ConvolutionType::Pointer conv2 = ConvolutionType::New();
  conv2->SetInput( div->GetOutput() );
  conv2->SetOpticalTransferFunction( psf );
  conv2->SetTranspose( true );
  conv2->SetNumberOfThreads( this->GetNumberOfThreads() );
  conv2->SetReleaseDataFlag( true );
  m_Convolution2 = conv2;

  // multiply the result with the input
  
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::PoissonMaximumAPosteriori2< TInternalPrecision > >
                  RealMultType;
  typename RealMultType::Pointer mult = RealMultType::New();
  mult->SetInput( 0, conv2->GetOutput() );
  // mult->SetInput( 1, input );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  // can't be released, it is required by two filters on next iteration
  // mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  m_Multiplication = mult;  
  
  SetEstimate( firstEstimate );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
PoissonMaximumAPosterioriDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::SetEstimate( InternalImageType * estimate )
{
  Superclass::SetEstimate( estimate );
  m_Convolution->SetInput( estimate );
  m_Multiplication->SetInput( 1, estimate );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
typename PoissonMaximumAPosterioriDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>::InternalImageType::Pointer
PoissonMaximumAPosterioriDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::NewEstimate()
{
  return m_Multiplication->GetOutput();
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
PoissonMaximumAPosterioriDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::End()
{
  Superclass::End();

  m_Convolution = NULL;
  m_Division = NULL;
  m_Convolution2 = NULL;
  m_Multiplication = NULL;
}

}// end namespace itk
#endif
