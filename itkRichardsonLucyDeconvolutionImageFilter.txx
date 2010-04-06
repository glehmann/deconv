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
#include "itkFFTConvolveByOpticalTransferFunctionImageFilter.h"
#include "itkDivideOrZeroOutImageFilter.h"
#include "itkMultiplyImageFilter.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Init()
{
  Superclass::Init();

  InternalImagePointerType input;
  ComplexImagePointerType psf;
  
  this->Superclass::Superclass::Init( input, psf, 0 );

  // first convolve the input image by the psf
  typedef itk::FFTConvolveByOpticalTransferFunctionImageFilter< InternalPrecisionType, ImageDimension > ConvolutionType;
  typename ConvolutionType::Pointer conv = ConvolutionType::New();
  // conv->SetInput( input );
  conv->SetOpticalTransferFunction( psf );
  conv->SetNumberOfThreads( this->GetNumberOfThreads() );
  conv->SetReleaseDataFlag( true );
  m_Convolution = conv;
  
  // divide the input by (the convolved image + epsilon)
  typedef itk::DivideOrZeroOutImageFilter< InternalImageType > DivideType;
  typename DivideType::Pointer div = DivideType::New();
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
  typedef itk::MultiplyImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > RealMultType;
  typename RealMultType::Pointer mult = RealMultType::New();
  mult->SetInput( 0, conv2->GetOutput() );
  // mult->SetInput( 1, input );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  // can't be released, it is required by two filters on next iteration
  // mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  m_Multiplication = mult;  
  
  SetEstimate( input );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::SetEstimate( InternalImageType * estimate )
{
  Superclass::SetEstimate( estimate );
  m_Convolution->SetInput( estimate );
  m_Multiplication->SetInput( 1, estimate );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
typename RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>::InternalImageType::Pointer
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::NewEstimate()
{
  return m_Multiplication->GetOutput();
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
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
