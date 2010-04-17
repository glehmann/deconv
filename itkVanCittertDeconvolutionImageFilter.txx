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
#include "itkBinaryFunctorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkFFTConvolveByOpticalTransferFunctionImageFilter.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::VanCittertDeconvolutionImageFilter()
{
  m_Alpha = 1;
  m_NonNegativity = true;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
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
  
  // compute the residual
  typedef itk::SubtractImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > SubtractType;
  typename SubtractType::Pointer sub = SubtractType::New();
  sub->SetInput( 0, input );
  sub->SetInput( 1, conv->GetOutput() );
  sub->SetNumberOfThreads( this->GetNumberOfThreads() );
  sub->SetReleaseDataFlag( true );
  // don't run in place - we need to keep the input image
  // sub->SetInPlace( true );
  m_Subtract = sub;

  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                Functor::VanCittert<TInternalPrecision> >
                  VanCittertType;
  typename VanCittertType::Pointer add = VanCittertType::New();
  // add->SetInput( 1, input );
  add->SetInput( 0, sub->GetOutput() );
  add->GetFunctor().m_Alpha = m_Alpha;
  add->GetFunctor().m_NonNegativity = m_NonNegativity;
  add->SetNumberOfThreads( this->GetNumberOfThreads() );
  add->SetReleaseDataFlag( true );
  add->SetInPlace( true );
  m_Add = add;
  
  SetEstimate( firstEstimate );

}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::SetEstimate( InternalImageType * estimate )
{
  Superclass::SetEstimate( estimate );
  m_Convolution->SetInput( estimate );
  m_Add->SetInput( 1, estimate );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
typename VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>::InternalImageType::Pointer
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::NewEstimate()
{
  return m_Add->GetOutput();
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::End()
{
  Superclass::End();

  m_Convolution = NULL;
  m_Subtract = NULL;
  m_Add = NULL;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: "  << m_Alpha << std::endl;
  os << indent << "NonNegativity: "  << m_NonNegativity << std::endl;
}

}// end namespace itk
#endif
