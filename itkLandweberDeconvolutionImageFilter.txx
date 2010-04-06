/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLandweberDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLandweberDeconvolutionImageFilter_txx
#define __itkLandweberDeconvolutionImageFilter_txx

#include "itkLandweberDeconvolutionImageFilter.h"
#include "itkFFTConvolveByOpticalTransferFunctionImageFilter.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
LandweberDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Init()
{
  Superclass::Init();
  
  typedef itk::FFTConvolveByOpticalTransferFunctionImageFilter< InternalPrecisionType, ImageDimension > ConvolutionType;
  
  ConvolutionType * conv = dynamic_cast<ConvolutionType*>( this->m_Convolution.GetPointer() );
  
  typename ConvolutionType::Pointer conv2 = ConvolutionType::New();
  conv2->SetInput( this->m_Subtract->GetOutput() );
  conv2->SetOpticalTransferFunction( conv->GetOpticalTransferFunction() );
  conv2->SetNumberOfThreads( this->GetNumberOfThreads() );
  conv2->SetReleaseDataFlag( true );
  
  this->m_Add->SetInput( 0, conv2->GetOutput() );
  m_Convolution2 = conv2;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
LandweberDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::End()
{
  Superclass::End();
  m_Convolution2 = NULL;
}


}// end namespace itk
#endif
