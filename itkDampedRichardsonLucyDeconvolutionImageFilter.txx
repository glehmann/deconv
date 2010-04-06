/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDampedRichardsonLucyDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDampedRichardsonLucyDeconvolutionImageFilter_txx
#define __itkDampedRichardsonLucyDeconvolutionImageFilter_txx

#include "itkDampedRichardsonLucyDeconvolutionImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyByComplexConjugateImageFilter.h"
#include "itkRelativeChangeCalculator.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
DampedRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::DampedRichardsonLucyDeconvolutionImageFilter()
{
  m_Threshold = 2;
}


template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
DampedRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Init()
{
  Superclass::Init();
  
  // replace the simple division by our more complex implementation
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::DampedRichardsonLucy< TInternalPrecision > >
                  DampedRichardsonLucyType;
  typename DampedRichardsonLucyType::Pointer div = DampedRichardsonLucyType::New();
  div->SetInput( 0, this->m_Division->GetInput( 0 ) );
  div->SetInput( 1, this->m_Division->GetInput( 1 ) );
  div->GetFunctor().m_T = m_Threshold;
  div->GetFunctor().m_N = 10;
  div->SetNumberOfThreads( this->GetNumberOfThreads() );
  div->SetReleaseDataFlag( true );
  div->SetInPlace( true );
  
  this->m_Convolution2->SetInput( div->GetOutput() );
  this->m_Division = div;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
DampedRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Threshold: "  << m_Threshold << std::endl;
}

}// end namespace itk
#endif
