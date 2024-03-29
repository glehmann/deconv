/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter_txx
#define __itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter_txx

#include "itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "itkLaplacianImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyByComplexConjugateImageFilter.h"
#include "itkRelativeChangeCalculator.h"
#include "itkDivideOrZeroOutImageFilter.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
TikhonovMillerRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::TikhonovMillerRichardsonLucyDeconvolutionImageFilter()
{
  m_Lambda = 0.0001;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Init()
{
  Superclass::Init();
  
  // TODO: check the execution speed to see if it can be done faster in the
  // frequency domain -- we already have the input image, so we only have to precompute
  // the laplacian kernel (done only one time, at initialization) and do the mulitplication
  // and the ifft at each iteration.
  typedef LaplacianImageFilter< InternalImageType, InternalImageType > LaplacianType;
  typename LaplacianType::Pointer laplacian = LaplacianType::New();
  laplacian->SetInput( this->m_Convolution->GetInput() );  // take the estimate from another filter
  laplacian->SetNormalizeToOne( true );
  laplacian->SetNumberOfThreads( this->GetNumberOfThreads() );
  laplacian->SetReleaseDataFlag( true );
  this->m_Laplacian = laplacian;
  
  // multiply the result with the input
  typedef itk::TernaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::TikhonovMillerRegularization< TInternalPrecision > > RealMultType;
  typename RealMultType::Pointer mult = RealMultType::New();
  mult->SetInput( 0, this->m_Multiplication->GetInput( 0 ) );
  mult->SetInput( 1, this->m_Multiplication->GetInput( 1 ) );
  mult->SetInput( 2, laplacian->GetOutput() );
  mult->GetFunctor().m_Lambda = m_Lambda;
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  // can't be released, it is required by two filters on next iteration
  // mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  this->m_Multiplication = mult;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::SetEstimate( InternalImageType * estimate )
{
  Superclass::SetEstimate( estimate );
  // make sure that we are not trying to set the input of an uninitialized filter
  if( m_Laplacian.IsNotNull() )
    {
    m_Laplacian->SetInput( estimate );
    }
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::End()
{
  Superclass::End();
  m_Laplacian = NULL;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Lambda: "  << m_Lambda << std::endl;
}

}// end namespace itk
#endif
