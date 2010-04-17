/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaximumAPosterioriLeastSquaresDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaximumAPosterioriLeastSquaresDeconvolutionImageFilter_txx
#define __itkMaximumAPosterioriLeastSquaresDeconvolutionImageFilter_txx

#include "itkMaximumAPosterioriLeastSquaresDeconvolutionImageFilter.h"
#include "itkBinaryFunctorWithIndexImageFilter.h"

namespace itk {

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
MaximumAPosterioriLeastSquaresDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::MaximumAPosterioriLeastSquaresDeconvolutionImageFilter()
{
  m_Alpha = 0.001;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
MaximumAPosterioriLeastSquaresDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  ComplexImagePointerType input;
  ComplexImagePointerType psf;
  
  this->Init( input, psf, 0.4f );
  
  typedef itk::BinaryFunctorWithIndexImageFilter< ComplexImageType,
                                    ComplexImageType,
                                    ComplexImageType,
                typename Functor::MaximumAPosterioriLeastSquaresDeconvolution< RegionType, ComplexType > >
                  MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, input );
  mult->SetInput( 1, psf );
  mult->GetFunctor().m_Alpha = m_Alpha;
  // std::cout << input->GetLargestPossibleRegion() << std::endl;
  mult->GetFunctor().m_Region = this->GetPaddedRegion();
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  this->RegisterInternalFilter( mult, 0.1f );
  
  this->End( mult->GetOutput(), 0.3f );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
MaximumAPosterioriLeastSquaresDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: "  << m_Alpha << std::endl;
}

}// end namespace itk
#endif
