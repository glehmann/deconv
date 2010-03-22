/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkJanssonVanCittertDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkJanssonVanCittertDeconvolutionImageFilter_txx
#define __itkJanssonVanCittertDeconvolutionImageFilter_txx

#include "itkJanssonVanCittertDeconvolutionImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
JanssonVanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::InitFunctor( FunctorType & functor )
{
  Superclass::InitFunctor( functor );
  
  typedef itk::MinimumMaximumImageCalculator< TInputImage > CalculatorType;
  typename CalculatorType::Pointer calc = CalculatorType::New();
  calc->SetImage( this->GetInput() );
  calc->Compute();

  functor.m_A = (TInternalPrecision) calc->GetMaximum() / 2.0;
  // std::cout << "functor.m_A: " << functor.m_A << std::endl;
}


}// end namespace itk
#endif
