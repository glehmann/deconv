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
::Init()
{
  Superclass::Init();
  
  typedef itk::MinimumMaximumImageCalculator< TInputImage > CalculatorType;
  typename CalculatorType::Pointer calc = CalculatorType::New();
  calc->SetImage( this->GetInput() );
  calc->Compute();

  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                Functor::JanssonVanCittert<TInternalPrecision> >
                  VanCittertType;
  typename VanCittertType::Pointer add = VanCittertType::New();
  add->SetInput( 1, this->m_Add->GetInput( 1 ) );
  add->SetInput( 0, this->m_Add->GetInput( 0 ) );
  add->GetFunctor().m_Alpha = this->GetAlpha();
  add->GetFunctor().m_NonNegativity = this->GetNonNegativity();
  add->GetFunctor().m_B2 = (TInternalPrecision) calc->GetMaximum() / 2.0;
  // std::cout << "functor.m_B2: " << functor.m_B2 << std::endl;
  add->SetNumberOfThreads( this->GetNumberOfThreads() );
  add->SetReleaseDataFlag( true );
  add->SetInPlace( true );
  
  this->m_Add = add;
}


}// end namespace itk
#endif
