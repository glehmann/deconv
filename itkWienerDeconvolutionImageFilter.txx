/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWienerDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWienerDeconvolutionImageFilter_txx
#define __itkWienerDeconvolutionImageFilter_txx

#include "itkWienerDeconvolutionImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"

namespace itk {

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
WienerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::WienerDeconvolutionImageFilter()
{
  m_Gamma = 0.001;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
WienerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  ComplexImagePointerType input;
  ComplexImagePointerType psf;
  
  this->PrepareInputs( input, psf, 0.6f );
  
  typedef itk::BinaryFunctorImageFilter< ComplexImageType,
                                    ComplexImageType,
                                    ComplexImageType,
                typename Functor::WienerDeconvolution< ComplexType > >
                  MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, input );
  mult->SetInput( 1, psf );
  mult->GetFunctor().m_Gamma = m_Gamma;
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  this->RegisterInternalFilter( mult, 0.1f );
  
  this->ProduceOutput( mult->GetOutput(), 0.3f );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
WienerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Gamma: "  << m_Gamma << std::endl;
}

}// end namespace itk
#endif
