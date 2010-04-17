/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTikhonovMillerDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTikhonovMillerDeconvolutionImageFilter_txx
#define __itkTikhonovMillerDeconvolutionImageFilter_txx

#include "itkTikhonovMillerDeconvolutionImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "itkLaplacianImageFilter.h"

namespace itk {

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
TikhonovMillerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::TikhonovMillerDeconvolutionImageFilter()
{
  m_Gamma = 0.1;
}

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InternalImageType * input = const_cast<InternalImageType *>(this->GetRegularizationImage());
  if ( !input )
    { 
    return;
    }
  input->SetRequestedRegion( input->GetLargestPossibleRegion() );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  ComplexImagePointerType input;
  ComplexImagePointerType psf;
  ComplexImagePointerType reg;
  
  this->Init( input, psf, 0.4f );
  this->PrepareImage( reg, this->InternalInput( 2 ), true, false, 0.2f );
  if( reg.IsNull() )
    {
    // something went wrong - lets use the default regularization image
    this->PrepareImage( reg, this->Laplacian(), true, false, 0.2f );
    }
  
  typedef itk::TernaryFunctorImageFilter< ComplexImageType,
                                    ComplexImageType,
                                    ComplexImageType,
                                    ComplexImageType,
                typename Functor::TikhonovMillerDeconvolution< ComplexType > >
                  MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, input );
  mult->SetInput( 1, psf );
  mult->SetInput( 2, reg );
  mult->GetFunctor().m_Gamma = m_Gamma;
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  this->RegisterInternalFilter( mult, 0.1f );
  
  this->End( mult->GetOutput(), 0.3f );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
typename TikhonovMillerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>::InternalImagePointerType
TikhonovMillerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Laplacian()
{
  InternalImagePointerType dot = InternalImageType::New();
  dot->SetSpacing( this->GetInput()->GetSpacing() );
  SizeType s;
  s.Fill( 3 );
  dot->SetRegions( s );
  dot->Allocate();
  dot->FillBuffer( 0 );
  IndexType idx;
  idx.Fill( 1 );
  dot->SetPixel( idx, 1 );
  
  typedef LaplacianImageFilter< InternalImageType, InternalImageType > LaplacianType;
  typename LaplacianType::Pointer laplacian = LaplacianType::New();
  laplacian->SetInput( dot );
  laplacian->SetNormalizeToOne( true );
  laplacian->Update();
  return laplacian->GetOutput();
  
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
TikhonovMillerDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Gamma: "  << m_Gamma << std::endl;
}

}// end namespace itk
#endif
