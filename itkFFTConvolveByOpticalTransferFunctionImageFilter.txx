/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolveByOpticalTransferFunctionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolveByOpticalTransferFunctionImageFilter_txx
#define __itkFFTConvolveByOpticalTransferFunctionImageFilter_txx

#include "itkFFTConvolveByOpticalTransferFunctionImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMultiplyByComplexConjugateImageFilter.h"

namespace itk {

template <class TPixel, unsigned int NDimension>
FFTConvolveByOpticalTransferFunctionImageFilter<TPixel, NDimension>
::FFTConvolveByOpticalTransferFunctionImageFilter()
{
  m_Transpose = false;
  this->SetNumberOfRequiredInputs(2);
}

template <class TPixel, unsigned int NDimension>
void 
FFTConvolveByOpticalTransferFunctionImageFilter<TPixel, NDimension>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput());
  OpticalTransferFunctionType * input1 = const_cast<OpticalTransferFunctionType *>(this->GetOpticalTransferFunction());
  if ( !input0 || !input1 )
    { 
    return;
    }
  input0->SetRequestedRegion( input0->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( input1->GetLargestPossibleRegion() );
}


template <class TPixel, unsigned int NDimension>
void
FFTConvolveByOpticalTransferFunctionImageFilter<TPixel, NDimension>
::GenerateData()
{
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter( this );

  // Allocate the output
  this->AllocateOutputs();
  
  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  fft->SetInput( this->GetInput() );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter(fft, .45f);
  
  typedef MultiplyImageFilter< OpticalTransferFunctionType,
                                OpticalTransferFunctionType,
                                OpticalTransferFunctionType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, this->GetOpticalTransferFunction() );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  
  typedef MultiplyByComplexConjugateImageFilter< OpticalTransferFunctionType > ConjMultType;
  typename ConjMultType::Pointer cmult = ConjMultType::New();
  cmult->SetInput( fft->GetOutput() );
  cmult->SetInput( 1, this->GetOpticalTransferFunction() );
  cmult->SetNumberOfThreads( this->GetNumberOfThreads() );
  cmult->SetReleaseDataFlag( true );
  cmult->SetInPlace( true );
  
  typename IFFTFilterType::Pointer ifft = IFFTFilterType::New();
  if( m_Transpose )
    {
    ifft->SetInput( cmult->GetOutput() );
    }
  else
    {
    ifft->SetInput( mult->GetOutput() );
    }
  ifft->SetActualXDimensionIsOdd( this->GetOutput()->GetLargestPossibleRegion().GetSize()[0] % 2 );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter(ifft, .45f);  

  ifft->GraftOutput( this->GetOutput() );
  ifft->Update();
  this->GraftOutput( ifft->GetOutput() );
}


template <class TPixel, unsigned int NDimension>
void
FFTConvolveByOpticalTransferFunctionImageFilter<TPixel, NDimension>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Transpose: "  << m_Transpose << std::endl;
}

}// end namespace itk
#endif
