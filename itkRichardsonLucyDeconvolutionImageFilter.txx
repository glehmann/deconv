/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRichardsonLucyDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRichardsonLucyDeconvolutionImageFilter_txx
#define __itkRichardsonLucyDeconvolutionImageFilter_txx

#include "itkRichardsonLucyDeconvolutionImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkFlipImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRelativeChangeCalculator.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  // members used to monitor the iterations
  this->SetIteration( 0 );
  this->SetRelativeChange( 0.0 );
  
  InternalImagePointerType input;
  ComplexImagePointerType psf;
  bool xIsOdd;
  
  this->Init( input, psf, xIsOdd, 0 );

  // iterated code from here
  
  // first convolve the input image by the psf
  
  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  fft->SetInput( input );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );

  typedef typename FFTFilterType::OutputImagePixelType ComplexType;
  
  typedef itk::MultiplyImageFilter< ComplexImageType,
                ComplexImageType,
                ComplexImageType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, psf );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  
  typedef itk::FFTComplexConjugateToRealImageFilter< InternalPrecisionType, ImageDimension > IFFTFilterType;
  typename IFFTFilterType::Pointer ifft = IFFTFilterType::New();
  ifft->SetInput( mult->GetOutput() );
  ifft->SetActualXDimensionIsOdd( xIsOdd );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft->SetReleaseDataFlag( true );

  // input convolution completed
  
  // do we have to do some work on the residual?
  typedef itk::SubtractImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > SubtractType;
  typename SubtractType::Pointer sub;
  typedef itk::AddImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > AddType;
  typename AddType::Pointer add;
  if( this->GetRegularizationFilter() != NULL )
    {
    // yes, we do, so we have to compute it
    sub = SubtractType::New();
    sub->SetInput( 0, input );
    sub->SetInput( 1, ifft->GetOutput() );
    sub->SetNumberOfThreads( this->GetNumberOfThreads() );
    sub->SetReleaseDataFlag( true );
    // don't run in place - we need to keep the input image
    // sub->SetInPlace( true );
    // connect the regularization filter
    this->GetRegularizationFilter()->SetInput( sub->GetOutput() );
    // and recompute I at iteration n
    add = AddType::New();
    add->SetInput( 0, input );
    add->SetInput( 1, this->GetRegularizationFilter()->GetOutput() );
    add->SetNumberOfThreads( this->GetNumberOfThreads() );
    add->SetReleaseDataFlag( true );
    add->SetInPlace( true );
    }
  
  // divide the input by (the convolved image + epsilon)
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::RichardsonLucy< TInternalPrecision > >
                  RichardsonLucyType;
  typename RichardsonLucyType::Pointer ediv = RichardsonLucyType::New();
  ediv->SetInput( 1, input );
  if( this->GetRegularizationFilter() == NULL )
    {
    ediv->SetInput( 0, ifft->GetOutput() );
    }
  else
    {
    ediv->SetInput( 0, add->GetOutput() );
    }
  ediv->SetNumberOfThreads( this->GetNumberOfThreads() );
  ediv->SetReleaseDataFlag( true );
  ediv->SetInPlace( true );
  
  // convolve the divided image by the psf
  
  typename FFTFilterType::Pointer fft2 = FFTFilterType::New();
  fft2->SetInput( ediv->GetOutput() );
  fft2->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft2->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( fft, 0.25f );

  typedef itk::BinaryFunctorImageFilter< ComplexImageType,
                ComplexImageType,
                ComplexImageType,
                typename Functor::MultiplyByComplexConjugate< ComplexType > >
                  MultiplyByComplexConjugateType;
  typename MultiplyByComplexConjugateType::Pointer cmult = MultiplyByComplexConjugateType::New();
  cmult->SetInput( 0, fft2->GetOutput() );
  cmult->SetInput( 1, psf );
  cmult->SetNumberOfThreads( this->GetNumberOfThreads() );
  cmult->SetReleaseDataFlag( true );
  cmult->SetInPlace( true );
  
  typename IFFTFilterType::Pointer ifft2 = IFFTFilterType::New();
  ifft2->SetInput( cmult->GetOutput() );
  ifft2->SetActualXDimensionIsOdd( xIsOdd );
  ifft2->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft2->SetReleaseDataFlag( true );
  
  // divided image convolution completed
  // multiply the result with the input
  
  typedef itk::MultiplyImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType > RealMultType;
  typename RealMultType::Pointer rmult = RealMultType::New();
  rmult->SetInput( 0, ifft2->GetOutput() );
  rmult->SetInput( 1, input );
  rmult->SetNumberOfThreads( this->GetNumberOfThreads() );
  // can't be released, it is required by two filters on next iteration
  // rmult->SetReleaseDataFlag( true );
  rmult->SetInPlace( true );
  
  // begin the iterations
  typedef typename itk::RelativeChangeCalculator< InternalImageType > ChangeType;
  typename ChangeType::Pointer change = ChangeType::New();
  typename InternalImageType::Pointer img = input;
  for( int i=1; i<=this->GetNumberOfIterations(); i++ )
    {
    this->SetIteration( i );
    // should we use smoothing filter? -- tested in the iteration on purpose, to be able to
    // change the filter by looking at the iteration event
    typename InternalFilterType::Pointer last = rmult.GetPointer();
    if( this->GetSmoothingFilter() != NULL && i % this->GetSmoothingPeriod() == 0 )
      {
      this->GetSmoothingFilter()->SetInput( rmult->GetOutput() );
      last = this->GetSmoothingFilter();
      }
    last->Update();
    
    // do we have to stop the iterations based on the relative change?
    change->SetImage( img );
    change->SetNewImage( last->GetOutput() );
    change->Compute();
    // std::cout << change->GetOutput() << std::endl;
    this->SetRelativeChange( change->GetOutput() );
    if( this->GetRelativeChangeThreshold() > 0 && this->GetRelativeChange() < this->GetRelativeChangeThreshold() )
      {
      break;
      }
    else
      {
      // ok, lets go for another round
      img = last->GetOutput();
      img->DisconnectPipeline();
      fft->SetInput( img );
      rmult->SetInput( 1, img );
      this->UpdateProgress( i/(float)this->GetNumberOfIterations() );
      this->InvokeEvent( IterationEvent() );
      }
    }
  img->SetReleaseDataFlag( true );
  this->UpdateProgress( 1.0 );
   
  this->End( img, 0 );
}

}// end namespace itk
#endif
