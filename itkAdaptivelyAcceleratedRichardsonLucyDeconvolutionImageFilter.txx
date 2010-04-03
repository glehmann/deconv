/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter_txx
#define __itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter_txx

#include "itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyByComplexConjugateImageFilter.h"
#include "itkRelativeChangeCalculator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkForwardDifferenceOperator.h"
#include "itkDivideOrZeroOutImageFilter.h"

namespace itk {

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  double m_V1 = 1;
  double m_V2 = 1;
  double m_Vk_1 = 1;
  double m_Vk = 1;

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
  
  // divide the input by (the convolved image + epsilon)
  typedef itk::DivideOrZeroOutImageFilter< InternalImageType > DivideType;
  typename DivideType::Pointer ediv = DivideType::New();
  ediv->SetInput( 1, input );
  ediv->SetInput( 0, ifft->GetOutput() );
  ediv->SetNumberOfThreads( this->GetNumberOfThreads() );
  ediv->SetReleaseDataFlag( true );
  ediv->SetInPlace( true );
  
  // convolve the divided image by the psf
  
  typename FFTFilterType::Pointer fft2 = FFTFilterType::New();
  fft2->SetInput( ediv->GetOutput() );
  fft2->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft2->SetReleaseDataFlag( true );
  // progress->RegisterInternalFilter( fft, 0.25f );

  typedef itk::MultiplyByComplexConjugateImageFilter< ComplexImageType > MultiplyByComplexConjugateType;
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
  
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::AdaptivelyAcceleratedRichardsonLucy2< TInternalPrecision > > RealMultType;
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
    
    // update the acceleration parameter
    if( this->GetIteration() < 3 )
      {
      rmult->GetFunctor().m_Q = 2;
      }
    else
      {
      rmult->GetFunctor().m_Q = vcl_exp( m_Vk / m_Vk_1 ) - m_V2 / m_V1;
      }
    // std::cout << "Q: " << rmult->GetFunctor().m_Q << std::endl;

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

      // compute the L2 norm of the derivative and store it so it can be reused later
      double v = ComputeL2NormOfFirstOrderDerivative( img );
      if( i == 1 )
        {
        m_V1 = v;
        }
      else if( i == 2 )
        {
        m_V2 = v;
        }
      m_Vk_1 = m_Vk;
      m_Vk = v;
      // std::cout << "v: " << v << std::endl;
    
      this->UpdateProgress( i/(float)this->GetNumberOfIterations() );
      this->InvokeEvent( IterationEvent() );
      }
    }
  img->SetReleaseDataFlag( true );
  this->UpdateProgress( 1.0 );
   
  this->End( img, 0 );
}


template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
double
AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::ComputeL2NormOfFirstOrderDerivative( const InternalImageType * input )
{
  double res = 0;

  typedef ForwardDifferenceOperator< InternalPrecisionType, ImageDimension > OperatorType;
  OperatorType op;
  
  typename InternalImageType::SpacingType spacing = input->GetSpacing();

//  for( int i=0; i<1; i++ )
  for( int i=0; i<ImageDimension; i++ )
    {
    op.SetDirection( i );
    op.CreateDirectional();

    typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InternalImageType>
      BFC;
  
    typedef typename BFC::FaceListType FaceListType;
  
    NeighborhoodInnerProduct<InternalImageType, InternalPrecisionType> smartInnerProduct;
    BFC faceCalculator;
    FaceListType faceList;
   
    // Break the input into a series of regions.  The first region is free
    // of boundary conditions, the rest with boundary conditions. Note,
    // we pass in the input image and the OUTPUT requested region. We are
    // only concerned with centering the neighborhood operator at the
    // pixels that correspond to output pixels.
    faceList = faceCalculator(input, input->GetLargestPossibleRegion(), op.GetRadius());
  
    typename FaceListType::iterator fit;
    // Process non-boundary region and each of the boundary faces.
    // These are N-d regions which border the edge of the buffer.
    ConstNeighborhoodIterator<InternalImageType> bit;
    for (fit=faceList.begin(); fit != faceList.end(); ++fit)
      {
      bit = ConstNeighborhoodIterator<InternalImageType>(op.GetRadius(), input, *fit);
      ZeroFluxNeumannBoundaryCondition<InternalImageType> nbc;
      bit.OverrideBoundaryCondition( &nbc );
      bit.GoToBegin();
      while ( ! bit.IsAtEnd() )
        {
        double v = smartInnerProduct(bit, op);
        v /= spacing[i];
        res += v * v;
        ++bit;
        }
      }
    } 
  return vcl_sqrt( res );
}

}// end namespace itk
#endif
