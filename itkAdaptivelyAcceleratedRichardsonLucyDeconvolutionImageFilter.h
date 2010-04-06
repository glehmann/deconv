/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter_h
#define __itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter_h

#include "itkRichardsonLucyDeconvolutionImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

namespace Functor {  

template< class TReal>
class AdaptivelyAcceleratedRichardsonLucy
{
public:
  AdaptivelyAcceleratedRichardsonLucy()
    {
    m_Q = 1;
    };
  ~AdaptivelyAcceleratedRichardsonLucy() {};
  bool operator!=( const AdaptivelyAcceleratedRichardsonLucy & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const AdaptivelyAcceleratedRichardsonLucy & other ) const
    {
    return true;
    }
  inline TReal operator()( const TReal & v, const TReal & i )
    {
    if( v < 1e-5 )
      {
      return 0;
      }
    assert( !std::isnan( i ) );
    assert( !std::isnan( v ) );
    assert( !std::isnan( m_Q ) );
//     if( std::isnan( vcl_pow( v, m_Q ) ) )
//       std::cout << "v: " << v << "  m_Q: " << m_Q << std::endl;
    assert( !std::isnan( vcl_pow( v, m_Q ) ) );
    assert( !std::isnan( i * vcl_pow( v, m_Q ) ) );
    return std::max( i * vcl_pow( v, m_Q ),
                     NumericTraits< TReal >::Zero );
    }
  TReal m_Q;
};

}

/** \class AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter
 * \brief 
 *
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TPointSpreadFunction=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter : 
    public RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision> 
{
public:
  /** Standard class typedefs. */
  typedef AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter Self;

  typedef RichardsonLucyDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>  Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TPointSpreadFunction                             PointSpreadFunctionType;
  typedef TOutputImage                             OutputImageType;
  typedef TInternalPrecision                            InternalPrecisionType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename PointSpreadFunctionType::Pointer        PointSpreadFunctionPointer;
  typedef typename PointSpreadFunctionType::ConstPointer   PointSpreadFunctionConstPointer;
  typedef typename PointSpreadFunctionType::PixelType      PointSpreadFunctionPixelType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename InputImageType::RegionType      RegionType;
  typedef typename InputImageType::IndexType       IndexType;
  typedef typename InputImageType::SizeType        SizeType;
  
  
  typedef typename Superclass::FFTFilterType       FFTFilterType;
  typedef typename Superclass::IFFTFilterType      IFFTFilterType;
  typedef typename Superclass::ComplexImageType    ComplexImageType;
  typedef typename ComplexImageType::Pointer       ComplexImagePointerType;
  typedef typename ComplexImageType::PixelType     ComplexType;

  typedef typename Superclass::InternalImageType   InternalImageType;
  typedef typename InternalImageType::Pointer      InternalImagePointerType;
  typedef typename Superclass::InternalFilterType  InternalFilterType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter, RichardsonLucyDeconvolutionImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter() {};
  ~AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter() {};

  virtual void Init();
  virtual void BeforeIteration();

  double ComputeL2NormOfFirstOrderDerivative( const InternalImageType * img );
  
  typedef itk::BinaryFunctorImageFilter< InternalImageType,
                InternalImageType,
                InternalImageType,
                typename Functor::AdaptivelyAcceleratedRichardsonLucy< TInternalPrecision > > MultType;

private:
  AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  double m_V1;
  double m_V2;
  double m_Vk_1;
  double m_Vk;
  
}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter.txx"
#endif

#endif
