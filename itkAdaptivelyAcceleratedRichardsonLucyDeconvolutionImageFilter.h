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

#include "itkIterativeDeconvolutionImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

namespace Functor {  
  
template< class TReal>
class AdaptivelyAcceleratedRichardsonLucy
{
public:
  AdaptivelyAcceleratedRichardsonLucy() {};
  ~AdaptivelyAcceleratedRichardsonLucy() {};
  bool operator!=( const AdaptivelyAcceleratedRichardsonLucy & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const AdaptivelyAcceleratedRichardsonLucy & other ) const
    {
    return true;
    }
  // arguments are in an unusual order to allow to run in place
  // of the denominator image
  inline TReal operator()( const TReal & d, const TReal & n )
    {
    if( d < 1e-5 )
      {
      return 0;
      }
    return n / d;
    }
};

template< class TReal>
class AdaptivelyAcceleratedRichardsonLucy2
{
public:
  AdaptivelyAcceleratedRichardsonLucy2()
    {
    m_Q = 1;
    };
  ~AdaptivelyAcceleratedRichardsonLucy2() {};
  bool operator!=( const AdaptivelyAcceleratedRichardsonLucy2 & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const AdaptivelyAcceleratedRichardsonLucy2 & other ) const
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
    return i * vcl_pow( v, m_Q );
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
    public IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision> 
{
public:
  /** Standard class typedefs. */
  typedef AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter Self;

  typedef IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>  Superclass;

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
  itkTypeMacro(AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter, IterativeDeconvolutionImageFilter);

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

  /** Single-threaded version of GenerateData.  This filter delegates
   * to other filters. */
  void GenerateData();

  double ComputeL2NormOfFirstOrderDerivative( const InternalImageType * img );
  
private:
  AdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAdaptivelyAcceleratedRichardsonLucyDeconvolutionImageFilter.txx"
#endif

#endif
