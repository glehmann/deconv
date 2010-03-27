/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter_h
#define __itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter_h

#include "itkIterativeDeconvolutionImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

namespace Functor {  
  
template< class TReal>
class TikhonovMillerRichardsonLucy
{
public:
  TikhonovMillerRichardsonLucy() {};
  ~TikhonovMillerRichardsonLucy() {};
  bool operator!=( const TikhonovMillerRichardsonLucy & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const TikhonovMillerRichardsonLucy & other ) const
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
class TikhonovMillerRegularization
{
public:
  TikhonovMillerRegularization()
   {
   m_Lambda = NumericTraits< TReal >::min();
   };
  ~TikhonovMillerRegularization() {};
  bool operator!=( const TikhonovMillerRegularization & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const TikhonovMillerRegularization & other ) const
    {
    return true;
    }
  inline TReal operator()( const TReal & e, const TReal & x, const TReal & l )
    {
    return std::max( e * x / ( 1 - m_Lambda * x ), NumericTraits< TReal >::Zero );
    }
  TReal m_Lambda;
};

}

/** \class TikhonovMillerRichardsonLucyDeconvolutionImageFilter
 * \brief 
 *
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TPointSpreadFunction=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT TikhonovMillerRichardsonLucyDeconvolutionImageFilter : 
    public IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision> 
{
public:
  /** Standard class typedefs. */
  typedef TikhonovMillerRichardsonLucyDeconvolutionImageFilter Self;

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
  itkTypeMacro(TikhonovMillerRichardsonLucyDeconvolutionImageFilter, IterativeDeconvolutionImageFilter);

  itkGetConstMacro(Lambda, InternalPrecisionType);
  itkSetMacro(Lambda, InternalPrecisionType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  TikhonovMillerRichardsonLucyDeconvolutionImageFilter();
  ~TikhonovMillerRichardsonLucyDeconvolutionImageFilter() {};

  /** Single-threaded version of GenerateData.  This filter delegates
   * to other filters. */
  void GenerateData();

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  TikhonovMillerRichardsonLucyDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  InternalPrecisionType m_Lambda;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTikhonovMillerRichardsonLucyDeconvolutionImageFilter.txx"
#endif

#endif
