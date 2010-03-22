/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkJanssonVanCittertDeconvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkJanssonVanCittertDeconvolutionImageFilter_h
#define __itkJanssonVanCittertDeconvolutionImageFilter_h

#include "itkVanCittertDeconvolutionImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

namespace Functor {  
  
template< class TReal>
class JanssonVanCittert
{
public:
  JanssonVanCittert() {};
  ~JanssonVanCittert() {};
  bool operator!=( const JanssonVanCittert & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const JanssonVanCittert & other ) const
    {
    return true;
    }
  inline TReal operator()( const TReal & residual, const TReal & img )
    {
    TReal gamma = ( img - m_A ) / m_A;
    TReal result = residual * ( 1 - gamma*gamma ) * m_Alpha + img;
    if( m_NonNegativity && result < 0 )
      {
      result = 0;
      }
    return result;
    }
  TReal m_A;
  TReal m_Alpha;
  bool  m_NonNegativity;
};

}

/** \class JanssonVanCittertDeconvolutionImageFilter
 * \brief 
 *
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TPointSpreadFunction=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT JanssonVanCittertDeconvolutionImageFilter : 
    public VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision, Functor::JanssonVanCittert<TInternalPrecision> >
{
public:
  /** Standard class typedefs. */
  typedef JanssonVanCittertDeconvolutionImageFilter Self;

  typedef VanCittertDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision, Functor::JanssonVanCittert<TInternalPrecision> > Superclass;

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
  typedef typename Superclass::FunctorType         FunctorType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef typename itk::Image< InternalPrecisionType, ImageDimension > InternalImageType;
  typedef typename itk::ImageToImageFilter< InternalImageType, InternalImageType > InternalFilterType;
  
  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(JanssonVanCittertDeconvolutionImageFilter, ImageToImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  JanssonVanCittertDeconvolutionImageFilter() {};
  ~JanssonVanCittertDeconvolutionImageFilter() {};

  virtual void InitFunctor( FunctorType & functor );

private:
  JanssonVanCittertDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkJanssonVanCittertDeconvolutionImageFilter.txx"
#endif

#endif
