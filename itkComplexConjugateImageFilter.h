/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComplexConjugateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-01 14:36:10 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComplexConjugateImageFilter_h
#define __itkComplexConjugateImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk
{
  
/** \class ComplexConjugateImageFilter
 * \brief Computes pixel-wise the complex conjugate.
 * 
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Function {  
  
template< class TComplex>
class ComplexConjugate
{
public:
  ComplexConjugate() {}
  ~ComplexConjugate() {}
  bool operator!=( const ComplexConjugate & ) const
    {
    return false;
    }
  bool operator==( const ComplexConjugate & other ) const
    {
    return !(*this != other);
    }
  inline TComplex operator()( const TComplex & A ) const
    {
    return std::conj(A);
    }
}; 
}

template <class TComplexImage>
class ITK_EXPORT ComplexConjugateImageFilter :
    public
UnaryFunctorImageFilter<TComplexImage,TComplexImage, 
                        Function::ComplexConjugate<typename TComplexImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef ComplexConjugateImageFilter Self;
  typedef UnaryFunctorImageFilter<
    TComplexImage,TComplexImage, 
    Function::ComplexConjugate< typename TComplexImage::PixelType > >
                                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ComplexConjugateImageFilter, 
               UnaryFunctorImageFilter);

  typedef typename TComplexImage::PixelType                     InputPixelType;
  typedef typename TComplexImage::PixelType                     OutputPixelType;
  typedef typename NumericTraits< InputPixelType >::ValueType   InputPixelValueType;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<InputPixelValueType, OutputPixelType>));
  /** End concept checking */
#endif


protected:
  ComplexConjugateImageFilter() {}
  virtual ~ComplexConjugateImageFilter() {}

private:
  ComplexConjugateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
