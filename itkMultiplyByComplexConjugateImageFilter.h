/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiplyByComplexConjugateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-01 14:36:31 $
  Version:   $Revision: 1.27 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiplyByComplexConjugateImageFilter_h
#define __itkMultiplyByComplexConjugateImageFilter_h

#include "itkBinaryFunctorImageFilter.h"

namespace itk
{
  
/** \class MultiplyByComplexConjugateImageFilter
 * \brief Implements an operator for pixel-wise multiplication of the first image by the complex 
 * conjugate of the second.
 *
 * Pixel types must be complex.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  

template< class TComplex>
class MultiplyByComplexConjugate
{
public:
  MultiplyByComplexConjugate() {};
  ~MultiplyByComplexConjugate() {};
  bool operator!=( const MultiplyByComplexConjugate & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const MultiplyByComplexConjugate & other ) const
    {
    return true;
    }
  inline TComplex operator()( const TComplex & c1, const TComplex & c2 )
    {
    return c1 * std::conj(c2);
    }
};
}

template <class TComplexImage>
class ITK_EXPORT MultiplyByComplexConjugateImageFilter :
    public
BinaryFunctorImageFilter<TComplexImage,TComplexImage,TComplexImage, 
                         Functor::MultiplyByComplexConjugate< typename TComplexImage::PixelType > >
{
public:
  /** Standard class typedefs. */
  typedef MultiplyByComplexConjugateImageFilter       Self;
  typedef BinaryFunctorImageFilter<TComplexImage,TComplexImage,TComplexImage, 
                         Functor::MultiplyByComplexConjugate< typename TComplexImage::PixelType> >
                                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(MultiplyByComplexConjugateImageFilter, 
               BinaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
#endif

protected:
  MultiplyByComplexConjugateImageFilter() {}
  virtual ~MultiplyByComplexConjugateImageFilter() {}

private:
  MultiplyByComplexConjugateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
