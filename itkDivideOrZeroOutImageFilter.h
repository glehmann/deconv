/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDivideOrZeroOutImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-01 14:36:25 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDivideOrZeroOutImageFilter_h
#define __itkDivideOrZeroOutImageFilter_h

#include "itkBinaryFunctorImageFilter.h"

namespace itk
{

namespace Functor {  
 
template< class TInput1, class TInput2, class TOutput >
class DivideOrZeroOut
{
public:
  DivideOrZeroOut()
    {
    m_Threshold = 1e-5;
    };
  ~DivideOrZeroOut() {};
  bool operator!=( const DivideOrZeroOut & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const DivideOrZeroOut & other ) const
    {
    return true;
    }
  // arguments are in an unusual order to allow to run in place
  // of the denominator image
  inline TOutput operator()( const TInput1 & d, const TInput2 & n )
    {
    if( d < m_Threshold )
      {
      return NumericTraits< TOutput >::Zero;
      }
    return static_cast< TOutput >( n / d );
    }
  double m_Threshold;
};

}  // end namespace functor


/** \class DivideOrZeroOutImageFilter
 * \brief 
 *
 *
 * \author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * \sa IntensityWindowingImageFilter ShiftScaleImageFilter
 * \ingroup IntensityImageFilters  Multithreaded
 *
 */
template <typename  TInputImage1, typename  TInputImage2=TInputImage1, typename  TOutputImage=TInputImage1>
class ITK_EXPORT DivideOrZeroOutImageFilter :
    public
    BinaryFunctorImageFilter<TInputImage1, TInputImage2,TOutputImage, 
                            Functor::DivideOrZeroOut< 
                                typename TInputImage1::PixelType, 
                                typename TInputImage2::PixelType, 
                                typename TOutputImage::PixelType > >
{
public:
  /** Standard class typedefs. */
  typedef DivideOrZeroOutImageFilter    Self;
  typedef BinaryFunctorImageFilter<TInputImage1, TInputImage2, TOutputImage, 
                            Functor::DivideOrZeroOut< 
                                typename TInputImage1::PixelType, 
                                typename TInputImage2::PixelType, 
                                typename TOutputImage::PixelType > > Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Runtime information support. */
  itkTypeMacro(DivideOrZeroOutImageFilter, 
               BinaryFunctorImageFilter);

  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "Threshold: "  << GetThreshold() << std::endl;
    }
  
  void SetThreshold( double v )
    {
    this->GetFunctor().m_Threshold = v;
    this->Modified();
    }
  double GetThreshold() const
    {
    return this->GetFunctor().m_Threshold;
    }
    
protected:
  DivideOrZeroOutImageFilter() {};
  virtual ~DivideOrZeroOutImageFilter() {};

private:
  DivideOrZeroOutImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};


  
} // end namespace itk
#endif
