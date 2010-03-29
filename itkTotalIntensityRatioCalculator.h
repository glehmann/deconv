/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTotalIntensityRatioCalculator.h,v $
  Language:  C++
  Date:      $Date: 2004/04/25 23:59:26 $
  Version:   $Revision: 1.37 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTotalIntensityRatioCalculator_h
#define __itkTotalIntensityRatioCalculator_h

#include "itkMacro.h"
#include "itkImage.h"

namespace itk
{

/** \class TotalIntensityRatioCalculator
 * \brief Compute the PSNR of a noisy image
 *
 * \ingroup Operators
 */
template < class TInputImage >
class ITK_EXPORT TotalIntensityRatioCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef TotalIntensityRatioCalculator Self;
  typedef Object                             Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef SmartPointer<const Self>           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TotalIntensityRatioCalculator, Object);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard image type within this class. */
  typedef TInputImage    InputImageType;

  /** Standard image type pointer within this class. */
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::PixelType       InputPixelType;

  /** Set the input image. */
  virtual void SetImage( const InputImageType * image )
    {
    if ( m_Image != image )
      {
      m_Image = image;
      this->Modified();
      m_Valid = false;
      }
    }

  virtual void SetNewImage( const InputImageType * image )
    {
    if ( m_NewImage != image )
      {
      m_NewImage = image;
      this->Modified();
      m_Valid = false;
      }
    }

  void Compute( void );
  
  const double & GetOutput() const;

protected:
  TotalIntensityRatioCalculator();
  virtual ~TotalIntensityRatioCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  TotalIntensityRatioCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool           m_Valid;                      // Have moments been computed yet?
  double         m_Output;

  InputImageConstPointer    m_Image;
  InputImageConstPointer    m_NewImage;

};  // class TotalIntensityRatioCalculator

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTotalIntensityRatioCalculator.txx"
#endif

#endif /* __itkTotalIntensityRatioCalculator_h */
