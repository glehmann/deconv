/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTotalIntensityRatioCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2005/08/10 16:32:57 $
  Version:   $Revision: 1.52 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTotalIntensityRatioCalculator_txx
#define __itkTotalIntensityRatioCalculator_txx
#include "itkTotalIntensityRatioCalculator.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{ 


template < class TInputImage>
TotalIntensityRatioCalculator<TInputImage>
::TotalIntensityRatioCalculator(void) 
{
  m_Valid = false;
  m_Image = NULL;
  m_NewImage = NULL;
  m_Output = NumericTraits< InputPixelType >::Zero;
}


template < class TInputImage >
void
TotalIntensityRatioCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
  os << indent << "NewImage: " << m_NewImage.GetPointer() << std::endl;
  os << indent << "Valid: " << m_Valid << std::endl;
  os << indent << "Output: " << m_Output << std::endl;
}


template < class TInputImage >
void
TotalIntensityRatioCalculator<TInputImage>
::Compute()
{

  typedef typename InputImageType::IndexType IndexType;

  if( !m_Image || !m_NewImage ) 
    {
    return;
    }

  ImageRegionConstIteratorWithIndex< InputImageType > iIt( m_Image,
                                                     m_Image->GetRequestedRegion() ); 
  iIt.GoToBegin();
  ImageRegionConstIteratorWithIndex< InputImageType > nIt( m_NewImage,
                                                     m_NewImage->GetRequestedRegion() ); 
  nIt.GoToBegin();

  // init the values
  double isum = 0;
  double nsum = 0;
  
  while( !iIt.IsAtEnd() )
    {
    isum += iIt.Get();
    nsum += nIt.Get();
    ++iIt;
    ++nIt;
    }
  
  m_Output = nsum / isum;
  m_Valid = true;

}


template < class TInputImage >
const double &
TotalIntensityRatioCalculator<TInputImage>
::GetOutput() const
{
  if (!m_Valid)
    {
    itkExceptionMacro( << "GetOutput() invoked, but the output have not been computed. Call Compute() first.");
    }
  return m_Output;
}

} // end namespace itk

#endif
