/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkIterativeDeconvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkIterativeDeconvolutionImageFilter_txx
#define __itkIterativeDeconvolutionImageFilter_txx

#include "itkIterativeDeconvolutionImageFilter.h"

namespace itk {

template <class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::IterativeDeconvolutionImageFilter()
{
  m_NumberOfIterations = 10;
  m_RelativeChangeThreshold = 0;
  m_SmoothingFilter = NULL;
  m_Iteration = 0;
  m_RelativeChange = 0;
  m_SmoothingPeriod = 1;
  m_RegularizationFilter = NULL;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Iteration: "  << m_Iteration << std::endl;
  os << indent << "NumberOfIterations: "  << m_NumberOfIterations << std::endl;
  os << indent << "RelativeChange: "  << m_RelativeChange << std::endl;
  os << indent << "RelativeChangeThreshold: "  << m_RelativeChangeThreshold << std::endl;
  os << indent << "SmoothingPeriod: "  << m_SmoothingPeriod << std::endl;
//  os << indent << "SmoothingFilter: "  << m_SmoothingFilter << std::endl;
  os << indent << "SmoothingFilter: ";
  if( m_SmoothingFilter.IsNull() )
    {
    std::cout << "NULL" << std::endl;
    }
  else
    {
    m_SmoothingFilter->Print( os, indent.GetNextIndent() );
    }
  os << indent << "RegularizationFilter: ";
  if( m_RegularizationFilter.IsNull() )
    {
    std::cout << "NULL" << std::endl;
    }
  else
    {
    m_RegularizationFilter->Print( os, indent.GetNextIndent() );
    }
}

}// end namespace itk
#endif
