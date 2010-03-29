/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConchelloIntensityPenalizationImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 16:45:11 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkConchelloIntensityPenalizationImageFilter_txx
#define __itkConchelloIntensityPenalizationImageFilter_txx

#include "itkConchelloIntensityPenalizationImageFilter.h"

namespace itk
{

/**
 *
 */
template <class TInputImage, class TOutputImage>
ConchelloIntensityPenalizationImageFilter<TInputImage, TOutputImage>
::ConchelloIntensityPenalizationImageFilter()
{
  m_Lambda = 1.0/255;
}

template <class TInputImage, class TOutputImage>
void
ConchelloIntensityPenalizationImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  this->GetFunctor().SetLambda(m_Lambda);
}

/**
 *
 */
template <class TInputImage, class TOutputImage>
void 
ConchelloIntensityPenalizationImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Lambda: "
     << static_cast<typename NumericTraits<InputPixelType>::PrintType>(m_Lambda)
     << std::endl;
}


} // end namespace itk

#endif
