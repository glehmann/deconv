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
  m_RelativeChangeCalculator = NULL;
  m_Estimate = NULL;
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::GenerateData()
{
  Init();
  Iterate();
  End();
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Init()
{
  m_RelativeChangeCalculator = ChangeType::New();
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::Iterate()
{
  for( m_Iteration=1; m_Iteration<=this->GetNumberOfIterations(); m_Iteration++ )
    {
    BeforeIteration();
    // should we use smoothing filter? -- tested in the iteration on purpose, to be able to
    // m_RelativeChangeCalculator the filter by looking at the iteration event
    typename InternalImageType::Pointer newEstimate = NewEstimate();
    if( m_SmoothingFilter.IsNotNull() && m_Iteration % m_SmoothingPeriod == 0 )
      {
      m_SmoothingFilter->SetInput( newEstimate );
      newEstimate = m_SmoothingFilter->GetOutput();
      }
    newEstimate->Update();
    newEstimate->DisconnectPipeline();
    
    // do we have to stop the iterations based on the relative m_RelativeChangeCalculator?
    m_RelativeChangeCalculator->SetNewImage( newEstimate );
    m_RelativeChangeCalculator->Compute();
    // std::cout << m_RelativeChangeCalculator->GetOutput() << std::endl;
    m_RelativeChange = m_RelativeChangeCalculator->GetOutput();
    if( m_RelativeChangeThreshold > 0 && m_RelativeChange < m_RelativeChangeThreshold )
      {
      break;
      }
    else
      {
      // ok, lets go for another round
      SetEstimate( newEstimate );
      AfterIteration();
      this->UpdateProgress( m_Iteration/(float)this->GetNumberOfIterations() );
      this->InvokeEvent( IterationEvent() );
      }
    }
  // to keep the actual number of iteration in m_Iteration
  m_Iteration--;
  // nothing more to do - just make sure that the data will be released
  m_Estimate->SetReleaseDataFlag( true );
}

template<class TInputImage, class TPointSpreadFunction, class TOutputImage, class TInternalPrecision>
void
IterativeDeconvolutionImageFilter<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
::End()
{
  this->Superclass::End( m_Estimate, 0 );
  this->UpdateProgress( 1.0 );
  // destroy now useless objects
  m_RelativeChangeCalculator = NULL;
  m_Estimate = NULL;
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
}

}// end namespace itk
#endif
