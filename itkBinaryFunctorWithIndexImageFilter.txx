/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryFunctorWithIndexImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-07 17:31:02 $
  Version:   $Revision: 1.40 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryFunctorWithIndexImageFilter_txx
#define __itkBinaryFunctorWithIndexImageFilter_txx

#include "itkBinaryFunctorWithIndexImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage1, class TInputImage2, 
          class TOutputImage, class TFunction  >
BinaryFunctorWithIndexImageFilter<TInputImage1,TInputImage2,TOutputImage,TFunction>
::BinaryFunctorWithIndexImageFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->InPlaceOff();
}


/**
 * Connect one of the operands for pixel-wise addition
 */
template <class TInputImage1, class TInputImage2, 
          class TOutputImage, class TFunction  >
void
BinaryFunctorWithIndexImageFilter<TInputImage1,TInputImage2,TOutputImage,TFunction>
::SetInput1( const TInputImage1 * image1 ) 
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage1 *>( image1 ));
}


/**
 * Connect one of the operands for pixel-wise addition
 */
template <class TInputImage1, class TInputImage2, 
          class TOutputImage, class TFunction  >
void
BinaryFunctorWithIndexImageFilter<TInputImage1,TInputImage2,TOutputImage,TFunction>
::SetInput2( const TInputImage2 * image2 ) 
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TInputImage2 *>( image2 ));
}

/**
 * ThreadedGenerateData Performs the pixel-wise addition
 */
template <class TInputImage1, class TInputImage2, class TOutputImage, class TFunction  >
void
BinaryFunctorWithIndexImageFilter<TInputImage1, TInputImage2, TOutputImage, TFunction>
::ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread,
                        int threadId)
{
  // We use dynamic_cast since inputs are stored as DataObjects.  The
  // ImageToImageFilter::GetInput(int) always returns a pointer to a
  // TInputImage1 so it cannot be used for the second input.
  Input1ImagePointer inputPtr1
    = dynamic_cast<const TInputImage1*>(ProcessObject::GetInput(0));
  Input2ImagePointer inputPtr2
    = dynamic_cast<const TInputImage2*>(ProcessObject::GetInput(1));
  OutputImagePointer outputPtr = this->GetOutput(0);
  
  ImageRegionConstIterator<TInputImage1> inputIt1(inputPtr1, outputRegionForThread);
  ImageRegionConstIterator<TInputImage2> inputIt2(inputPtr2, outputRegionForThread);

  ImageRegionIteratorWithIndex<TOutputImage> outputIt(outputPtr, outputRegionForThread);

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  inputIt1.GoToBegin();
  inputIt2.GoToBegin();
  outputIt.GoToBegin();

  while( !inputIt1.IsAtEnd() ) 
    {
    outputIt.Set( m_Functor( outputIt.GetIndex(), inputIt1.Get(), inputIt2.Get() ) );
    ++inputIt2;
    ++inputIt1;
    ++outputIt;
    progress.CompletedPixel(); // potential exception thrown here
    }
}

} // end namespace itk

#endif
