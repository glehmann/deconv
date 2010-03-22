/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkIterativeDeconvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkIterativeDeconvolutionImageFilter_h
#define __itkIterativeDeconvolutionImageFilter_h

#include "itkFFTConvolutionImageFilterBase.h"
#include "itkConceptChecking.h"

namespace itk {

/** \class IterativeDeconvolutionImageFilter
 * \brief 
 *
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TPointSpreadFunction=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT IterativeDeconvolutionImageFilter : 
    public FFTConvolutionImageFilterBase<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
{
public:
  /** Standard class typedefs. */
  typedef IterativeDeconvolutionImageFilter Self;

  typedef FFTConvolutionImageFilterBase<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision> Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TPointSpreadFunction                     PointSpreadFunctionType;
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
  
  typedef typename Superclass::ComplexImageType    ComplexImageType;
  typedef typename ComplexImageType::Pointer       ComplexImagePointerType;
  typedef typename ComplexImageType::PixelType     ComplexType;

  typedef typename Superclass::InternalImageType                                   InternalImageType;
  typedef typename itk::ImageToImageFilter< InternalImageType, InternalImageType > InternalFilterType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(IterativeDeconvolutionImageFilter, FFTConvolutionImageFilterBase);

  /**
   * Set/Get the number of iterations to run. If RelativeChangeThreshold is set to a value
   * greater than zero, then this parameter is the maximum number of iterations which can
   * be exceeded even if the RelativeChangeThreshold has not been reached.
   * Defaults to 10.
   */
  itkGetConstMacro(NumberOfIterations, int);
  itkSetMacro(NumberOfIterations, int);
  
  /**
   * Set/Get the relative change threshold between two iterations to stop the iteration.
   * A value lower or equal to 0 mean that this feature is not used. A usual value is 
   * between 10^-3 and 10^-5.
   * Defaults to 0 (not used).
   */
  itkGetConstMacro(RelativeChangeThreshold, double);
  itkSetMacro(RelativeChangeThreshold, double);

  /**
   * Set/Get the smoothing which filter is applied at the end of an iteration, each
   * SmoothingPeriod iterations. Default is NULL (no smoothing).
   */
  itkGetConstObjectMacro(SmoothingFilter, InternalFilterType);
  itkGetObjectMacro(SmoothingFilter, InternalFilterType);
  itkSetObjectMacro(SmoothingFilter, InternalFilterType);

  /**
   * Set/Get how often a smoothing is applied with the SmoothingFilter.
   * By default, the smoothing is applied every iterations.
   */
  itkGetConstMacro(SmoothingPeriod, int);
  itkSetMacro(SmoothingPeriod, int);

  /**
   * Set/Get the regularization filter which is applied during each iteration on
   * the residual image. This filter should keep only the noise in the image.
   * Default is NULL (no regularization).
   */
  itkGetConstObjectMacro(RegularizationFilter, InternalFilterType);
  itkGetObjectMacro(RegularizationFilter, InternalFilterType);
  itkSetObjectMacro(RegularizationFilter, InternalFilterType);

  itkGetConstMacro(Iteration, int);
  itkGetConstMacro(RelativeChange, double);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  IterativeDeconvolutionImageFilter();
  ~IterativeDeconvolutionImageFilter() {};

  void PrintSelf(std::ostream& os, Indent indent) const;

  itkSetMacro(Iteration, int);
  itkSetMacro(RelativeChange, double);

private:
  IterativeDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int                                        m_NumberOfIterations;
  double                                     m_RelativeChangeThreshold;
  typename InternalFilterType::Pointer       m_SmoothingFilter;
  int                                        m_Iteration;
  double                                     m_RelativeChange;
  int                                        m_SmoothingPeriod;
  typename InternalFilterType::Pointer       m_RegularizationFilter;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIterativeDeconvolutionImageFilter.txx"
#endif

#endif
