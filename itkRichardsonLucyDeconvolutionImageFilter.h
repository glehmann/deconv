/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRichardsonLucyDeconvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRichardsonLucyDeconvolutionImageFilter_h
#define __itkRichardsonLucyDeconvolutionImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

namespace Functor {  
  
template< class TReal>
class EpsilonDivide
{
public:
  EpsilonDivide() {};
  ~EpsilonDivide() {};
  bool operator!=( const EpsilonDivide & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const EpsilonDivide & other ) const
    {
    return true;
    }
  // arguments are in an unusual order to allow to run in place
  // of the denominator image
  inline TReal operator()( const TReal & d, const TReal & n )
    {
    return n / ( d + NumericTraits< TReal >::min() + 1 );
    }
};

template< class TComplex>
class MultiplyByComplexConjugate
{
public:
  typedef typename TComplex::value_type RealType;
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

/** \class RichardsonLucyDeconvolutionImageFilter
 * \brief 
 *
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TPointSpreadFunction=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT RichardsonLucyDeconvolutionImageFilter : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef RichardsonLucyDeconvolutionImageFilter Self;

  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TPointSpreadFunction                             PointSpreadFunctionType;
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
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef typename itk::Image< InternalPrecisionType, ImageDimension > InternalImageType;
  typedef typename itk::ImageToImageFilter< InternalImageType, InternalImageType > SmoothingFilterType;
  
  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(RichardsonLucyDeconvolutionImageFilter, ImageToImageFilter);

   /** Set the kernel image */
  void SetPointSpreadFunction(const TPointSpreadFunction *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<TPointSpreadFunction *>(input) );
    }

  /** Get the kernel image */
  const PointSpreadFunctionType * GetPointSpreadFunction() const
    {
    return static_cast<PointSpreadFunctionType*>(
      const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

  /** Set the input image */
  void SetInput1(const TInputImage *input)
    {
    this->SetInput( input );
    }

  /** Set the kernel image */
  void SetInput2(const TPointSpreadFunction *input)
    {
    this->SetPointSpreadFunction( input );
    }

  /**
   * Set/Get the greatest prime factor allowed on the size of the padded image.
   * The filter increase the size of the image to reach a size with the greatest
   * prime factor smaller or equal to the specified value. The default value is
   * 13, which is the greatest prime number for which the FFT are precomputed
   * in FFTW, and thus gives very good performance.
   * A greatest prime factor of 2 produce a size which is a power of 2, and thus
   * is suitable for vnl base fft filters.
   * A greatest prime factor of 1 or less - typically 0 - disable the extra padding.
   *
   * Warning: this parameter is not used (and useful) only when ITK is built with
   * FFTW support.
   */
  itkGetConstMacro(GreatestPrimeFactor, int);
  itkSetMacro(GreatestPrimeFactor, int);
  
  /**
   * Set/Get the padding method.
   */
  typedef enum { NO_PADDING=0, ZERO_FLUX_NEUMANN=1, ZERO=2, MIRROR=3, WRAP=4 } PadMethod;
  itkGetConstMacro(PadMethod, int);
  itkSetMacro(PadMethod, int);
  
  /**
   * Set/Get whether the kernel should be normalized to one or not.
   * Default is true.
   */
  itkGetConstMacro(Normalize, bool);
  itkSetMacro(Normalize, bool);
  itkBooleanMacro(Normalize);
  
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
   * Set/Get the smoothing which is applied at the end of an iteration, each
   * SmoothingPeriod iterations. Default is NULL (no smoothing).
   */
  itkGetConstObjectMacro(SmoothingFilter, SmoothingFilterType);
  itkGetObjectMacro(SmoothingFilter, SmoothingFilterType);
  itkSetObjectMacro(SmoothingFilter, SmoothingFilterType);

  /**
   * Set/Get how often a smoothing is applied with the SmoothingFilter.
   * By default, the smoothing is applied every iterations.
   */
  itkGetConstMacro(SmoothingPeriod, int);
  itkSetMacro(SmoothingPeriod, int);

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
  RichardsonLucyDeconvolutionImageFilter();
  ~RichardsonLucyDeconvolutionImageFilter() {};

  void GenerateInputRequestedRegion();
  
  /** Single-threaded version of GenerateData.  This filter delegates
   * to other filters. */
  void GenerateData();
  
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  RichardsonLucyDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int                                        m_GreatestPrimeFactor;
  int                                        m_PadMethod;
  bool                                       m_Normalize;
  int                                        m_NumberOfIterations;
  double                                     m_RelativeChangeThreshold;
  typename SmoothingFilterType::Pointer      m_SmoothingFilter;
  int                                        m_Iteration;
  double                                     m_RelativeChange;
  int                                        m_SmoothingPeriod;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRichardsonLucyDeconvolutionImageFilter.txx"
#endif

#endif
