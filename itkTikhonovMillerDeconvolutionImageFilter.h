/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTikhonovMillerDeconvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTikhonovMillerDeconvolutionImageFilter_h
#define __itkTikhonovMillerDeconvolutionImageFilter_h

#include "itkFFTConvolutionImageFilterBase.h"
#include "itkConceptChecking.h"

namespace itk {

namespace Functor {  
  
template< class TComplex>
class TikhonovMillerDeconvolution
{
public:
  typedef typename TComplex::value_type RealType;
  TikhonovMillerDeconvolution()
   {
   m_Gamma = NumericTraits< RealType >::min();
   };
  ~TikhonovMillerDeconvolution() {};
  bool operator!=( const TikhonovMillerDeconvolution & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const TikhonovMillerDeconvolution & other ) const
    {
    return true;
    }
  inline TComplex operator()( const TComplex & input, const TComplex & psf, const TComplex & reg )
    {
    return input * std::conj(psf) / ( std::norm(psf) + m_Gamma * std::norm(reg) );
    }
  RealType m_Gamma;
};
}

/** \class TikhonovMillerDeconvolutionImageFilter
 * \brief 
 *
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TPointSpreadFunction=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT TikhonovMillerDeconvolutionImageFilter : 
    public FFTConvolutionImageFilterBase<TInputImage, TPointSpreadFunction, TOutputImage, TInternalPrecision>
{
public:
  /** Standard class typedefs. */
  typedef TikhonovMillerDeconvolutionImageFilter Self;

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

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef typename itk::Image< InternalPrecisionType, ImageDimension > InternalImageType;
  typedef typename InternalImageType::Pointer                          InternalImagePointerType;

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(TikhonovMillerDeconvolutionImageFilter, ImageToImageFilter);

  /**
   * Set/Get the smoothing parameter. Usual values are in the range 0.001-0.1.
   * Higher value produce smoother results.
   * Defaults to 0.001.
   */
  itkGetConstMacro(Gamma, InternalPrecisionType);
  itkSetMacro(Gamma, InternalPrecisionType);

  /** Set the kernel image */
  void SetRegularizationImage(const InternalImageType *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<InternalImageType *>(input) );
    }

  /** Get the kernel image */
  const InternalImageType * GetRegularizationImage() const
    {
    return static_cast<InternalImageType*>(
      const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  TikhonovMillerDeconvolutionImageFilter();
  ~TikhonovMillerDeconvolutionImageFilter() {};

  /** Single-threaded version of GenerateData.  This filter delegates
   * to other filters. */
  void GenerateData();
  
  void GenerateInputRequestedRegion();

  void PrintSelf(std::ostream& os, Indent indent) const;
  
  InternalImagePointerType Laplacian();

private:
  TikhonovMillerDeconvolutionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  InternalPrecisionType m_Gamma;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTikhonovMillerDeconvolutionImageFilter.txx"
#endif

#endif
