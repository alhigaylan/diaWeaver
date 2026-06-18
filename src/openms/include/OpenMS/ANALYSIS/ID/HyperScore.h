// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <vector>

namespace OpenMS
{

/**
 *  @brief An implementation of the X!Tandem HyperScore PSM scoring function
 */               

struct OPENMS_DLLAPI HyperScore
{
  typedef std::pair<Size, double> IndexScorePair; 

  /** @brief compute the (ln transformed) X!Tandem HyperScore 
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param[in] fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param[in] fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param[in] exp_spectrum measured spectrum
   * @param[in] theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   */
//  static double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const RichPeakSpectrum& theo_spectrum);

  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum);

  /** @brief compute the (ln transformed) X!Tandem HyperScore 
   *  overload that returns some additional information on the match
   */
  struct PSMDetail
  {
    size_t matched_prefix_ions = 0;  ///< N-terminal ions (a, b, c)
    size_t matched_suffix_ions = 0;  ///< C-terminal ions (x, y, z)
    double mean_error = 0.0;
    /// Index into the experimental spectrum for each matched theoretical ion (parallel to matched_ion_types).
    /// Populated by computeWithDetail to allow callers to reconstruct the exact match set without re-running scoring.
    std::vector<Size> matched_exp_idxs;
    /// Effective ion-type character for each match: 'a'/'b'/'c' (prefix) or 'x'/'y'/'z' (suffix); '\0' for other.
    /// Parallel to matched_exp_idxs.
    std::vector<char> matched_ion_types;
  };

  static double computeWithDetail(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum,
                        PSMDetail& d
                       );

  /* @brief compute the (ln transformed) X!Tandem HyperScore only matching peaks that match in charge
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param[in] fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param[in] fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param[in] exp_spectrum measured spectrum
   * @param[in] exp_charges charges of measured peaks
   * @param[in] theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   * @param[in] theo_charges charges of theoretical peaks
  */
  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        const PeakSpectrum& theo_spectrum,
                        const DataArrays::IntegerDataArray& theo_charges);

  /* @brief compute the (ln transformed) X!Tandem HyperScore only matching peaks that match in charge
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param[in] fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param[in] fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param[in] exp_spectrum measured spectrum
   * @param[in] exp_charges charges of measured peaks
   * @param[in] theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   * @param[in] theo_charges charges of theoretical peaks
   * @param[in] intensity_sum summed intensity for observed bond indices (e.g., b3=123 -> intensity_sum[2]=123)
   * Note: intensity_sum must be zeroed and of size #AA in peptide
  */
  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        const PeakSpectrum& theo_spectrum,
                        const DataArrays::IntegerDataArray& theo_charges,
                        std::vector<double>& intensity_sum);

  private:
    /// helper to compute the log factorial
    static double logfactorial_(const int x, int base = 2);
};

}


