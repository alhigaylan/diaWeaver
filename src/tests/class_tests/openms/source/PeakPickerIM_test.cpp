// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <cmath>

///////////////////////////
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerIM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerIM* ptr = nullptr;
PeakPickerIM* nullPointer = nullptr;

START_SECTION((PeakPickerIM()))
  ptr = new PeakPickerIM();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerIM()))
  delete ptr;
END_SECTION


MSSpectrum input;

// two m/z profile peaks at 300 and 800 m/z, each with 5 ion mobility scans (0.4 to 1.2 Vs/cm2) and Gaussian intensity distribution along IM
// each m/z peak has 4 points per IM scan (0.01 m/z spacing, 0.1 Vs/cm2 spacing)

///std::vector<double> mzs = {299.91, 299.97, 300.03, 300.09, 299.91, 299.97, 300.03, 300.09, 299.91, 299.97, 300.03, 300.09, 299.91, 299.97, 300.03, 300.09, 299.91, 299.97, 300.03, 300.09, 799.88, 799.96, 800.04, 800.12, 799.88, 799.96, 800.04, 800.12, 799.88, 799.96, 800.04, 800.12, 799.88, 799.96, 800.04, 800.12, 799.88, 799.96, 800.04, 800.12};
///std::vector<double> IMs = {0.4469780242779817, 0.4469780242779817, 0.4469780242779817, 0.4469780242779817, 0.5669780242779817, 0.5669780242779817, 0.5669780242779817, 0.5669780242779817, 0.6869780242779817, 0.6869780242779817, 0.6869780242779817, 0.6869780242779817, 0.8069780242779817, 0.8069780242779817, 0.8069780242779817, 0.8069780242779817, 0.9269780242779817, 0.9269780242779817, 0.9269780242779817, 0.9269780242779817, 0.21943921987602621, 0.21943921987602621, 0.21943921987602621, 0.21943921987602621, 0.36943921987602624, 0.36943921987602624, 0.36943921987602624, 0.36943921987602624, 0.5194392198760263, 0.5194392198760263, 0.5194392198760263, 0.5194392198760263, 0.6694392198760263, 0.6694392198760263, 0.6694392198760263, 0.6694392198760263, 0.8194392198760263, 0.8194392198760263, 0.8194392198760263, 0.8194392198760263};
///std::vector<double> Intensities = {0.1234098040869882, 6.737946999091595, 6.737946999091595, 0.1234098040869882, 3.606563136024751, 196.91167520437313, 196.91167520437313, 3.606563136024751, 11.108996538270091, 606.530659713185, 606.530659713185, 11.108996538270091, 3.606563136024751, 196.91167520437313, 196.91167520437313, 3.606563136024751, 0.1234098040869882, 6.737946999091595, 6.737946999091595, 0.1234098040869882, 0.6170490204331873, 33.689734995457975, 33.689734995457975, 0.6170490204331873, 18.032815680072503, 984.5583760218657, 984.5583760218657, 18.032815680072503, 55.54498269119259, 3032.6532985659255, 3032.6532985659255, 55.54498269119259, 18.032815680072503, 984.5583760218657, 984.5583760218657, 18.032815680072503, 0.6170490204331873, 33.689734995457975, 33.689734995457975, 0.6170490204331873};

// Synthetic 2D Gaussian peak around m/z 300 with two IM peaks (0.8, 1.0)
std::vector<double> mzs;
std::vector<double> IMs;
std::vector<double> Intensities;
// --- Parameters for m/z Gaussian ---
const double mz_center = 300.0;
const double mz_sigma  = 0.0025;
const double mz_min    = 299.99;
const double mz_max    = 300.01;
const double mz_step   = 0.0005;
// --- Parameters for IM Gaussians ---
const double im_center1 = 0.8;
const double im_center2 = 1.0;
const double im_sigma   = 0.03;
const double im_min     = 0.7;
const double im_max     = 1.1;
const double im_step    = 0.01;

// Build the 2D grid: for each m/z, duplicate across all IM scans
for (double mz = mz_min; mz <= mz_max + 1e-9; mz += mz_step)
{
  // 1D Gaussian in m/z
  const double mz_int = std::exp(-0.5 * std::pow((mz - mz_center) / mz_sigma, 2.0));
  for (double im = im_min; im <= im_max + 1e-9; im += im_step)
  {
    // Two Gaussians in IM dimension
    const double im_int1 = std::exp(-0.5 * std::pow((im - im_center1) / im_sigma, 2.0));
    const double im_int2 = std::exp(-0.5 * std::pow((im - im_center2) / im_sigma, 2.0));
    const double im_int  = im_int1 + im_int2;
    // Final separable 2D intensity
    const double intensity = mz_int * im_int;
    mzs.push_back(mz);
    IMs.push_back(im);
    Intensities.push_back(intensity);
  }
}

for (size_t i = 0; i < mzs.size(); ++i)
{
    input.emplace_back(mzs[i], Intensities[i]);
}
input.getFloatDataArrays().resize(1);
input.getFloatDataArrays()[0].setName(Constants::UserParam::ION_MOBILITY);
for (size_t i = 0; i < IMs.size(); ++i)
{
    input.getFloatDataArrays()[0].push_back(IMs[i]);
}
input.setIMFormat(IMFormat::IM_PEAK);

/* TODO make it more robust to also handle sparse data (e.g. asymetrical or even only one m/z value per IM scan) 
// create dummy ion mobility spectrum with concatenated scans
{
    const size_t points_per_im = 10;

    // For each ion mobility value (1.0 to 10.0)
    for (double im = 1.0; im <= 10.0; im += 1.0)
    {
        double baseIntensity = (im <= 5.0) ? 
                            200.0 + ((im - 1.0) * 200.0) :    // Increasing intensity from IM 1.0 to 5.0
                            1000.0 - ((im - 5.0) * 200.0);    // Decreasing intensity from IM 6.0 to 10.0
        
        // Create 10 data points for each ion mobility value
        for (size_t i = 0; i < points_per_im; i++)
        {
            double mz = 100.0 + (i * 0.01);
            double intensity = baseIntensity + (i * 20.0);
            input.emplace_back(mz, intensity);
        }
    }

    // Set up the ion mobility array
    input.getFloatDataArrays().resize(1);
    input.getFloatDataArrays()[0].setName(Constants::UserParam::ION_MOBILITY);

    // Add ion mobility values (each repeated 10 times)
    for (double im = 1.0; im <= 10.0; im += 1.0)
    {
        for (size_t i = 0; i < points_per_im; i++)
        {
            input.getFloatDataArrays()[0].push_back(im);
        }
    }
}
*/

START_SECTION(void pickIMTraces(MSSpectrum& spectrum))
{
    PeakPickerIM pp_im;
    
    pp_im.pickIMTraces(input);

    TEST_EQUAL(input.size(), 2)
    TEST_REAL_SIMILAR(input[0].getIntensity(), 5.70646)
    TEST_REAL_SIMILAR(input[0].getMZ(), 300.0)

    TEST_EQUAL(input.getFloatDataArrays().size(), 1)
    TEST_EQUAL(input.getFloatDataArrays()[0].getName(), Constants::UserParam::ION_MOBILITY_CENTROID)
    TEST_REAL_SIMILAR(input.getFloatDataArrays()[0][0], 0.8)
}
END_SECTION

START_SECTION(static double aggregationRtFwhm(std::vector<double> rts, int frames_per_side, double cutoff, double rt_fwhm))
{
    // number of frames on each side of the centre (index centre) whose weight is at least the cutoff
    auto frames_pooled_per_side = [](const std::vector<double>& rts, size_t centre, double fwhm, double cutoff)
    {
        size_t n = 0;
        for (size_t j = centre + 1; j < rts.size(); ++j)
        {
            const double d = rts[j] - rts[centre];
            if (std::exp(-4.0 * std::log(2.0) * d * d / (fwhm * fwhm)) < cutoff) break;
            ++n;
        }
        return n;
    };

    const double cutoff = 0.01;
    const double half_width_per_fwhm = std::sqrt(-std::log(cutoff) / (4.0 * std::log(2.0)));  // 1.2888
    const double spacing = 0.738;
    std::vector<double> grid;
    for (int i = 0; i < 40; ++i) grid.push_back(i * spacing);

    // negative frames_per_side: rt_FWHM is returned as given
    TEST_REAL_SIMILAR(PeakPickerIM::aggregationRtFwhm(grid, -1, cutoff, 3.0), 3.0)

    // ±1 frame at 0.738 s spacing: half-width 1.5 * 0.738 s, divided by sqrt(-ln(0.01) / (4 ln 2)) = 1.2888
    TEST_REAL_SIMILAR(PeakPickerIM::aggregationRtFwhm(grid, 1, cutoff, 3.0), 1.5 * spacing / half_width_per_fwhm)

    // a regular grid pools exactly k frames on each side
    for (int k : {0, 1, 2, 5})
    {
        const double fwhm = PeakPickerIM::aggregationRtFwhm(grid, k, cutoff, 3.0);
        TEST_EQUAL(frames_pooled_per_side(grid, 20, fwhm, cutoff), static_cast<size_t>(k))
    }

    // the same count at a different spacing, and with the retention times unsorted
    std::vector<double> coarse;
    for (int i = 0; i < 40; ++i) coarse.push_back(i * 1.786);
    std::vector<double> reversed(coarse.rbegin(), coarse.rend());
    const double fwhm_coarse = PeakPickerIM::aggregationRtFwhm(reversed, 1, cutoff, 3.0);
    TEST_REAL_SIMILAR(fwhm_coarse, 1.5 * 1.786 / half_width_per_fwhm)
    TEST_EQUAL(frames_pooled_per_side(coarse, 20, fwhm_coarse, cutoff), 1u)

    // jitter (gaps of 0.737 to 0.914 s) does not change the count: the median gap sets the width
    std::vector<double> jittered;
    double t = 0.0;
    for (int i = 0; i < 40; ++i)
    {
        jittered.push_back(t);
        t += (i % 10 == 3) ? 0.914 : 0.738;
    }
    const double fwhm_jitter = PeakPickerIM::aggregationRtFwhm(jittered, 1, cutoff, 3.0);
    for (size_t centre = 1; centre + 1 < jittered.size(); ++centre)
    {
        TEST_EQUAL(frames_pooled_per_side(jittered, centre, fwhm_jitter, cutoff) >= 1u, true)
        TEST_EQUAL(frames_pooled_per_side(jittered, centre, fwhm_jitter, cutoff) <= 1u, true)
    }

    // spacing cannot be determined
    TEST_EXCEPTION(Exception::InvalidValue, PeakPickerIM::aggregationRtFwhm({1.0}, 1, cutoff, 3.0))
    TEST_EXCEPTION(Exception::InvalidValue, PeakPickerIM::aggregationRtFwhm({1.0, 1.0, 1.0}, 1, cutoff, 3.0))
    TEST_EXCEPTION(Exception::InvalidValue, PeakPickerIM::aggregationRtFwhm(grid, 1, 1.5, 3.0))
}
END_SECTION

START_SECTION(double aggregationRtFwhm(const std::vector<double>& rts) const)
{
    std::vector<double> grid;
    for (int i = 0; i < 10; ++i) grid.push_back(i * 0.738);

    const double half_width_per_fwhm = std::sqrt(-std::log(0.01) / (4.0 * std::log(2.0)));
    PeakPickerIM pp_im;
    TEST_REAL_SIMILAR(pp_im.aggregationRtFwhm(grid), 5.0)  // default: rt_FWHM as given

    Param p = pp_im.getParameters();
    p.setValue("aggregation:frames_per_side", 1);
    pp_im.setParameters(p);
    TEST_REAL_SIMILAR(pp_im.aggregationRtFwhm(grid), 1.5 * 0.738 / half_width_per_fwhm)
}
END_SECTION

END_TEST