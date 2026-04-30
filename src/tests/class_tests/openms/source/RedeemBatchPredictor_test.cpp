// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/RedeemBatchPredictor.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>

using namespace OpenMS;
using namespace std;

START_TEST(RedeemBatchPredictor, "$Id$")

START_SECTION(static String canonicalizeModifiedPeptide(const String& modified_peptide))
{
  TEST_STRING_EQUAL(
    RedeemBatchPredictor::canonicalizeModifiedPeptide("PEPC(UniMod:4)PEPM(UniMod:35)PEPR"),
    "PEPC(UniMod:4)PEPM(UniMod:35)PEPR")

  TEST_STRING_EQUAL(
    RedeemBatchPredictor::canonicalizeModifiedPeptide(".(UniMod:1)PEPTIDE"),
    "(UniMod:1)PEPTIDE")

  TEST_EXCEPTION(
    Exception::InvalidValue,
    RedeemBatchPredictor::canonicalizeModifiedPeptide("PEPTIDE.(UniMod:2)"))
}
END_SECTION

START_SECTION(static double ccsToMobilityBruker(double ccs, Int charge, double precursor_mz))
{
  const double ccs = 400.0;
  const Int charge = 2;
  const double precursor_mz = 500.0;
  const double reduced_mass = (precursor_mz * charge * 28.0) / ((precursor_mz * charge) + 28.0);
  const double expected = ccs * std::sqrt(reduced_mass) / (charge * 1059.62245);

  TEST_REAL_SIMILAR(
    RedeemBatchPredictor::ccsToMobilityBruker(ccs, charge, precursor_mz),
    expected)
}
END_SECTION

START_SECTION([EXTRA] batch prediction returns RT optional CCS MS2 matrix and annotated fragments)
{
  const String rt_model = File::getOpenMSDataPath() + "/redeem/models/properties/20251205_100_epochs_min_max_rt_cnn_tf.safetensors";
  const String ccs_model = File::getOpenMSDataPath() + "/redeem/models/properties/20251205_500_epochs_early_stopped_100_min_max_ccs_cnn_tf.safetensors";
  const String ms2_model = File::getOpenMSDataPath() + "/redeem/models/properties/ms2.pth";

  if (!File::exists(rt_model) || !File::exists(ccs_model) || !File::exists(ms2_model))
  {
    STATUS("Skipping RedeemBatchPredictor integration test because slim ReDeem model files are not present under share/OpenMS/redeem/models/properties.");
    TEST_EQUAL(true, true)
  }
  else
  {
    RedeemBatchPredictor predictor;
    vector<RedeemPrecursorInput> inputs{
      {"PEPTIDE", 2, std::nullopt, std::nullopt},
      {"MGC(UniMod:4)AAR", 3, 27, String("QE")}
    };

    const auto predictions = predictor.predict(inputs);
    TEST_EQUAL(predictions.size(), 2)

    const auto& peptide_prediction = predictions[0];
    TEST_EQUAL(peptide_prediction.ccs.has_value(), true)
    TEST_EQUAL(peptide_prediction.mobility_1k0.has_value(), true)
    TEST_EQUAL(peptide_prediction.ms2_matrix.size(), 6)
    TEST_EQUAL(peptide_prediction.annotated_fragments.size(), 24)
    TEST_REAL_SIMILAR(peptide_prediction.precursor_mz, AASequence::fromString("PEPTIDE").getMZ(2))

    const auto& first_fragment = peptide_prediction.annotated_fragments.front();
    TEST_STRING_EQUAL(first_fragment.ion_type, "b")
    TEST_EQUAL(first_fragment.ordinal, 1)
    TEST_EQUAL(first_fragment.charge, 1)
    TEST_REAL_SIMILAR(first_fragment.mz, AASequence::fromString("PEPTIDE").getPrefix(1).getMZ(1, Residue::BIon))

    RedeemBatchPredictor::Config no_ccs_config;
    no_ccs_config.enable_ccs = false;
    RedeemBatchPredictor no_ccs_predictor(no_ccs_config);
    const auto no_ccs_predictions = no_ccs_predictor.predict({inputs.front()});
    TEST_EQUAL(no_ccs_predictions.size(), 1)
    TEST_EQUAL(no_ccs_predictions.front().ccs.has_value(), false)
    TEST_EQUAL(no_ccs_predictions.front().mobility_1k0.has_value(), false)
  }
}
END_SECTION

END_TEST
