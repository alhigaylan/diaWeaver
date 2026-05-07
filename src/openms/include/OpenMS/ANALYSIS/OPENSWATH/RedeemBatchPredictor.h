// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace OpenMS
{
  struct OPENMS_DLLAPI RedeemPrecursorInput
  {
    String modified_peptide;
    Int precursor_charge = 0;
    std::optional<Int> nce;
    std::optional<String> instrument;
  };

  struct OPENMS_DLLAPI RedeemFragmentPrediction
  {
    String ion_type;
    Size ordinal = 0;
    Int charge = 0;
    double mz = 0.0;
    float intensity = 0.0f;
  };

  struct OPENMS_DLLAPI RedeemPrediction
  {
    double precursor_mz = 0.0;
    double rt = 0.0;
    std::optional<double> ccs;
    std::optional<double> mobility_1k0;
    std::vector<std::array<float, 8>> ms2_matrix;
    std::vector<RedeemFragmentPrediction> annotated_fragments;
  };

  struct OPENMS_DLLAPI RedeemFineTuneConfig
  {
    String training_tsv_path;
    std::optional<String> validation_tsv_path;
    double validation_fraction = 0.2;
    Size batch_size = 64;
    Size epochs = 10;
    Size early_stopping_patience = 5;
    double learning_rate = 1e-4;
    double warmup_fraction = 0.0;
    Int default_nce = 30;
    String default_instrument = "Lumos";
    bool enable_rt = true;
    bool enable_ccs = true;
    bool enable_ms2 = true;
    std::optional<String> rt_model_output_path;
    std::optional<String> ccs_model_output_path;
    std::optional<String> ms2_model_output_path;
  };

  class OPENMS_DLLAPI RedeemBatchPredictor
  {
  public:
    struct Config
    {
      std::optional<String> rt_model_path;
      std::optional<String> ccs_model_path;
      std::optional<String> ms2_model_path;
      bool enable_ccs = true;
      String device_preference = "cpu";
    };

    RedeemBatchPredictor();
    explicit RedeemBatchPredictor(Config config);
    ~RedeemBatchPredictor();

    RedeemBatchPredictor(const RedeemBatchPredictor&) = delete;
    RedeemBatchPredictor& operator=(const RedeemBatchPredictor&) = delete;
    RedeemBatchPredictor(RedeemBatchPredictor&&) noexcept;
    RedeemBatchPredictor& operator=(RedeemBatchPredictor&&) noexcept;

    std::vector<RedeemPrediction> predict(const std::vector<RedeemPrecursorInput>& inputs) const;
    void fineTuneFromTransitionTsv(const RedeemFineTuneConfig& config);

    static String canonicalizeModifiedPeptide(const String& modified_peptide);
    static double ccsToMobilityBruker(double ccs, Int charge, double precursor_mz);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}
