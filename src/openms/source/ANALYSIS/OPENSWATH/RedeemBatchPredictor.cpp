// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/RedeemBatchPredictor.h>

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  namespace
  {
    constexpr double IM_GAS_MASS = 28.0;
    constexpr double CCS_IM_COEF = 1059.62245;

    const char* DEFAULT_RT_MODEL = "redeem/models/properties/20251205_100_epochs_min_max_rt_cnn_tf.safetensors";
    const char* DEFAULT_CCS_MODEL = "redeem/models/properties/20251205_500_epochs_early_stopped_100_min_max_ccs_cnn_tf.safetensors";
    const char* DEFAULT_MS2_MODEL = "redeem/models/properties/ms2.pth";

    struct OpenMsRedeemPredictorConfigFFI
    {
      const char* rt_model_path;
      const char* ccs_model_path;
      const char* ms2_model_path;
      const char* device_preference;
    };

    struct OpenMsRedeemPredictionInputFFI
    {
      const char* modified_peptide;
      int precursor_charge;
      int nce;
      const char* instrument;
    };

    struct OpenMsRedeemBatchOutputFFI
    {
      size_t count;
      int has_ccs;
      float* rt_values;
      float* ccs_values;
      size_t* ms2_row_counts;
      float* ms2_values;
      size_t ms2_value_count;
    };

    extern "C"
    {
      void* openms_redeem_predictor_create(const OpenMsRedeemPredictorConfigFFI* config);
      void openms_redeem_predictor_destroy(void* predictor);
      int openms_redeem_predict_batch(
        const void* predictor,
        const OpenMsRedeemPredictionInputFFI* inputs,
        size_t input_count,
        OpenMsRedeemBatchOutputFFI* output);
      void openms_redeem_batch_output_free(OpenMsRedeemBatchOutputFFI* output);
      const char* openms_redeem_last_error();
    }

    struct PreparedInput
    {
      AASequence sequence;
      std::string canonical_modified_peptide;
      std::string instrument;
      Int precursor_charge = 0;
      Int nce = 30;
      double precursor_mz = 0.0;
    };

    String lastRedeemError_()
    {
      const char* error = openms_redeem_last_error();
      return error == nullptr ? String("redeem returned an unknown error") : String(error);
    }

    [[noreturn]] void throwRedeemError_(const String& context, const String& detail = "")
    {
      String message = context;
      if (!detail.empty())
      {
        message += ": " + detail;
      }
      throw Exception::InvalidValue(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        lastRedeemError_(),
        message);
    }

    String resolveModelPath_(const std::optional<String>& explicit_path, const char* default_relative_path)
    {
      return explicit_path.has_value() ? *explicit_path : File::find(default_relative_path);
    }

    String canonicalizeModifiedPeptide_(const AASequence& sequence, const String& original)
    {
      if (sequence.hasCTerminalModification())
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          original,
          "C-terminal peptide modifications are not supported by redeem v1.");
      }

      String canonical = sequence.toUniModString();
      if (!canonical.empty() && canonical[0] == '.')
      {
        canonical.erase(0, 1);
      }
      if (canonical.hasSubstring(".(") || canonical.hasSubstring(".["))
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          original,
          "Unexpected terminal dot notation remained after canonicalization.");
      }
      return canonical;
    }

    PreparedInput prepareInput_(const RedeemPrecursorInput& input)
    {
      AASequence sequence = AASequence::fromString(input.modified_peptide);
      String canonical = canonicalizeModifiedPeptide_(sequence, input.modified_peptide);

      if (input.precursor_charge <= 0)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          String(input.precursor_charge),
          "Redeem precursor_charge must be positive.");
      }

      PreparedInput prepared;
      prepared.sequence = std::move(sequence);
      prepared.canonical_modified_peptide = canonical;
      prepared.instrument = input.instrument.has_value() ? std::string(*input.instrument) : std::string("Lumos");
      prepared.precursor_charge = input.precursor_charge;
      prepared.nce = input.nce.has_value() ? *input.nce : 30;
      prepared.precursor_mz = prepared.sequence.getMZ(prepared.precursor_charge);
      return prepared;
    }

    double reducedMass_(double precursor_mz, Int charge)
    {
      const double ion_mass = precursor_mz * static_cast<double>(charge);
      return ion_mass * IM_GAS_MASS / (ion_mass + IM_GAS_MASS);
    }

    std::vector<RedeemFragmentPrediction> annotateFragments_(
      const AASequence& sequence,
      const std::vector<std::array<float, 8>>& matrix)
    {
      std::vector<RedeemFragmentPrediction> fragments;

      const Size peptide_length = sequence.size();
      if (matrix.empty())
      {
        return fragments;
      }
      if (peptide_length < 2)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          sequence.toString(),
          "Cannot annotate MS2 fragments for peptides shorter than length 2.");
      }
      const Size expected_rows = peptide_length - 1;
      if (matrix.size() < expected_rows)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          String(matrix.size()),
          "Redeem MS2 row count was smaller than peptide_length - 1.");
      }
      fragments.reserve(expected_rows * 4);

      // Some model outputs include a trailing padded row. Fragment annotation only
      // uses the first peptide_length - 1 cleavage positions.
      for (Size row = 0; row < expected_rows; ++row)
      {
        const Size cleavage = row + 1;
        const Size y_ordinal = peptide_length - cleavage;
        const AASequence b_fragment = sequence.getPrefix(cleavage);
        const AASequence y_fragment = sequence.getSuffix(y_ordinal);

        for (Int charge = 1; charge <= 2; ++charge)
        {
          fragments.push_back(
            RedeemFragmentPrediction{
              "b",
              cleavage,
              charge,
              b_fragment.getMZ(charge, Residue::BIon),
              matrix[row][static_cast<Size>(charge - 1)]
            });
        }

        for (Int charge = 1; charge <= 2; ++charge)
        {
          fragments.push_back(
            RedeemFragmentPrediction{
              "y",
              y_ordinal,
              charge,
              y_fragment.getMZ(charge, Residue::YIon),
              matrix[row][static_cast<Size>(charge + 1)]
            });
        }
      }

      return fragments;
    }
  }

  struct RedeemBatchPredictor::Impl
  {
    void* predictor = nullptr;
  };

  RedeemBatchPredictor::RedeemBatchPredictor()
    : RedeemBatchPredictor(Config())
  {
  }

  RedeemBatchPredictor::RedeemBatchPredictor(Config config)
    : impl_(std::make_unique<Impl>())
  {
    const String rt_model_path = resolveModelPath_(config.rt_model_path, DEFAULT_RT_MODEL);
    const String ms2_model_path = resolveModelPath_(config.ms2_model_path, DEFAULT_MS2_MODEL);
    const std::optional<String> ccs_model_path = config.enable_ccs
      ? std::optional<String>(config.ccs_model_path.has_value() ? *config.ccs_model_path : File::find(DEFAULT_CCS_MODEL))
      : std::nullopt;

    OpenMsRedeemPredictorConfigFFI ffi_config
    {
      rt_model_path.c_str(),
      ccs_model_path.has_value() ? ccs_model_path->c_str() : nullptr,
      ms2_model_path.c_str(),
      config.device_preference.c_str()
    };

    impl_->predictor = openms_redeem_predictor_create(&ffi_config);
    if (impl_->predictor == nullptr)
    {
      throwRedeemError_("Failed to create Redeem predictor");
    }
  }

  RedeemBatchPredictor::~RedeemBatchPredictor()
  {
    if (impl_ && impl_->predictor != nullptr)
    {
      openms_redeem_predictor_destroy(impl_->predictor);
      impl_->predictor = nullptr;
    }
  }

  RedeemBatchPredictor::RedeemBatchPredictor(RedeemBatchPredictor&&) noexcept = default;
  RedeemBatchPredictor& RedeemBatchPredictor::operator=(RedeemBatchPredictor&&) noexcept = default;

  std::vector<RedeemPrediction> RedeemBatchPredictor::predict(const std::vector<RedeemPrecursorInput>& inputs) const
  {
    if (impl_ == nullptr || impl_->predictor == nullptr)
    {
      throw Exception::InvalidValue(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "null predictor",
        "Redeem predictor has not been initialized.");
    }
    if (inputs.empty())
    {
      return {};
    }

    std::vector<PreparedInput> prepared_inputs;
    prepared_inputs.reserve(inputs.size());
    for (const auto& input : inputs)
    {
      prepared_inputs.push_back(prepareInput_(input));
    }

    std::vector<OpenMsRedeemPredictionInputFFI> ffi_inputs(prepared_inputs.size());
    for (Size i = 0; i < prepared_inputs.size(); ++i)
    {
      ffi_inputs[i] = OpenMsRedeemPredictionInputFFI
      {
        prepared_inputs[i].canonical_modified_peptide.c_str(),
        prepared_inputs[i].precursor_charge,
        prepared_inputs[i].nce,
        prepared_inputs[i].instrument.c_str()
      };
    }

    OpenMsRedeemBatchOutputFFI ffi_output
    {
      0,
      0,
      nullptr,
      nullptr,
      nullptr,
      nullptr,
      0
    };

    struct OutputGuard
    {
      OpenMsRedeemBatchOutputFFI* output;
      ~OutputGuard()
      {
        if (output != nullptr)
        {
          openms_redeem_batch_output_free(output);
        }
      }
    } output_guard{&ffi_output};

    const int status = openms_redeem_predict_batch(
      impl_->predictor,
      ffi_inputs.data(),
      ffi_inputs.size(),
      &ffi_output);

    if (status != 1)
    {
      throwRedeemError_("Redeem batch prediction failed");
    }
    if (ffi_output.count != prepared_inputs.size())
    {
      throw Exception::InvalidValue(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        String(ffi_output.count),
        "Redeem batch prediction returned an unexpected number of outputs.");
    }

    const float* rt_values = ffi_output.rt_values;
    const float* ccs_values = ffi_output.ccs_values;
    const size_t* row_counts = ffi_output.ms2_row_counts;
    const float* flat_ms2 = ffi_output.ms2_values;

    std::vector<RedeemPrediction> predictions;
    predictions.reserve(prepared_inputs.size());

    Size flat_offset = 0;
    for (Size i = 0; i < prepared_inputs.size(); ++i)
    {
      const Size row_count = row_counts[i];
      std::vector<std::array<float, 8>> matrix;
      matrix.reserve(row_count);

      for (Size row = 0; row < row_count; ++row)
      {
        std::array<float, 8> row_values {};
        for (Size channel = 0; channel < row_values.size(); ++channel)
        {
          if (flat_offset >= ffi_output.ms2_value_count)
          {
            throw Exception::InvalidValue(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              String(flat_offset),
              "Redeem returned a truncated flattened MS2 tensor.");
          }
          row_values[channel] = flat_ms2[flat_offset++];
        }
        matrix.push_back(row_values);
      }

      RedeemPrediction prediction;
      prediction.precursor_mz = prepared_inputs[i].precursor_mz;
      prediction.rt = rt_values[i];
      prediction.ms2_matrix = matrix;
      prediction.annotated_fragments = annotateFragments_(prepared_inputs[i].sequence, prediction.ms2_matrix);

      if (ffi_output.has_ccs == 1 && ccs_values != nullptr)
      {
        prediction.ccs = ccs_values[i];
        prediction.mobility_1k0 = ccsToMobilityBruker(*prediction.ccs, prepared_inputs[i].precursor_charge, prediction.precursor_mz);
      }

      predictions.push_back(std::move(prediction));
    }

    if (flat_offset != ffi_output.ms2_value_count)
    {
      throw Exception::InvalidValue(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        String(flat_offset),
        "Redeem returned extra flattened MS2 values beyond the reconstructed matrices.");
    }

    return predictions;
  }

  String RedeemBatchPredictor::canonicalizeModifiedPeptide(const String& modified_peptide)
  {
    const AASequence sequence = AASequence::fromString(modified_peptide);
    return canonicalizeModifiedPeptide_(sequence, modified_peptide);
  }

  double RedeemBatchPredictor::ccsToMobilityBruker(double ccs, Int charge, double precursor_mz)
  {
    return ccs * std::sqrt(reducedMass_(precursor_mz, charge)) / (static_cast<double>(charge) * CCS_IM_COEF);
  }
}
