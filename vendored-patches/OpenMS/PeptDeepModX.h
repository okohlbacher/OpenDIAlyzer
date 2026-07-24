// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// Vendored patch P2 for OpenDIAlyzer (github.com/okohlbacher/OpenDIAlyzer):
// modification featurization for the PeptDeep ONNX predictors. Stock OpenMS
// zero-fills mod_x (unmodified only); this fills it from an AASequence's
// modifications, matching AlphaPeptDeep (peptdeep 1.5.0) exactly.
//
// The 109-element vocabulary order and the [batch, nAA+2, 109] layout
// (N-term = position 0, residue r (0-based) = position r+1, C-term = position
// nAA+1) were extracted from peptdeep and are pinned by
// src/testdata/peptdeep_modx_reference.json. Header-only so it can be validated
// without rebuilding libOpenMS.
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
  namespace ML
  {
    /// AlphaPeptDeep mod-feature element vocabulary, exact order
    /// (peptdeep.settings.model_const['mod_elements']). Index into mod_x's last
    /// dimension. Length == PEPTDEEP_MOD_ELEMENTS (109).
    inline const std::array<const char*, 109>& peptDeepModElements()
    {
      static const std::array<const char*, 109> elements = {{
    "C", "H", "N", "O", "P", "S", "B", "F", "I", "K",
    "U", "V", "W", "X", "Y", "Ac", "Ag", "Al", "Am", "Ar",
    "As", "At", "Au", "Ba", "Be", "Bi", "Bk", "Br", "Ca", "Cd",
    "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Dy", "Er",
    "Es", "Eu", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "He", "Hf",
    "Hg", "Ho", "In", "Ir", "Kr", "La", "Li", "Lr", "Lu", "Md",
    "Mg", "Mn", "Mo", "Na", "Nb", "Nd", "Ne", "Ni", "No", "Np",
    "Os", "Pa", "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra",
    "Rb", "Re", "Rh", "Rn", "Ru", "Sb", "Sc", "Se", "Si", "Sm",
    "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Th", "Ti", "Tl", "Tm",
    "Xe", "Yb", "Zn", "Zr", "2H", "13C", "15N", "18O", "?",
      }};
      return elements;
    }

    /// Element symbol -> vocabulary index. Unknown symbols map to '?' (index 108),
    /// matching AlphaPeptDeep's catch-all slot.
    inline int peptDeepModElementIndex(const std::string& symbol)
    {
      static const std::unordered_map<std::string, int> lut = [] {
        std::unordered_map<std::string, int> m;
        const auto& e = peptDeepModElements();
        for (int i = 0; i < static_cast<int>(e.size()); ++i) m[e[i]] = i;
        return m;
      }();
      auto it = lut.find(symbol);
      return it != lut.end() ? it->second : 108; // '?'
    }

    /// @brief Fill the mod_x slice for one peptide from its AASequence modifications.
    ///
    /// @param seq the (possibly modified) peptide
    /// @param add_terminal_tokens whether aa_indices reserves position 0 for the
    ///        N-term token (must match the aa_indices layout; PeptDeep uses true)
    /// @param sequence_length padded length of the mod_x slice (== aa_indices row length)
    /// @param out pointer to a sequence_length*109 float slice; ZEROED then filled
    ///
    /// Placement (verified against peptdeep): N-terminal mod -> position 0,
    /// residue r (0-based) -> position r+1, C-terminal mod -> position nAA+1.
    /// Values are the modification diff-formula element counts (additive if
    /// several mods share a position).
    inline void fillPeptideModX(const AASequence& seq, bool add_terminal_tokens,
                                size_t sequence_length, float* out)
    {
      const size_t E = static_cast<size_t>(PEPTDEEP_MOD_ELEMENTS);
      std::fill(out, out + sequence_length * E, 0.0f);
      const size_t offset = add_terminal_tokens ? 1 : 0; // residue r -> position r+offset
      const size_t nAA = seq.size();

      auto add_mod = [&](size_t position, const ResidueModification* mod) {
        if (mod == nullptr || position >= sequence_length) return;
        for (const auto& [symbol, count] : mod->getDiffFormula().toMap())
          out[position * E + peptDeepModElementIndex(symbol)] += static_cast<float>(count);
      };

      // N-terminal modification -> position 0 (only meaningful with a terminal token).
      if (add_terminal_tokens) add_mod(0, seq.getNTerminalModification());

      for (size_t r = 0; r < nAA; ++r)
        if (seq[r].isModified()) add_mod(r + offset, seq[r].getModification());

      // C-terminal modification -> the C-term token position (nAA+1 with tokens).
      if (add_terminal_tokens) add_mod(nAA + 1, seq.getCTerminalModification());
    }

  } // namespace ML
} // namespace OpenMS
