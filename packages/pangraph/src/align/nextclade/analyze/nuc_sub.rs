use crate::align::nextclade::alphabet::letter::Letter;
use crate::align::nextclade::alphabet::nuc::{from_nuc, Nuc};
use crate::align::nextclade::analyze::abstract_mutation::{AbstractMutation, MutParams, Pos, QryLetter, RefLetter};
use crate::align::nextclade::analyze::nuc_del::NucDel;
use crate::align::nextclade::coord::position::NucRefGlobalPosition;
use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize, Hash)]
#[serde(rename_all = "camelCase")]
pub struct NucSub {
  pub pos: NucRefGlobalPosition,
  pub ref_nuc: Nuc,
  pub qry_nuc: Nuc,
}

impl AbstractMutation<NucRefGlobalPosition, Nuc> for NucSub {
  fn clone_with(&self, params: MutParams<NucRefGlobalPosition, Nuc>) -> Self {
    Self {
      pos: params.pos,
      ref_nuc: params.ref_letter,
      qry_nuc: params.qry_letter,
    }
  }
}

impl QryLetter<Nuc> for NucSub {
  fn qry_letter(&self) -> Nuc {
    self.qry_nuc
  }
}

impl RefLetter<Nuc> for NucSub {
  fn ref_letter(&self) -> Nuc {
    self.ref_nuc
  }
}

impl Pos<NucRefGlobalPosition> for NucSub {
  fn pos(&self) -> NucRefGlobalPosition {
    self.pos
  }
}

impl NucSub {
  /// Checks whether this substitution is a deletion (substitution of letter `Gap`)
  pub fn is_del(&self) -> bool {
    self.qry_nuc.is_gap()
  }

  #[must_use]
  pub const fn invert(&self) -> Self {
    Self {
      ref_nuc: self.qry_nuc,
      pos: self.pos,
      qry_nuc: self.ref_nuc,
    }
  }
}

impl Display for NucSub {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    // NOTE: by convention, in bioinformatics, nucleotides are numbered starting from 1, however our arrays are 0-based
    write!(
      f,
      "{}{}{}",
      from_nuc(self.ref_nuc),
      self.pos + 1,
      from_nuc(self.qry_nuc)
    )
  }
}

impl From<&NucDel> for NucSub {
  fn from(del: &NucDel) -> Self {
    Self {
      pos: del.pos,
      ref_nuc: del.ref_nuc,
      qry_nuc: Nuc::Gap,
    }
  }
}
