use crate::align::alignment::AnchorBlock;
use crate::pangraph::pangraph_block::PangraphBlock;

pub fn select_anchor_block(block_qry: &PangraphBlock, block_ref: &PangraphBlock) -> AnchorBlock {
  if block_qry.depth() >= block_ref.depth() {
    AnchorBlock::Ref
  } else {
    AnchorBlock::Qry
  }
}
