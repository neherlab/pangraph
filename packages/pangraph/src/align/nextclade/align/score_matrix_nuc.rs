use crate::align::nextclade::alphabet::nuc::Nuc;

const NUM_COLS: usize = 16;
const SCORING_MATRIX_NUC_SIZE: usize = NUM_COLS * NUM_COLS;

#[rustfmt::skip]
static SCORING_MATRIX_NUC: &[i32; SCORING_MATRIX_NUC_SIZE] = &[
  /*           01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16 */
  /*            T   A   W   C   Y   M   H   G   K   R   D   S   B   V   N   - */
  /* 01   T */  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,
  /* 02   A */  0,  1,  1,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  1,  1,  0,
  /* 03   W */  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1,  0,
  /* 04   C */  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  1,  1,  1,  1,  0,
  /* 05   Y */  1,  0,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  1,  1,  1,  0,
  /* 06   M */  0,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1,  0,
  /* 07   H */  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 08   G */  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 09   K */  1,  0,  1,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 10   R */  0,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 11   D */  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 12   S */  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 13   B */  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 14   V */  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
  /* 15   N */  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 16   - */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,
  ];

pub fn lookup_nuc_scoring_matrix(x: Nuc, y: Nuc) -> i32 {
  SCORING_MATRIX_NUC[x as usize * NUM_COLS + y as usize]
}
