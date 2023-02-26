import { RefObject } from 'react'
import { useMouse } from 'react-use'

/**
 * Get relative position of mouse cursor over the element given by `ref`
 */
export function useMouseRel<T extends Element>(ref: RefObject<T>) {
  const { elX, elY, elW, elH } = useMouse(ref)

  return {
    x: elX / elW - 0.5,
    y: elY / elH - 0.5,
  }
}
