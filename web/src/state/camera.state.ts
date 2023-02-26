import { atom } from 'recoil'

export interface CameraPreset {
  polar: number
  azimuthal: number
}

const cameraPresets: Record<string, CameraPreset> = {
  Lateral: { polar: Math.PI * 0.5, azimuthal: 0 },
  Ventral: { polar: Math.PI, azimuthal: 0 },
  Dorsal: { polar: 0, azimuthal: 0 },
  Anterior: { polar: Math.PI * 0.5, azimuthal: -Math.PI * 0.5 },
  Posterior: { polar: Math.PI * 0.5, azimuthal: Math.PI * 0.5 },
  Perspective: { polar: Math.PI * 0.4, azimuthal: -0.15 * Math.PI },
}

export const cameraAtom = atom({
  key: 'Camera',
  default: {
    cameraPreset: cameraPresets.Lateral,
    cameraPresets,
  },
})
