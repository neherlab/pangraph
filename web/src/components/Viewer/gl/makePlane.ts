import { DoubleSide, Mesh, MeshBasicMaterial, PlaneGeometry } from 'three'

export interface PlaneParams {
  width: number
  height: number
  x?: number
  y?: number
  z?: number
  rotX?: number
  rotY?: number
  rotZ?: number
  color: string
}

export function makePlane({ width, height, x, y, z, rotX, rotY, rotZ, color }: PlaneParams) {
  const plane = new Mesh(new PlaneGeometry(width, height), new MeshBasicMaterial({ color, side: DoubleSide }))
  plane.position.x = x ?? 0
  plane.position.y = y ?? 0
  plane.position.z = z ?? 0
  plane.rotation.x = rotX ?? 0
  plane.rotation.y = rotY ?? 0
  plane.rotation.z = rotZ ?? 0

  plane.geometry.boundingBox = null
  plane.geometry.boundingSphere = null

  plane.name = 'plane'

  return plane
}

export function makeCoordinatePlanes() {
  const colorX = '#200'
  const colorY = '#020'
  const colorZ = '#002'
  const angle = Math.PI * 0.5
  const width = 1000
  const height = 1000

  return [
    makePlane({ width, height, x: -500, rotY: angle, rotZ: angle, color: colorX }),
    makePlane({ width, height, x: +500, rotY: angle, rotZ: angle, color: colorX }),

    makePlane({ width, height, y: -500, rotX: angle, rotZ: angle, color: colorY }),
    makePlane({ width, height, y: +500, rotX: angle, rotZ: angle, color: colorY }),

    makePlane({ width, height, z: -500, color: colorZ }),
    makePlane({ width, height, z: +500, color: colorZ }),
  ]
}
