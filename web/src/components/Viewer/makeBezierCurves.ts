import { CubicBezierCurve3, Vector2, Vector3 } from 'three'
import { Line2, LineGeometry, LineMaterial } from 'src/vendor/three-fatline/src'
import Color from 'color'

export interface Curve {
  color: string
  points: [Vector3, Vector3, Vector3, Vector3]
}

const CURVES: Curve[] = [
  {
    color: '#f00',
    points: [
      // prettier-ignore
      new Vector3(-10, 0, 0),
      new Vector3(-5, 15, 0),
      new Vector3(20, 15, 0),
      new Vector3(10, 0, 0),
    ],
  },
  {
    color: '#00f',
    points: [
      // prettier-ignore
      new Vector3(-10, 0, 0),
      new Vector3(0, -5, 5),
      new Vector3(0, -8, 5),
      new Vector3(10, 0, 0),
    ],
  },
  {
    color: '#0f0',
    points: [
      // prettier-ignore
      new Vector3(-10, 0, 0),
      new Vector3(0, 0, 2),
      new Vector3(0, 2, 0),
      new Vector3(10, 0, 0),
    ],
  },
  {
    color: '#0ff',
    points: [
      // prettier-ignore
      new Vector3(-10, 0, 0),
      new Vector3(-15, 5, 2),
      new Vector3(-20, 2, 5),
      new Vector3(-10, 0, 0),
    ],
  },
  {
    color: '#ff0',
    points: [
      // prettier-ignore
      new Vector3(10, 0, 0),
      new Vector3(20, 5, 2),
      new Vector3(15, 2, 5),
      new Vector3(10, 0, 0),
    ],
  },
]

export function makeBezierCurve(
  bezierPoints: [Vector3, Vector3, Vector3, Vector3],
  color: string,
  resolution: Vector2,
) {
  const points = new CubicBezierCurve3(...bezierPoints).getPoints(50)

  const geometry = new LineGeometry()
  geometry.setPositions(points.flatMap(({ x, y, z }) => [x, y, z]))

  const material = new LineMaterial({
    color: new Color(color).rgbNumber(),
    linewidth: 10,
    resolution,
  })

  const curve = new Line2(geometry, material)
  curve.geometry.computeBoundingSphere()

  return curve
}

export function makeBezierCurves(resolution: Vector2) {
  return CURVES.map((curve) => makeBezierCurve(curve.points, curve.color, resolution))
}
