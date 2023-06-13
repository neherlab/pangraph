import { AxesHelper, Object3D, PerspectiveCamera, Scene, Vector2, WebGLRenderer } from 'three'
import { Easing, Tween } from '@tweenjs/tween.js'
import { OrbitControls } from 'three-stdlib'
import type { RenderSettings } from 'src/components/Viewer/Scene'
import type { UseViewportResult } from 'src/hooks/useViewport'
import { CameraPreset } from 'src/state/camera.state'
import { FrameCounter, RenderStats } from './gl/FrameCounter'
import { makeCoordinatePlanes } from './gl/makePlane'
import { UserInput } from './gl/UserInput'
import { makeBezierCurves } from './makeBezierCurves'

const CAMERA_POS_Z = 50

export function getViewportSize(canvas: HTMLCanvasElement) {
  if (!canvas.parentElement) {
    throw new Error('Canvas element parent is invalid')
  }
  const rect = canvas.parentElement?.getBoundingClientRect()
  const { width, height } = rect
  return { width, height }
}

export class Renderer {
  private readonly renderer: WebGLRenderer
  private readonly scene: Scene
  private readonly camera: PerspectiveCamera
  private readonly controls: OrbitControls
  private readonly input: UserInput

  private cameraPresetTarget: CameraPreset | undefined
  private cameraTween: Tween<CameraPreset> | undefined

  private objects: Object3D[] = []

  private rafHandle: number | undefined
  private shouldRender = true

  public readonly frameCounter: FrameCounter

  public constructor(
    canvas: HTMLCanvasElement,
    renderSettings: RenderSettings,
    cameraPreset: CameraPreset,
    color: string,
    setRenderStats: (stats: RenderStats) => void,
  ) {
    if (!canvas) {
      throw new Error('Canvas element is invalid')
    }

    const { width, height } = getViewportSize(canvas)
    const aspect = width / height

    this.renderer = new WebGLRenderer({
      canvas,
      antialias: true,
      precision: 'highp',
      powerPreference: 'high-performance',
    })
    this.renderer.setClearColor(renderSettings.clearColor)
    this.renderer.setPixelRatio(window.devicePixelRatio)
    this.renderer.setSize(width, height)

    this.scene = new Scene()

    this.camera = new PerspectiveCamera(75, aspect, 0.1, 100_000)
    this.camera.position.z = CAMERA_POS_Z
    this.controls = new OrbitControls(this.camera, canvas)
    this.controls.minDistance = 10
    this.controls.maxDistance = 500
    this.controls.zoomSpeed = 2
    this.onCameraPreset(cameraPreset, false)

    this.frameCounter = new FrameCounter(setRenderStats)

    this.input = new UserInput(canvas)

    void this.init(new Vector2(width, height)) // eslint-disable-line no-void
  }

  public async init(resolution: Vector2) {
    const axesHelper = new AxesHelper(1000)
    this.scene.add(axesHelper)

    const coordPlanes = makeCoordinatePlanes()
    this.scene.add(...coordPlanes)

    this.objects = makeBezierCurves(resolution)
    this.scene.add(...this.objects)
  }

  public destroy() {
    if (this.rafHandle) {
      cancelAnimationFrame(this.rafHandle)
    }

    // this.objects.forEach(object => object.geometry.dispose())
    // this.material?.dispose()
    this.camera.clear()
    this.controls.dispose()
    this.scene.clear()
    this.scene.removeFromParent()
  }

  public onResize({ width, height }: UseViewportResult) {
    if (width === 0 || height === 0) {
      return
    }
    this.camera.aspect = width / height
    this.camera.updateProjectionMatrix()
    this.renderer.setSize(width, height)
  }

  public onCameraPreset(cameraPreset: CameraPreset, animated = true) {
    if (animated) {
      this.cameraTween = new Tween<CameraPreset>({
        polar: this.controls.getPolarAngle(),
        azimuthal: this.controls.getAzimuthalAngle(),
      })
        .to(cameraPreset, 250)
        .easing(Easing.Quadratic.Out)
        .onComplete((_) => {
          this.cameraPresetTarget = undefined
        })
        .onUpdate((cameraPreset) => {
          this.cameraPresetTarget = cameraPreset
        })
        .start()
    } else {
      const { polar, azimuthal } = cameraPreset
      this.controls.setPolarAngle(polar)
      this.controls.setAzimuthalAngle(azimuthal)
    }
  }

  public startRenderLoop(_time: number) {
    const { deltaMs } = this.frameCounter.update()
    this.cameraTween?.update()

    this.update(deltaMs)
    this.render(deltaMs)

    this.rafHandle = requestAnimationFrame((time) => this.startRenderLoop(time))
  }

  public update(_deltaMs: number) {
    if (this.cameraPresetTarget) {
      this.controls.setPolarAngle(this.cameraPresetTarget.polar)
      this.controls.setAzimuthalAngle(this.cameraPresetTarget.azimuthal)
    }

    this.controls.update()
    this.camera.updateProjectionMatrix()

    if (this.input.isCtrlDown) {
      this.controls.enableRotate = false
      this.controls.enablePan = false
    } else {
      this.controls.enableRotate = true
      this.controls.enablePan = true
    }
  }

  public render(_deltaMs: number) {
    if (this.shouldRender) {
      this.frameCounter.frame()
      this.renderer.clear()
      this.renderer.render(this.scene, this.camera)
    }
  }
}
