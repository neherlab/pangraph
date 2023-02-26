import { useCallback, useEffect, useRef, useState } from 'react'
import { useRecoilValue } from 'recoil'
import styled from 'styled-components'
import { cameraAtom } from 'src/state/camera.state'
import { useViewport } from 'src/hooks/useViewport'
import { RenderStats } from './gl/FrameCounter'
import { Renderer } from './Renderer'

export const Viewport = styled.div`
  width: 100%;
  height: 100%;
  overflow: hidden;
`

export const Canvas = styled.canvas`
  width: 100%;
  height: 100%;
  overflow: hidden;
`

export interface RenderSettings {
  clearColor: string
  wireframe: boolean
}

const renderSettingsDefaults: RenderSettings = {
  clearColor: '#111',
  wireframe: false,
}

export interface SceneProps {
  name: string
  color: string
}

export default function Scene({ name, color }: SceneProps) {
  const [renderer, setRenderer] = useState<Renderer | null>(null)
  const { cameraPreset, cameraPresets } = useRecoilValue(cameraAtom)

  const [renderSettings] = useState(renderSettingsDefaults)
  const viewportElem = useRef<HTMLDivElement>(null)
  const viewport = useViewport({ targetRef: viewportElem })

  const [, setStats] = useState<RenderStats>({ fps: 0, frameTime: 0, ups: 0, updateTime: 0 })

  const onCanvasRef = useCallback(
    (canvasNode: HTMLCanvasElement) => {
      if (renderer && !canvasNode) {
        renderer?.destroy()
        setRenderer(null)
      }

      if (renderer || !canvasNode) {
        return
      }

      const cameraPreset = cameraPresets.Lateral
      const newRendererGeo = new Renderer(canvasNode, renderSettings, cameraPreset, color, setStats)

      setRenderer(newRendererGeo)
      newRendererGeo.startRenderLoop(0)
    },
    [renderer], // eslint-disable-line react-hooks/exhaustive-deps
  )

  /** Handle viewport resize */
  useEffect(() => renderer?.onResize(viewport), [renderer, viewport])

  useEffect(() => renderer?.onCameraPreset(cameraPreset), [renderer, cameraPreset])

  return (
    <Viewport ref={viewportElem}>
      <Canvas ref={onCanvasRef} />
    </Viewport>
  )
}
