export function now(): number {
  return (performance || Date).now()
}

export interface RenderStats {
  fps: number
  frameTime: number
  ups: number
  updateTime: number
}

export class FrameCounter {
  private timeMsPrev = now()
  private nFrames = 0
  private nUpdates = 0
  private nFps = 0
  private nUps = 0

  public readonly setRenderStats: (stats: RenderStats) => void

  public constructor(setRenderStats: (stats: RenderStats) => void | undefined) {
    this.setRenderStats = setRenderStats
  }

  public update() {
    this.nUpdates += 1

    const timeMs = now()
    const deltaMs = timeMs - this.timeMsPrev

    if (deltaMs > 1000) {
      this.nFps = (this.nFrames * 1000) / deltaMs
      this.nUps = (this.nUpdates * 1000) / deltaMs

      this.nFrames = 0
      this.nUpdates = 0
      this.timeMsPrev = timeMs

      this.setRenderStats({
        fps: this.fps,
        frameTime: this.frameTime,
        ups: this.ups,
        updateTime: this.updateTime,
      })
    }

    return { deltaMs, timeMs }
  }

  public frame() {
    this.nFrames += 1
  }

  public get fps() {
    return this.nFps
  }

  public get frameTime() {
    return this.fps === 0 ? 0 : 1000 / this.fps
  }

  public get ups() {
    return this.nUps
  }

  public get updateTime() {
    return this.ups === 0 ? 0 : 1000 / this.ups
  }
}
