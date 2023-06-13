import { Vector2 } from 'three'

export type MouseEventNames = 'mousemove' | 'mousedown' | 'mouseup'

export type KeyboardEventNames = 'keydown' | 'keyup'

export class UserInput {
  private readonly mouse: Vector2
  private domElement: HTMLElement

  private mouseButtonDown = {
    left: false,
    middle: false,
    right: false,
  }

  private ctrlKey = false
  private shiftKey = false
  private altKey = false
  private metaKey = false

  private mouseHandlers = new Map<MouseEventNames, (event: MouseEvent) => void>([
    [
      'mousemove',
      (event: MouseEvent) => {
        const { left, top, width, height } = this.domElement.getBoundingClientRect()
        this.mouse.x = ((event.clientX - left) / width) * 2 - 1
        this.mouse.y = -((event.clientY - top) / height) * 2 + 1
      },
    ],
    [
      'mousedown',
      (event: MouseEvent) => {
        // prettier-ignore
        switch (event.button) {
          case 0: {
            this.mouseButtonDown.left = true;
            break
          }
          case 1: {
            this.mouseButtonDown.middle = true;
            break
          }
          case 2: {
            this.mouseButtonDown.right = true;
            break
          }
          default: {
            break
          }
        }
      },
    ],
    [
      'mouseup',
      (event: MouseEvent) => {
        // prettier-ignore
        switch (event.button) {
          case 0: {
            this.mouseButtonDown.left = false;
            break
          }
          case 1: {
            this.mouseButtonDown.middle = false;
            break
          }
          case 2: {
            this.mouseButtonDown.right = false;
            break
          }
          default: {
            break
          }
        }
      },
    ],
  ])

  private keyboardHandlers = new Map<KeyboardEventNames, (event: KeyboardEvent) => void>([
    [
      'keydown',
      (event: KeyboardEvent) => {
        if (event.ctrlKey) {
          this.ctrlKey = true
        }
        if (event.shiftKey) {
          this.shiftKey = true
        }
        if (event.altKey) {
          this.altKey = true
        }
        if (event.metaKey) {
          this.metaKey = true
        }
      },
    ],
    [
      'keyup',
      (event: KeyboardEvent) => {
        if (!event.ctrlKey) {
          this.ctrlKey = false
        }
        if (!event.shiftKey) {
          this.shiftKey = false
        }
        if (!event.altKey) {
          this.altKey = false
        }
        if (!event.metaKey) {
          this.metaKey = false
        }
      },
    ],
  ])

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private subscribe(handlers: any) {
    // eslint-disable-next-line prefer-const,no-loops/no-loops
    for (let [name, handler] of handlers.entries()) {
      handler = handler.bind(this)
      window.addEventListener(name, handler, false)
    }
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private unsubscribe(handlers: any) {
    // eslint-disable-next-line no-loops/no-loops
    for (const [name, handler] of this.mouseHandlers.entries()) {
      window.removeEventListener(name, handler, false)
    }
  }

  public constructor(domElement: HTMLElement) {
    this.subscribe(this.mouseHandlers)
    this.subscribe(this.keyboardHandlers)

    this.mouse = new Vector2(-1, -1)
    this.domElement = domElement
  }

  public destroy() {
    this.unsubscribe(this.mouseHandlers)
    this.unsubscribe(this.keyboardHandlers)
  }

  public get coords() {
    return this.mouse
  }

  public get isLeftDown() {
    return this.mouseButtonDown.left
  }

  public get isMiddleDown() {
    return this.mouseButtonDown.left
  }

  public get isRightDown() {
    return this.mouseButtonDown.left
  }

  public get isCtrlDown() {
    return this.ctrlKey
  }

  public get isShiftDown() {
    return this.shiftKey
  }

  public get isAltDown() {
    return this.altKey
  }

  public get isMetaDown() {
    return this.metaKey
  }
}
