/* eslint-disable no-loops/no-loops,no-param-reassign,no-plusplus */
import fs from 'fs-extra'
import path from 'path'

const THIS_DIR = new URL('.', import.meta.url).pathname

export function findModuleRoot(maxDepth = 10) {
  let moduleRoot = THIS_DIR
  while (--maxDepth) {
    moduleRoot = path.resolve(moduleRoot, '..')
    const file = path.join(moduleRoot, 'package.json')
    if (fs.existsSync(file)) {
      const pkg = fs.readJsonSync(file)
      return { moduleRoot, pkg }
    }
  }
  throw new Error('Module root not found')
}
