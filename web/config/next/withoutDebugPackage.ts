import path from 'path'
import type { NextConfig } from 'next'
import { addWebpackLoader } from './lib/addWebpackLoader'

const THIS_DIR = new URL('.', import.meta.url).pathname

export default function withoutDebugPackage(nextConfig: NextConfig) {
  return addWebpackLoader(nextConfig, (_webpackConfig, _context) => ({
    test: /\.(ts|tsx|ctx|mtx|js|jsx|cjs|mjs|)$/i,
    use: [
      {
        loader: path.resolve(THIS_DIR, 'loaders', 'removeDebugPackageLoader.cjs'),
      },
    ],
  }))
}
