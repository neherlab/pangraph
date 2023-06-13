import path from 'path'
import { readJSONSync } from 'fs-extra/esm'
import ForkTsCheckerWebpackPlugin from 'fork-ts-checker-webpack-plugin'
import ESLintWebpackPlugin from 'eslint-webpack-plugin'
import type { NextConfig } from 'next'
import { addWebpackPlugin } from './lib/addWebpackPlugin'
import { findModuleRoot } from '../../lib/findModuleRoot'

const { moduleRoot } = findModuleRoot()

export interface GetWithTypeCheckingParams {
  eslint: boolean
  typeChecking: boolean
  memoryLimit?: number
  exclude?: string[]
}

const getWithTypeChecking =
  ({ eslint, typeChecking, memoryLimit = 512, exclude }: GetWithTypeCheckingParams) =>
  (nextConfig: NextConfig) => {
    let config = nextConfig

    if (!typeChecking && !eslint) {
      return config
    }

    if (typeChecking) {
      const tsConfig = readJSONSync('tsconfig.json')

      config = addWebpackPlugin(
        config,
        new ForkTsCheckerWebpackPlugin({
          issue: {
            exclude: exclude?.map((file) => ({ origin: 'typescript', file })),
          },

          typescript: {
            configFile: path.join(moduleRoot, 'tsconfig.json'),
            memoryLimit,
            mode: 'write-references',
            diagnosticOptions: { syntactic: true, semantic: true, declaration: true, global: true },
            configOverwrite: {
              compilerOptions: {
                ...tsConfig.compilerOptions,
                allowJs: false,
                skipLibCheck: true,
                sourceMap: false,
                inlineSourceMap: false,
                declarationMap: false,
                tsBuildInfoFile: '.cache/.tsbuildinfo.webpackplugin',
              },
              include: [
                'lib/**/*.js',
                'lib/**/*.jsx',
                'lib/**/*.ts',
                'lib/**/*.tsx',
                'src/**/*.js',
                'src/**/*.jsx',
                'src/**/*.ts',
                'src/**/*.tsx',
              ],
              exclude: [...tsConfig.exclude, ...(exclude ?? [])],
            },
          },

          formatter: 'codeframe',
        }),
      )
    }

    if (eslint) {
      config = addWebpackPlugin(
        config,
        new ESLintWebpackPlugin({
          threads: true,
          files: [path.join(moduleRoot, 'src/**/*.{js,jsx,ts,tsx}')],
          cache: false,
          formatter: 'codeframe',
        }),
      )
    }

    return config
  }

export default getWithTypeChecking
