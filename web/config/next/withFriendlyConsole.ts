import FriendlyErrorsWebpackPlugin from '@nuxt/friendly-errors-webpack-plugin'
import type { NextConfig } from 'next'
import { addWebpackPlugin } from './lib/addWebpackPlugin'

interface FriendlyErrorsWebpackPluginError {
  message: string
  file: string
}

function cleanup() {
  return (error: FriendlyErrorsWebpackPluginError) => ({
    ...error,
    message: error.message.replace(/.*ERROR in.*\n/, '').replace(/.*WARNING in.*\n/, ''),
  })
}

function stripProjectRoot(projectRoot: string) {
  return (error: FriendlyErrorsWebpackPluginError) => ({
    ...error,
    message: error && error.message && error.message.replace(`${projectRoot}/`, ''),
    file: error && error.file && error.file.replace(`${projectRoot}/`, ''),
  })
}

export interface WithFriendlyConsoleParams {
  clearConsole: boolean
  projectRoot: string
}

const getWithFriendlyConsole =
  ({ clearConsole, projectRoot }: WithFriendlyConsoleParams) =>
  (nextConfig: NextConfig) => {
    return addWebpackPlugin(
      nextConfig,
      new FriendlyErrorsWebpackPlugin({
        clearConsole,
        additionalTransformers: [cleanup(), stripProjectRoot(projectRoot)],
        additionalFormatters: [],
      }),
    )
  }

export default getWithFriendlyConsole
