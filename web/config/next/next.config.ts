import type { NextConfig } from 'next'
import path from 'path'
import { uniq } from 'lodash-es'
import getWithMDX from '@next/mdx'
import remarkBreaks from 'remark-breaks'
import remarkMath from 'remark-math'
import remarkToc from 'remark-toc'
import remarkSlug from 'remark-slug'
import remarkImages from 'remark-images'

import { findModuleRoot } from '../../lib/findModuleRoot'
import { getGitBranch } from '../../lib/getGitBranch'
import { getBuildNumber } from '../../lib/getBuildNumber'
import { getBuildUrl } from '../../lib/getBuildUrl'
import { getGitCommitHash } from '../../lib/getGitCommitHash'
import { getEnvVars } from './lib/getEnvVars'

import getWithExtraWatch from './withExtraWatch'
import getWithFriendlyConsole from './withFriendlyConsole'
import { getWithRobotsTxt } from './withRobotsTxt'
import getWithTypeChecking from './withTypeChecking'
import withoutDebugPackage from './withoutDebugPackage'
import withSvg from './withSvg'
import withIgnore from './withIgnore'
import withoutMinification from './withoutMinification'
import withFriendlyChunkNames from './withFriendlyChunkNames'
import withResolve from './withResolve'
import withWebpackWatchPoll from './withWebpackWatchPoll'
import withUrlAsset from './withUrlAsset'

const {
  PRODUCTION,
  PROFILE,
  ENABLE_SOURCE_MAPS,
  ENABLE_ESLINT,
  ENABLE_TYPE_CHECKS,
  DOMAIN,
  DOMAIN_STRIPPED,
  WATCH_POLL,
  DATA_ROOT_URL,
} = getEnvVars()

const BRANCH_NAME = getGitBranch()

const { pkg, moduleRoot } = findModuleRoot()

const clientEnv = {
  BRANCH_NAME,
  PACKAGE_VERSION: pkg.version ?? '',
  BUILD_NUMBER: getBuildNumber(),
  TRAVIS_BUILD_WEB_URL: getBuildUrl(),
  COMMIT_HASH: getGitCommitHash(),
  DOMAIN,
  DOMAIN_STRIPPED,
  DATA_ROOT_URL,
}

const transpilationListDev = [
  // prettier-ignore
  'd3-scale',
]

const transpilationListProd = uniq([
  // prettier-ignore
  ...transpilationListDev,
  'debug',
  'lodash',
  'react-share',
  'recharts',
  'semver',
])

const nextConfig: NextConfig = {
  distDir: `.build/${process.env.NODE_ENV}/tmp`,
  pageExtensions: ['js', 'jsx', 'ts', 'tsx', 'md', 'mdx', 'all-contributorsrc'],
  onDemandEntries: {
    maxInactiveAge: 180 * 1000,
    pagesBufferLength: 5,
  },
  reactStrictMode: true,
  cleanDistDir: true,
  experimental: {
    legacyBrowsers: true,
    newNextLinkBehavior: true,
    scrollRestoration: true,
    swcMinify: true,
  },
  swcMinify: true,
  productionBrowserSourceMaps: ENABLE_SOURCE_MAPS,
  excludeDefaultMomentLocales: true,
  devIndicators: {
    buildActivity: false,
  },
  typescript: {
    ignoreBuildErrors: true,
  },
  eslint: {
    ignoreDuringBuilds: true,
  },
  compiler: {
    styledComponents: true,
  },
  env: clientEnv,
  poweredByHeader: false,
  webpack(config) {
    config.experiments.topLevelAwait = true
    return config
  },
  transpilePackages: PRODUCTION ? transpilationListProd : transpilationListDev,
  async rewrites() {
    return [{ source: '/:any*', destination: '/' }]
  },
}

const withMDX = getWithMDX({
  extension: /\.mdx?$/,
  options: {
    remarkPlugins: [
      remarkBreaks,
      remarkImages,
      remarkMath,
      remarkSlug,
      [remarkToc, { tight: true }],

      // [
      //   require('remark-autolink-headings'),
      //   {
      //     behavior: 'prepend',
      //     content: {
      //       type: 'element',
      //       tagName: 'i',
      //       properties: { className: ['bi', 'bi-link-45deg', 'mdx-link-icon'] },
      //     },
      //   },
      // ],
    ],
    rehypePlugins: [],
  },
})

const withFriendlyConsole = getWithFriendlyConsole({
  clearConsole: false,
  projectRoot: path.resolve(moduleRoot),
})

const withExtraWatch = getWithExtraWatch({
  files: [path.join(moduleRoot, 'src/types/**/*.d.ts')],
  dirs: [],
})

const withTypeChecking = getWithTypeChecking({
  typeChecking: ENABLE_TYPE_CHECKS,
  eslint: ENABLE_ESLINT,
  memoryLimit: 4096,
})

const withRobotsTxt = getWithRobotsTxt(`User-agent: *\nDisallow:${BRANCH_NAME === 'release' ? '' : ' *'}\n`)

export default function config(phase: string, defaultConfig: NextConfig) {
  const plugins = [
    withIgnore,
    withExtraWatch,
    withSvg,
    withFriendlyConsole,
    withMDX,
    withTypeChecking,
    PROFILE && withoutMinification,
    WATCH_POLL && withWebpackWatchPoll,
    withFriendlyChunkNames,
    withResolve,
    withRobotsTxt,
    withUrlAsset,
    PRODUCTION && withoutDebugPackage,
  ].filter(Boolean)

  return plugins.reduce(
    (acc, plugin) => {
      const update = plugin(acc)
      return typeof update === 'function' ? update(phase, defaultConfig) : update
    },
    { ...nextConfig },
  )
}
