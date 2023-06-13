/* eslint-disable no-template-curly-in-string,sonarjs/no-duplicate-string */

/**
 *
 * This dangerously and unreliably patches some node_modules. Mostly cosmetic stuff.
 * Do not use this to fix bugs or introduce features. Consider contributing to the upstream project instead.
 *
 */
import { concurrent } from 'fasy'
import fs from 'fs-extra'
import { isArray } from 'lodash-es'

export async function replace(filename: string | string[], searchValue: string | RegExp, replaceValue = '') {
  const filenames = !isArray(filename) ? [filename] : filename
  await concurrent.forEach(async (filename) => {
    const content = await fs.readFile(filename, 'utf8')
    const newContent = content.replace(searchValue, replaceValue)
    await fs.writeFile(filename, newContent, { encoding: 'utf8' })
  }, filenames)
}

export async function main() {
  // Removes reminder about upgrading caniuse database. Nice, but not that important. Will be handled along with
  // routine package updates.
  // Reason: too noisy
  await replace(
    ['node_modules/browserslist/node.js'],
    `      console.warn(
        'Browserslist: caniuse-lite is outdated. Please run:\\n' +
          '  npx browserslist@latest --update-db\\n' +
          '  Why you should do it regularly: ' +
          'https://github.com/browserslist/browserslist#browsers-data-updating'
      )`,
  )

  await replace(
    [
      'node_modules/next/dist/compiled/browserslist/index.js',
      'node_modules/next/dist/compiled/cssnano-simple/index.js',
    ],
    `console.warn("Browserslist: caniuse-lite is outdated. Please run:\\n"+"  npx browserslist@latest --update-db\\n"+"  Why you should do it regularly: "+"https://github.com/browserslist/browserslist#browsers-data-updating")`,
  )

  // Fast refresh messages in browser console
  await replace(
    'node_modules/next/dist/client/dev/error-overlay/hot-dev-client.js',
    "console.log('[Fast Refresh] rebuilding');",
  )
  await replace(
    'node_modules/next/dist/client/dev/error-overlay/hot-dev-client.js',
    'console.log(`[Fast Refresh] done in ${latency}ms`);',
  )

  await replace(
    'node_modules/next/dist/server/base-server.js',
    'Log.warn(`You have added a custom /_error page without a custom /404 page. This prevents the 404 page from being auto statically optimized.\\nSee here for info: https://nextjs.org/docs/messages/custom-error-no-custom-404`);',
  )

  // Removes warning "<title> should not be used in _document.js".
  // Reason: We want title and other SEO tags to be pre-rendered, so that crawlers could find them.
  await replace(
    ['node_modules/next/dist/pages/_document.js'],
    `console.warn("Warning: <title> should not be used in _document.js's <Head>. https://nextjs.org/docs/messages/no-document-title");`,
  )

  await replace(
    ['node_modules/next/dist/build/index.js'],
    "`${Log.prefixes.info} ${ignoreTypeScriptErrors ? 'Skipping validation of types' : 'Checking validity of types'}`",
    '""',
  )

  // More useless messages from Next.js
  await replace(['node_modules/next/dist/server/config.js'], 'console.warn();')

  await replace(
    ['node_modules/next/dist/server/config.js'],
    "Log.warn('SWC minify release candidate enabled. https://nextjs.org/docs/messages/swc-minify-enabled');",
  )

  await replace(
    [
      'node_modules/next/dist/compiled/next-server/next-server.js',
      'node_modules/next/dist/server/config.js',
      'node_modules/next/dist/esm/server/config.js',
    ],
    'Log.warn(_chalk.default.bold(`You have enabled experimental feature${s} (${features.join(", ")}) in ${configFileName}.`));',
  )

  await replace(
    [
      'node_modules/next/dist/compiled/next-server/next-server.js',
      'node_modules/next/dist/server/config.js',
      'node_modules/next/dist/esm/server/config.js',
    ],
    'Log.warn(chalk.bold(`You have enabled experimental feature${s} (${features.join(", ")}) in ${configFileName}.`));',
  )

  await replace(
    [
      'node_modules/next/dist/compiled/next-server/next-server.js',
      'node_modules/next/dist/esm/server/config.js',
      'node_modules/next/dist/server/config.js',
    ],
    'Wt.warn(`Experimental features are not covered by semver, and may cause unexpected or broken application behavior. `+`Use at your own risk.`);',
  )

  await replace(
    [
      'node_modules/next/dist/compiled/next-server/next-server.js',
      'node_modules/next/dist/esm/server/config.js',
      'node_modules/next/dist/server/config.js',
    ],
    'Log.warn(`Experimental features are not covered by semver, and may cause unexpected or broken application behavior. ` + `Use at your own risk.`);',
  )

  await replace(
    ['node_modules/next/dist/build/webpack-config.js', 'node_modules/next/dist/esm/build/webpack-config.js'],
    'Log.info(`automatically enabled Fast Refresh for ${injections} custom loader${injections > 1 ? "s" : ""}`);',
  )

  await replace(
    ['node_modules/@next/env/dist/index.js', 'node_modules/next/dist/compiled/next-server/next-server.js'],
    'n.info(`Loaded env from ${t.join(r||"",o.path)}`)',
  )

  await replace(
    ['node_modules/next/dist/build/output/store.js', 'node_modules/next/dist/esm/build/output/store.js'],
    'Log.wait("compiling...");',
  )

  await replace(
    ['node_modules/next/dist/build/output/store.js', 'node_modules/next/dist/esm/build/output/store.js'],
    'Log.wait(`compiling ${state.trigger}...`);',
  )

  await replace(
    'node_modules/next/dist/build/output/store.js',
    'Log.info(`bundled${partialMessage} successfully${timeMessage}${modulesMessage}, waiting for typecheck results...`);',
  )

  await replace(
    'node_modules/next/dist/build/output/store.js',
    'Log.event(`compiled${partialMessage} successfully${timeMessage}${modulesMessage}`);',
  )

  await replace(
    ['node_modules/next/dist/server/image-optimizer.js', 'node_modules/next/dist/esm/server/image-optimizer.js'],
    'console.warn(chalk.yellow.bold("Warning: ") + `For production Image Optimization with Next.js, the optional \'sharp\' package is strongly recommended. Run \'yarn add sharp\', and Next.js will use it automatically for Image Optimization.\\n` + "Read more: https://nextjs.org/docs/messages/sharp-missing-in-production");',
  )

  await replace(
    ['node_modules/next/dist/server/image-optimizer.js', 'node_modules/next/dist/esm/server/image-optimizer.js'],
    'console.warn(_chalk.default.yellow.bold("Warning: ") + `For production Image Optimization with Next.js, the optional \'sharp\' package is strongly recommended. Run \'yarn add sharp\', and Next.js will use it automatically for Image Optimization.\\n` + "Read more: https://nextjs.org/docs/messages/sharp-missing-in-production");',
  )

  // Fix timestamp in @nuxt/friendly-errors-webpack-plugin
  await replace(
    'node_modules/@nuxt/friendly-errors-webpack-plugin/src/reporters/base.js',
    "return `${message}${' '.repeat(logSpace)}${dateString}`",
    "return `${message}${' '.repeat(logSpace)}`",
  )

  // Copy next.js cli such that it has file extension (otherwise ts-node + Node 18 throw an error)
  await fs.copy('node_modules/next/dist/bin/next', 'node_modules/next/dist/bin/next.js', {
    dereference: true,
    overwrite: true,
  })
}

await main()
