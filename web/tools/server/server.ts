/**
 * Serves production build artifacts.
 *
 * /!\ Only for development purposes, e.g. verifying that production build runs
 * on developer's machine.
 *
 * This server is very naive, slow and insecure. Real-world deployments should
 * use either a 3rd-party static hosting or a robust static server, such as
 * Nginx, instead.
 *
 */

import type { ServerResponse } from 'http'
import path from 'path'
import express, { Request, Response } from 'express'
import allowMethods from 'allow-methods'
import compression from 'compression'
import morgan from 'morgan'
import expressStaticGzip from 'express-static-gzip'
import { getenv } from '../../lib/getenv'
import { findModuleRoot } from '../../lib/findModuleRoot'
import { modifyHeaders } from '../../infra/web/lambda-at-edge/ViewerResponse.lambda'

const { moduleRoot } = findModuleRoot()

const buildDir = path.join(moduleRoot, '.build', 'production', 'web')
const nextDir = path.join(buildDir, '_next')

export interface NewHeaders {
  [key: string]: { key: string; value: string }[]
}

const expressStaticGzipOptions = { enableBrotli: true, serveStatic: { extensions: ['html'] } }

const cacheNone = {
  ...expressStaticGzipOptions,
  serveStatic: {
    ...expressStaticGzipOptions.serveStatic,
    setHeaders: (res: ServerResponse) => res.setHeader('Cache-Control', 'no-cache'),
  },
}
const cacheOneYear = {
  ...expressStaticGzipOptions,
  serveStatic: {
    ...expressStaticGzipOptions.serveStatic,
    maxAge: '31556952000',
    immutable: true,
  },
}

function notFound(req: Request, res: Response) {
  res.status(404)
  res.format({
    json() {
      res.json({ error: 'Not found' })
    },
    default() {
      res.type('txt').send('Not found')
    },
  })
}

function setHeaders(req: express.Request, res: express.Response, next: express.NextFunction) {
  const newHeaders = modifyHeaders({ response: res }) as NewHeaders
  Object.entries(newHeaders).forEach(([header, arr]) => {
    const [{ value }] = arr
    if (header.toLowerCase() === 'strict-transport-security') {
      return
    }
    res.set({ [header.toLowerCase()]: value })
  })
  next()
}

function main() {
  const app = express()

  app.use(morgan('dev'))
  app.use(compression())
  app.use(allowMethods(['GET', 'HEAD']))
  app.use(setHeaders)
  app.use('/_next', expressStaticGzip(nextDir, cacheOneYear))
  app.get('*', expressStaticGzip(buildDir, cacheNone))
  app.use(notFound)

  const port = getenv('WEB_PORT_PROD')
  app.listen(port, () => {
    console.info(`Server is listening on port ${port}`)
  })
}

main()
