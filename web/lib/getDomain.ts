import { isNil } from 'lodash-es'
import { getenv } from './getenv'

const WEB_PORT_DEV = getenv('WEB_PORT_DEV', null)
const WEB_PORT_PROD = getenv('WEB_PORT_PROD', null)
const devDomain = `http://localhost:${WEB_PORT_DEV}`
const prodDomain = `http://localhost:${WEB_PORT_PROD}`

const ENV_VARS = [
  // prettier-ignore
  'VERCEL_URL',
  'NOW_URL',
  'ZEIT_URL',
  'DEPLOY_PRIME_URL',
  'DEPLOY_URL',
  'URL',
]

export function getenvFirst(vars: string[]): string | undefined {
  return vars.map((v) => getenv(v, null)).find((v) => !isNil(v))
}

function sanitizeDomain(domain: string) {
  return domain.startsWith('http') ? domain : `https://${domain}`
}

export function getDomain() {
  const domain = getenv('FULL_DOMAIN')
  if (domain === 'autodetect') {
    if (process.env.NODE_ENV === 'development') {
      return devDomain
    }
    return sanitizeDomain(getenvFirst(ENV_VARS) ?? prodDomain)
  }
  return sanitizeDomain(domain)
}
