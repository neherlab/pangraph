import '../config/dotenv'

export default class EnvVarError extends TypeError {
  constructor(key, value) {
    super(`
      When reading an environment variable "${key}" (as \`process.env.${key}\`):
      it was expected to find a valid string, but found \`${value}\`.

      Have you followed the instructions in Developer's guide?

      There might have been additions to the list of environement variables required.
      Verify that your \`.env\` file has all the variables present in \`.env.example\`:

        diff --color .env.example .env

      In simple cases you might just copy the example:

        cp .env.example .env

      `)
  }
}

export function getenv(key, defaultValue) {
  const value = process.env[key]
  if (!value) {
    if (typeof defaultValue !== 'undefined') {
      return defaultValue
    }

    throw new EnvVarError(key, value)
  }
  return value
}

export function getbool(key, defaultValue) {
  const value = process.env[key]
  if (!value) {
    if (typeof defaultValue !== 'undefined') {
      return defaultValue
    }

    throw new EnvVarError(key, value)
  }

  return value === '1' || value === 'true' || value === 'yes'
}
