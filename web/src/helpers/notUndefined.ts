/* eslint-disable lodash/prefer-is-nil */

export function notUndefined<T>(x: T | undefined): x is T {
  return x !== undefined
}

export function notUndefinedOrNull<T>(x: T | undefined): x is T {
  return x !== undefined && x !== null
}

export function maybe<T, Y>(f: (t: T) => Y, param?: T): Y | undefined {
  return param ? f(param) : undefined
}
