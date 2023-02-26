import { DateTime, DurationLike } from 'luxon'
import { getLocaleWithKey, numbro } from 'src/i18n/i18n'

export const formatProportion = (locale: string) => (value: number) => {
  try {
    return new Intl.NumberFormat(getLocaleWithKey(locale).full, {
      maximumFractionDigits: 2,
      minimumFractionDigits: 2,
    }).format(value)
  } catch {
    return numbro(value).format({ totalLength: 3 })
  }
}

export const formatRange = (locale: string) => (begin: number, end: number) => {
  return `[${formatProportion(locale)(begin)}; ${formatProportion(locale)(end)}]`
}

export const formatInteger = (locale: string) => (value: number) => {
  void locale // eslint-disable-line no-void
  try {
    return new Intl.NumberFormat(getLocaleWithKey(locale).full).format(value)
  } catch {
    return numbro(value).format()
  }
}

export function timestampToDate(seconds: number) {
  return DateTime.fromSeconds(seconds)
}

export function formatDateRange(weekTimestamp: number, range: DurationLike) {
  const begin = DateTime.fromSeconds(weekTimestamp)
  const end = begin.plus(range)
  return `${begin.toFormat('dd MMM yyyy')} - ${end.toFormat('dd MMM yyyy')}`
}

export function formatDateWeekly(weekTimestamp: number) {
  return formatDateRange(weekTimestamp, { weeks: 1 })
}

export function formatDateBiweekly(weekTimestamp: number) {
  return formatDateRange(weekTimestamp, { weeks: 2 })
}

export const formatDateHumanely = (locale: string) => (date: number) => {
  void locale // eslint-disable-line no-void
  return DateTime.fromSeconds(date).toLocaleString({ month: 'short', year: 'numeric' }).replace(' ', '\n')
}

export function dateFromYmd(ymd: string): DateTime {
  return DateTime.fromFormat(ymd, 'yyyy-MM-dd')
}

export function ymdToTimestamp(ymd: string): number {
  return dateFromYmd(ymd).toSeconds()
}
