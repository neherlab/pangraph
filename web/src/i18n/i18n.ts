import { ReactNode } from 'react'
import type { StrictOmit } from 'ts-essentials'
import { get, isNil, mapValues } from 'lodash-es'
import i18nOriginal, { i18n as I18N, Resource } from 'i18next'
import { initReactI18next } from 'react-i18next'
import { Settings as LuxonSettings } from 'luxon'
import numbro from 'numbro'
import { languages } from 'countries-list'
import prettyBytesOriginal, { Options as PrettyBytesOptionsOriginal } from 'pretty-bytes'
import { getFlag } from 'src/components/Common/CountryFlag'

// eslint-disable-next-line @typescript-eslint/ban-ts-comment
// @ts-ignore
import numbroLanguages from 'numbro/dist/languages.min'

import ar from './resources/ar/common.json'
import de from './resources/de/common.json'
import el from './resources/el/common.json'
import en from './resources/en/common.json'
import es from './resources/es/common.json'
import fa from './resources/fa/common.json'
import fr from './resources/fr/common.json'
import he from './resources/he/common.json'
import hi from './resources/hi/common.json'
import id from './resources/id/common.json'
import it from './resources/it/common.json'
import ja from './resources/ja/common.json'
import ko from './resources/ko/common.json'
import nl from './resources/nl/common.json'
import pt from './resources/pt/common.json'
import ru from './resources/ru/common.json'
import ta from './resources/ta/common.json'
import th from './resources/th/common.json'
import tr from './resources/tr/common.json'
import ur from './resources/ur/common.json'
import vi from './resources/vi/common.json'
import zh from './resources/zh/common.json'

export const translations = {
  ar,
  de,
  el,
  en,
  es,
  fa,
  fr,
  he,
  hi,
  id,
  it,
  ja,
  ko,
  nl,
  pt,
  ru,
  ta,
  th,
  tr,
  ur,
  vi,
  zh,
}

export type LocaleKey = keyof typeof translations

export const DEFAULT_LOCALE_KEY: LocaleKey = 'en'
export const resources: Record<LocaleKey, Resource> = mapValues(translations, (value) => ({ translation: value }))

export interface Locale {
  readonly key: LocaleKey
  readonly full: string
  readonly name: string
  readonly native: string
  readonly flag: ReactNode
}

// prettier-ignore
export const locales: Record<LocaleKey, Locale> = {
  en: {
    key: 'en',
    full: 'en-US',
    name: languages.en.name,
    native: languages.en.native,
    flag: getFlag('United Kingdom'),
  },
  ar: { key: 'ar', full: 'ar-SA', name: languages.ar.name, native: languages.ar.native, flag: getFlag('Saudi Arabia') },
  de: { key: 'de', full: 'de-DE', name: languages.de.name, native: languages.de.native, flag: getFlag('Germany') },
  el: { key: 'el', full: 'el-GR', name: languages.el.name, native: languages.el.native, flag: getFlag('Greece') },
  es: { key: 'es', full: 'es-ES', name: languages.es.name, native: languages.es.native, flag: getFlag('Spain') },
  fa: { key: 'fa', full: 'fa-IR', name: languages.fa.name, native: languages.fa.native, flag: getFlag('Iran') },
  fr: { key: 'fr', full: 'fr-FR', name: languages.fr.name, native: languages.fr.native, flag: getFlag('France') },
  he: { key: 'he', full: 'he-IL', name: languages.he.name, native: languages.he.native, flag: getFlag('Israel') },
  hi: { key: 'hi', full: 'hi-IN', name: languages.hi.name, native: languages.hi.native, flag: getFlag('India') },
  id: { key: 'id', full: 'id-ID', name: languages.id.name, native: languages.id.native, flag: getFlag('Indonesia') },
  it: { key: 'it', full: 'it-IT', name: languages.it.name, native: languages.it.native, flag: getFlag('Italy') },
  ja: { key: 'ja', full: 'ja-JP', name: languages.ja.name, native: languages.ja.native, flag: getFlag('Japan') },
  ko: { key: 'ko', full: 'ko-KR', name: languages.ko.name, native: languages.ko.native, flag: getFlag('South Korea') },
  nl: { key: 'nl', full: 'nl-NL', name: languages.nl.name, native: languages.nl.native, flag: getFlag('Netherlands') },
  pt: { key: 'pt', full: 'pt-PT', name: languages.pt.name, native: languages.pt.native, flag: getFlag('Portugal') },
  ru: { key: 'ru', full: 'ru-RU', name: languages.ru.name, native: languages.ru.native, flag: getFlag('Russia') },
  ta: { key: 'ta', full: 'ta-IN', name: languages.ta.name, native: languages.ta.native, flag: getFlag('India') },
  th: { key: 'th', full: 'th-TH', name: languages.th.name, native: languages.th.native, flag: getFlag('Thailand') },
  tr: { key: 'tr', full: 'tr-TR', name: languages.tr.name, native: languages.tr.native, flag: getFlag('Turkey') },
  ur: { key: 'ur', full: 'ur-PK', name: languages.ur.name, native: languages.ur.native, flag: getFlag('Pakistan') },
  vi: { key: 'vi', full: 'vi-VN', name: languages.vi.name, native: languages.vi.native, flag: getFlag('Vietnam') },
  zh: { key: 'zh', full: 'zh-CN', name: languages.zh.name, native: languages.zh.native, flag: getFlag('China') },
} as const

export const localeKeys = Object.keys(locales)

export const localesArray: Locale[] = Object.values(locales)

export interface I18NInitParams {
  localeKey: LocaleKey
}

export type PrettyBytesOptions = StrictOmit<PrettyBytesOptionsOriginal, 'locale'>

export class PrettyBytes {
  private localeKey: string = DEFAULT_LOCALE_KEY as string

  public setLocale(localeKey: string) {
    this.localeKey = getLocaleWithKey(localeKey).key
  }

  public format(numBytes: number, options?: PrettyBytesOptions) {
    return prettyBytesOriginal(numBytes, { binary: true, ...options, locale: this.localeKey })
  }
}

const prettyBytes = new PrettyBytes()

export function i18nInit({ localeKey }: I18NInitParams) {
  const enUS = numbro.languages()['en-US']
  const allNumbroLanguages = numbroLanguages as Record<string, numbro.NumbroLanguage>

  Object.values(locales).forEach((locale) => {
    const numbroLocale = Object.values(allNumbroLanguages).find(
      (numbroLocale) => numbroLocale.languageTag === locale.full,
    )
    numbro.registerLanguage({ ...enUS, ...numbroLocale })
  })

  const i18n = i18nOriginal.use(initReactI18next).createInstance({
    resources,
    lng: localeKey,
    fallbackLng: DEFAULT_LOCALE_KEY,
    debug: process.env.DEV_ENABLE_I18N_DEBUG === '1',
    keySeparator: false, // Disable dots as key separators as we use dots in keys
    nsSeparator: false,
    interpolation: { escapeValue: false },
    initImmediate: true,
  })

  // eslint-disable-next-line no-void
  void i18n.init()

  const locale = locales[localeKey]
  LuxonSettings.defaultLocale = localeKey
  numbro.setLanguage(locale.full)
  numbro.setDefaults({ thousandSeparated: true })

  return i18n
}

export function getLocaleWithKey(key: string) {
  const locale = get(locales, key) as Locale
  if (isNil(locale)) {
    return { ...locales[DEFAULT_LOCALE_KEY], key: DEFAULT_LOCALE_KEY }
  }
  return locale
}

export async function changeLocale(i18n: I18N, localeKey: string) {
  if (localeKeys.includes(localeKey)) {
    const locale = getLocaleWithKey(localeKey)
    LuxonSettings.defaultLocale = localeKey
    numbro.setLanguage(locale.full)
    await i18n.changeLanguage(localeKey)
    prettyBytes.setLocale(localeKey)
    return true
  }
  return false
}

const i18n = i18nInit({ localeKey: DEFAULT_LOCALE_KEY })

export { prettyBytes }

export default i18n

export { default as numbro } from 'numbro'
