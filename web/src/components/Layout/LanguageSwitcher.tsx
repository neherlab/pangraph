import React, { useCallback, useMemo } from 'react'
import { useRecoilState } from 'recoil'
import styled from 'styled-components'
import { Locale, localesArray } from 'src/i18n/i18n'
import { localeAtom } from 'src/state/locale.state'
import { ErrorInternal } from 'src/helpers/ErrorInternal'
import { Dropdown as DropdownBase, DropdownEntry, DropdownProps } from 'src/components/Common/Dropdown'
import { FlagWrapper as FlagWrapperBase } from 'src/components/Common/CountryFlag'

export function LanguageSwitcher({
  ...restProps
}: Omit<DropdownProps, 'entries' | 'currentEntry' | 'setCurrentEntry'>) {
  const [currentLocale, setCurrentLocale] = useRecoilState(localeAtom)
  const setCurrentEntry = useCallback((entry: DropdownEntry) => setCurrentLocale(entry.key), [setCurrentLocale])

  const { entries } = useMemo(() => {
    const entries = localesArray.map((locale) => ({
      key: locale.key,
      value: <LanguageSwitcherMenuItem locale={locale} />,
    }))
    return { entries }
  }, [])

  const currentEntry = useMemo(() => {
    const locale = localesArray.find((entry) => entry.key === currentLocale)
    if (!locale) {
      throw new ErrorInternal(`Locale not found: '${currentLocale}'`)
    }
    return { key: locale.key, value: <LanguageSwitcherCurrentItem locale={locale} /> }
  }, [currentLocale])

  return <Dropdown {...restProps} entries={entries} currentEntry={currentEntry} setCurrentEntry={setCurrentEntry} />
}

export function LanguageSwitcherMenuItem({ locale }: { locale: Locale }) {
  const { flag, name, native, key } = locale

  const label = useMemo(() => {
    if (name === native) {
      return `${name}`
    }

    return `${name} (${native})`
  }, [name, native])

  return (
    <LanguageSwitcherItemWrapper>
      <FlagWrapper>{flag}</FlagWrapper>
      <span className="pl-2 text-monospace">{`${key.toUpperCase()}`}</span>
      <span className="pl-2">{label}</span>
    </LanguageSwitcherItemWrapper>
  )
}

export function LanguageSwitcherCurrentItem({ locale }: { locale: Locale }) {
  const { flag, key } = locale
  return (
    <LanguageSwitcherItemWrapper>
      <FlagWrapper>{flag}</FlagWrapper>
      <span className="pl-2 text-monospace">{`${key.toUpperCase()}`}</span>
    </LanguageSwitcherItemWrapper>
  )
}

const FlagWrapper = styled(FlagWrapperBase)`
  height: 1.4rem;
  width: 2rem;
`

const LanguageSwitcherItemWrapper = styled.span`
  display: flex;
`

const Dropdown = styled(DropdownBase)`
  border: none !important;
`
