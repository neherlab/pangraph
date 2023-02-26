import React, { ReactElement } from 'react'
import styled from 'styled-components'
import iso3311a2 from 'iso-3166-1-alpha-2'
import Flags from 'country-flag-icons/react/3x2'

export const FlagWrapper = styled.div`
  height: calc(1em + 2px);
  width: calc(1.33em + 2px);
  border: 0.5px solid #aaaa;
  display: flex;

  > * {
    width: 100%;
    height: 100%;
  }
`

export const missingCountryCodes: Record<string, string> = {
  'Bolivia': 'BO',
  'Bonaire': 'BQ',
  'Brunei': 'BN',
  'Cabo Verde': 'CV',
  'Curacao': 'CW',
  'Democratic Republic of the Congo': 'CD',
  'Eswatini': 'SZ',
  'Iran': 'IR',
  'Kosovo': 'XK',
  'Laos': 'LA',
  'Moldova': 'MD',
  'North Macedonia': 'MK',
  'Republic of the Congo': 'CD',
  'Russia': 'RU',
  'Saint Martin': 'SX',
  'Sint Maarten': 'SX',
  'South Korea': 'KR',
  'Taiwan': 'TW',
  'USA': 'US',
  'Venezuela': 'VE',
  'Vietnam': 'VN',
}

export function getFlag(country: string): ReactElement | null {
  const countryCode = missingCountryCodes[country] ?? iso3311a2.getCode(country) ?? '?'

  const Flag = Flags[countryCode]
  if (Flag) {
    return <Flag />
  }

  return null
}
