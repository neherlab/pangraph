import React, { ReactNode, useCallback, useMemo, useState } from 'react'
import styled from 'styled-components'
import { rgba, darken } from 'polished'
import {
  Dropdown as DropdownBase,
  DropdownToggle as DropdownToggleBase,
  DropdownMenu as DropdownMenuBase,
  DropdownItem as DropdownItemBase,
  DropdownProps as DropdownBaseProps,
} from 'reactstrap'
import { MdArrowDropDown } from 'react-icons/md'

export interface DropdownEntry {
  key: string
  value: ReactNode
}

export interface DropdownProps extends DropdownBaseProps {
  entries: DropdownEntry[]
  currentEntry: DropdownEntry
  setCurrentEntry(key: DropdownEntry): void
}

export function Dropdown({ currentEntry, setCurrentEntry, entries, ...restProps }: DropdownProps) {
  const [dropdownOpen, setDropdownOpen] = useState(false)
  const toggle = useCallback(() => setDropdownOpen((prevState) => !prevState), [])
  const onClick = useCallback((entry: DropdownEntry) => () => setCurrentEntry(entry), [setCurrentEntry])
  const menuItems = useMemo(() => {
    return entries.map((entry) => {
      return (
        <DropdownItem key={entry.key} active={entry.key === currentEntry.key} onClick={onClick(entry)}>
          {entry.value}
        </DropdownItem>
      )
    })
  }, [currentEntry, entries, onClick])

  return (
    <DropdownStyled isOpen={dropdownOpen} toggle={toggle} {...restProps}>
      <DropdownToggle nav>
        <DropdownToggleText>{currentEntry.value}</DropdownToggleText>
        <DropdownCaret />
      </DropdownToggle>
      <DropdownMenu positionFixed>{menuItems}</DropdownMenu>
    </DropdownStyled>
  )
}

const DropdownStyled = styled(DropdownBase)`
  display: flex;
  border: ${(props) => rgba(props.theme.primary, 0.5)} solid 1px;
  border-radius: 3px;
  box-shadow: ${(props) => props.theme.shadows.lighter};

  color: ${(props) => props.theme.bodyColor};

  :hover {
    color: ${(props) => props.theme.bodyColor};
  }

  & > a {
    // min-width: 200px;
  }
`

const DropdownToggle = styled(DropdownToggleBase)`
  display: flex;
  color: ${(props) => props.theme.bodyColor};
`

const DropdownToggleText = styled.span`
  flex: 1;
  color: ${(props) => props.theme.bodyColor};

  :hover {
    color: ${(props) => props.theme.bodyColor};
  }
`

const DropdownCaret = styled(MdArrowDropDown)`
  margin: 0 auto;
  width: 22px;
  height: 22px;
`

const DropdownMenu = styled(DropdownMenuBase)`
  color: ${(props) => props.theme.bodyColor};
  background-color: ${(props) => props.theme.bodyBg};
  box-shadow: ${(props) => props.theme.shadows.blurredMedium};
  max-height: 70vh;
  overflow-y: scroll;
`

const DropdownItem = styled(DropdownItemBase)<{ active: boolean }>`
  :hover {
    color: ${({ active, theme }) => (active ? theme.white : theme.bodyColor)};
    background: ${({ active, theme }) => (active ? darken(0.025)(theme.primary) : theme.gray200)};
  }
`
