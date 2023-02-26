import React from 'react'
import styled from 'styled-components'
import { Nav as NavBase, Navbar as NavbarBase } from 'reactstrap'
import { Link } from 'src/components/Common/Link'
import BrandLogoBase from 'src/assets/images/logo.svg'
import { LanguageSwitcher } from 'src/components/Layout/LanguageSwitcher'
import { PROJECT_NAME, THEME_COLOR_BG } from 'src/constants'

export const HEIGHT_NAVBAR = 50

export const Navbar = styled(NavbarBase)`
  height: ${HEIGHT_NAVBAR}px;
  width: 100%;
  box-shadow: ${(props) => props.theme.shadows.medium};
  z-index: 100;
  background-color: ${(props) => props.theme.white};
  opacity: 1;
`

export const Nav = styled(NavBase)`
  background-color: transparent !important;
`

export function NavigationBar() {
  return (
    <Navbar expand="xs" role="navigation">
      <Brand />

      <Nav className="ml-auto" navbar>
        <LanguageSwitcher />
      </Nav>
    </Navbar>
  )
}

const BrandWrapper = styled.div`
  display: flex;
  height: 100%;
  margin-bottom: 0.1rem;
`
export const BrandLogoSmall = styled(BrandLogoBase)`
  flex: 1;
  height: 36px;
  margin-left: 1rem;
`

const BrandNameStyled = styled.span`
  flex: 1;
  color: ${THEME_COLOR_BG};
  font-size: 24px;
  font-weight: bold;
  margin-left: 1rem;
  margin-right: 2rem;
`

export function Brand() {
  return (
    <Link href="/" className="text-decoration-none">
      <BrandWrapper>
        <BrandLogoSmall />
        <BrandNameStyled>{PROJECT_NAME}</BrandNameStyled>
      </BrandWrapper>
    </Link>
  )
}
