import React, { PropsWithChildren } from 'react'
import styled from 'styled-components'
import { HEIGHT_NAVBAR, NavigationBar } from 'src/components/Layout/NavigationBar'
import { Footer } from 'src/components/Layout/Footer'

const Header = styled.div`
  height: ${HEIGHT_NAVBAR}px;
`

const Body = styled.div`
  position: absolute;
  top: calc(${HEIGHT_NAVBAR}px);
  right: 0;
  bottom: 0;
  left: 0;
  display: flex;
  background-color: ${(props) => props.theme.bodyBg};
`

export interface LayoutProps {
  wide?: boolean
}

export function Layout({ children }: PropsWithChildren<LayoutProps>) {
  return (
    <div>
      <Header>
        <NavigationBar />
      </Header>

      <Body>{children}</Body>

      <Footer />
    </div>
  )
}
