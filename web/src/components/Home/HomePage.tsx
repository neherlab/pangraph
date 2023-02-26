import React from 'react'
import styled from 'styled-components'
import { PageContainerHorizontal } from 'src/components/Layout/PageContainer'
import Viewer from 'src/components/Viewer/Viewer'

export function HomePage() {
  return (
    <PageContainerHorizontal>
      <MainContent>
        <MainContentInner>
          <Viewer />
        </MainContentInner>
      </MainContent>
    </PageContainerHorizontal>
  )
}

const MainContent = styled.div`
  display: flex;
  flex-direction: row;
  flex: 1 1 100%;
  overflow: hidden;
`

const MainContentInner = styled.div`
  display: flex;
  flex-direction: column;
  flex: 1 1 100%;
  overflow: hidden;
`
