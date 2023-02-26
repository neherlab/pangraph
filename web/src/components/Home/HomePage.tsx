import React from 'react'
import styled from 'styled-components'
import { PROJECT_NAME } from 'src/constants'
import { PageContainerHorizontal } from 'src/components/Layout/PageContainer'

export function HomePage() {
  return (
    <PageContainerHorizontal>
      <MainContent>
        <MainContentInner>
          <h2>{PROJECT_NAME}</h2>
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
