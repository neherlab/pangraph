import React from 'react'
import styled from 'styled-components'
import Scene from 'src/components/Viewer/Scene'

export const Container = styled.div`
  margin: 0;
  padding: 0;
  height: 100%;
  overflow: auto;
`

export const ViewerWrapper = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  width: 100%;
  height: 100%;
`

export const Window = styled.div`
  flex: 1 0 50%;
  height: 100%;
  overflow: auto;
`

export default function Viewer() {
  return (
    <Container>
      <ViewerWrapper>
        <Window>
          <Scene name="1" color="#000000" />
        </Window>
      </ViewerWrapper>
    </Container>
  )
}
