import { Container } from 'reactstrap'
import styled from 'styled-components'

export const PageContainerNarrow = styled(Container)`
  display: flex;
  flex-direction: column;
  max-width: 1200px;
  margin: 0 auto;
`

export const PageContainer = styled(Container)`
  display: flex;
  flex-direction: column;
  margin: 0 auto;
`

export const PageContainerHorizontal = styled.div`
  display: flex;
  flex: 1;
  flex-direction: row;
`

export const PageMainWrapper = styled.main`
  flex: 1;
  max-height: 100%;
  overflow: hidden auto;
`
