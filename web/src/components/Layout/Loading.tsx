import React, { useEffect, useState } from 'react'
import styled from 'styled-components'

import { useTranslationSafe } from 'src/helpers/useTranslationSafe'
import Logo from 'src/assets/images/logo.svg'

const Container = styled.div`
  display: flex;
  width: 100%;
  height: 100%;
  overflow: hidden;
`

const SpinningLogo = styled(Logo)`
  margin: auto;
  width: 100px;
  height: 100px;

  box-shadow: 0 0 0 0 rgba(0, 0, 0, 1);
  transform: scale(1);
  animation: pulse 2s ease-out infinite;

  @keyframes pulse {
    0% {
      transform: scale(0.66);
    }
    70% {
      transform: scale(1);
    }
    100% {
      transform: scale(0.66);
    }
  }
`

function Loading() {
  const { t } = useTranslationSafe()

  const [show, setShow] = useState(false)

  useEffect(() => {
    const timer = setTimeout(() => setShow(true), 300)
    return () => {
      clearTimeout(timer)
    }
  })

  if (!show) {
    return null
  }

  return (
    <Container title={t('Loading...')}>
      <SpinningLogo />
    </Container>
  )
}

export default Loading
export const LOADING = <Loading />
