import React from 'react'

import { useTranslationSafe as useTranslation } from 'src/helpers/useTranslationSafe'
import { Col, Container, Row } from 'reactstrap'
import styled from 'styled-components'

import { PROJECT_NAME, COMPANY_NAME } from 'src/constants'
import { getCopyrightYearRange } from 'src/helpers/getCopyrightYearRange'
import { LinkExternal } from 'src/components/Common/Link'
import { getVersionString } from 'src/helpers/getVersionString'

import LogoBiozentrum from 'src/assets/images/biozentrum_square.svg'
import LogoNeherlab from 'src/assets/images/neherlab.svg'
import LogoNextjs from 'src/assets/images/nextjs.svg'
import LogoSib from 'src/assets/images/sib.svg'
import LogoUnibas from 'src/assets/images/unibas.svg'
import LogoVercel from 'src/assets/images/vercel.svg'

export const PAGE_FOOTER_HEIGHT = 36

const FooterContainer = styled(Container)`
  position: fixed;
  height: ${PAGE_FOOTER_HEIGHT}px;
  bottom: 0;
  padding: 6px 10px;
  box-shadow: ${(props) => props.theme.shadows.medium};
  z-index: 100;
  background-color: ${(props) => props.theme.white};
  opacity: 1;
`

const CopyrightText = styled.div`
  font-size: 0.75rem;
  flex-grow: 1;

  @media (max-width: 576px) {
    font-size: 0.5rem;
  }
`

const LogoContainer = styled.div`
  flex-grow: 1;
  text-align: right;
`

const VersionText = styled.div`
  flex-grow: 1;
  font-size: 0.75rem;
  text-align: right;

  @media (max-width: 992px) {
    display: none;
  }
`

const LogoLink = styled(LinkExternal)`
  padding: 10px 20px;

  svg {
    height: 22px;
  }

  @media (max-width: 768px) {
    padding: 5px 4px;
    flex-grow: 1;

    svg {
      height: 20px;
    }
  }

  @media (max-width: 576px) {
    padding: 5px 5px;
    flex-grow: 1;

    svg {
      height: 20px;
    }
  }
`

export function Footer() {
  const { t } = useTranslation()
  const copyrightYearRange = getCopyrightYearRange()

  return (
    <FooterContainer fluid tag="footer">
      <Row noGutters>
        <Col className="d-flex">
          <CopyrightText className="mr-auto my-auto">
            {t('{{PROJECT_NAME}} (c) {{copyrightYearRange}} {{COMPANY_NAME}}', {
              PROJECT_NAME,
              copyrightYearRange,
              COMPANY_NAME,
            })}
          </CopyrightText>

          <LogoContainer className="mx-auto">
            <LogoLink icon={null} href="https://neherlab.org" title="NeherLab">
              <LogoNeherlab />
            </LogoLink>

            <LogoLink icon={null} href="https://www.biozentrum.unibas.ch" title="Biozentrum">
              <LogoBiozentrum />
            </LogoLink>

            <LogoLink icon={null} href="https://www.unibas.ch" title="University of Basel">
              <LogoUnibas />
            </LogoLink>

            <LogoLink icon={null} href="https://www.sib.swiss" title="Swiss institute of Bioinformatics">
              <LogoSib />
            </LogoLink>

            <LogoLink icon={null} href={`https://https://nextjs.org/?utm_source=${COMPANY_NAME}`} title="Next.js">
              <LogoNextjs />
            </LogoLink>

            <LogoLink icon={null} href={`https://vercel.com/?utm_source=${COMPANY_NAME}`} title="Vercel">
              <LogoVercel />
            </LogoLink>
          </LogoContainer>

          <VersionText className="ml-auto my-auto">{getVersionString()}</VersionText>
        </Col>
      </Row>
    </FooterContainer>
  )
}
