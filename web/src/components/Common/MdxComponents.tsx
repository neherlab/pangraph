import styled, { css } from 'styled-components'
import type { MDXComponents } from 'mdx/types'
import { LinkSmart } from 'src/components/Common/Link'

export const CodeCommon = css`
  padding: 1px 5px;
  border-radius: 2px;
  background-color: ${(props) => props.theme.code.pre.background};
`

export const Pre = styled.pre`
  padding: 0.5rem 1rem;
  ${CodeCommon};
  overflow-x: auto;
  white-space: pre-wrap;
  word-wrap: break-word;
`

export const Code = styled.code`
  ${CodeCommon};
`

const P = styled.p`
  & > code {
    ${CodeCommon};
  }
`

export const mdxComponents = {
  a: LinkSmart,
  pre: Pre,
  code: Code,
  p: P,
}

export function getMdxComponents(components: MDXComponents): MDXComponents {
  return { ...components, ...mdxComponents }
}
