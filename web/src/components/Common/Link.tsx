import React, { HTMLAttributes, HTMLProps, PropsWithChildren, useMemo } from 'react'
import styled from 'styled-components'
import { isNil, noop } from 'lodash-es'
import isAbsoluteUrl from 'is-absolute-url'
import NextLink, { LinkProps as NextLinkProps } from 'next/link'
import { StrictOmit } from 'ts-essentials'
import { Button, ButtonProps } from 'reactstrap'
import { GoLinkExternal } from 'react-icons/go'

export interface LinkProps
  extends StrictOmit<PropsWithChildren<NextLinkProps & HTMLAttributes<HTMLAnchorElement>>, 'href'> {
  href?: string
  className?: string
}

export function Link({ className, children, href, onClick, ...restProps }: LinkProps) {
  if (isNil(href)) {
    return (
      // eslint-disable-next-line jsx-a11y/click-events-have-key-events,jsx-a11y/no-static-element-interactions
      <a {...restProps} className={className} onClick={noop}>
        {children}
      </a>
    )
  }

  return (
    <NextLink {...restProps} href={href} passHref={false} className={className} onClick={onClick}>
      {children}
    </NextLink>
  )
}

export interface LinkLookingLikeButtonProps extends PropsWithChildren<ButtonProps> {
  href?: string
}

export function LinkLookingLikeButton({ href, children, ...restProps }: LinkLookingLikeButtonProps) {
  return (
    <Link href={href} passHref legacyBehavior>
      <Button {...restProps}>{children}</Button>
    </Link>
  )
}

const LinkExternalIconWrapper = styled.span<{ $color?: string }>`
  color: ${(props) => props.$color ?? props.theme.link.dim.iconColor};
`

export interface LinkExternalProps extends Omit<HTMLProps<HTMLAnchorElement>, 'as' | 'ref'> {
  href?: string
  color?: string
  iconColor?: string
  icon?: React.ReactNode
}

const A = styled.a<{ $color?: string } & LinkExternalProps>`
  color: ${(props) => props.$color ?? undefined};
  text-decoration: none;

  &:hover {
    color: ${(props) => props.$color ?? undefined};
    text-decoration: none;
  }

  white-space: nowrap;
`

export const ContentWrapper = styled.span`
  white-space: normal;
`

export function LinkExternal({
  href,
  color,
  iconColor,
  icon,
  children,
  ...restProps
}: PropsWithChildren<LinkExternalProps>) {
  const Icon: React.ReactNode = icon === undefined ? <GoLinkExternal /> : icon

  return (
    <A target="_blank" rel="noopener noreferrer" href={href} $color={color} {...restProps}>
      <ContentWrapper>{children}</ContentWrapper>{' '}
      {Icon && <LinkExternalIconWrapper $color={iconColor}>{Icon}</LinkExternalIconWrapper>}
    </A>
  )
}

export interface LinkSmartProps extends StrictOmit<LinkProps, 'href' | 'as'> {
  href?: string
}

export function LinkSmart({ href, ...restProps }: LinkSmartProps & LinkExternalProps) {
  const external = useMemo(() => isAbsoluteUrl(href ?? ''), [href])

  if (!href) {
    return <span {...restProps} />
  }

  if (external) {
    return <LinkExternal href={href} {...restProps} />
  }

  return <Link href={href} {...restProps} />
}
