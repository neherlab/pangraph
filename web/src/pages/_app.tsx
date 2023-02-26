import 'reflect-metadata'
import dynamic from 'next/dynamic'
import React, { Suspense, useMemo } from 'react'
import { QueryClient, QueryClientConfig, QueryClientProvider } from '@tanstack/react-query'
import { RecoilRoot, RecoilEnv } from 'recoil'
import { ReactQueryDevtools } from '@tanstack/react-query-devtools'
import type { AppProps } from 'next/app'
import NextProgress from 'next-progress'
import { ThemeProvider } from 'styled-components'
import { MDXProvider } from '@mdx-js/react'
import { I18nextProvider } from 'react-i18next'
import i18n from 'src/i18n/i18n'
import { theme } from 'src/theme'
import { DOMAIN_STRIPPED } from 'src/constants'
import { Plausible } from 'src/components/Common/Plausible'
import { LOADING } from 'src/components/Layout/Loading'
import { getMdxComponents } from 'src/components/Common/MdxComponents'
import { ErrorBoundary } from 'src/components/Error/ErrorBoundary'
import { Layout } from 'src/components/Layout/Layout'
import 'src/styles/global.scss'

RecoilEnv.RECOIL_DUPLICATE_ATOM_KEY_CHECKING_ENABLED = false

const REACT_QUERY_OPTIONS: QueryClientConfig = {
  defaultOptions: {
    queries: {
      suspense: true,
      retry: 1,
      staleTime: Number.POSITIVE_INFINITY,
      refetchOnMount: false,
      refetchOnWindowFocus: false,
      refetchOnReconnect: true,
      refetchInterval: Number.POSITIVE_INFINITY,
    },
  },
}

const NEXT_PROGRESS_OPTIONS = { showSpinner: false }

function MyApp({ Component, pageProps }: AppProps) {
  const queryClient = useMemo(() => new QueryClient(REACT_QUERY_OPTIONS), [])

  return (
    <Suspense fallback={LOADING}>
      <ErrorBoundary>
        <QueryClientProvider client={queryClient}>
          <ThemeProvider theme={theme}>
            <I18nextProvider i18n={i18n}>
              <Plausible domain={DOMAIN_STRIPPED} />
              <NextProgress delay={100} options={NEXT_PROGRESS_OPTIONS} />
              <Suspense fallback={LOADING}>
                <RecoilRoot>
                  <MDXProvider components={getMdxComponents}>
                    <Layout>
                      <Suspense fallback={LOADING}>
                        <Component {...pageProps} />
                      </Suspense>
                    </Layout>
                  </MDXProvider>
                </RecoilRoot>
              </Suspense>
              <ReactQueryDevtools initialIsOpen={false} />
            </I18nextProvider>
          </ThemeProvider>
        </QueryClientProvider>
      </ErrorBoundary>
    </Suspense>
  )
}

export default dynamic(async () => MyApp, { ssr: false })
