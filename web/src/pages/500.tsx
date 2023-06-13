import React from 'react'
import { ErrorPage } from 'src/components/Error/ErrorPage'

export default function InternalServerErrorPage() {
  const error: Error & { statusCode?: number } = new Error('Internal server error')
  error.statusCode = 500
  return <ErrorPage statusCode={error.statusCode} title={error.message} error={error} showDetails={false} />
}
