import React from 'react';
import { translate } from '@docusaurus/Translate';
import { PageMetadata } from '@docusaurus/theme-common';
import Layout from '@theme/Layout';
import NotFoundContent from '@theme/NotFound/Content';

/**
 * Customized 404 page.
 *
 * Swizzled (ejected) from the original theme using:
 *
 * ```bash
 * bun run swizzle --typescript --eject @docusaurus/theme-classic NotFound
 * ```
 *
 * */
export default function Index(): JSX.Element {
  const title = translate({
    id: 'theme.NotFound.title',
    message: 'Page Not Found',
  });
  return (
    <>
      <PageMetadata title={title} />
      <Layout>
        <NotFoundContent />
      </Layout>
    </>
  );
}
