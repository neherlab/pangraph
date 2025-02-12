import React from 'react';
import clsx from 'clsx';
import Translate from '@docusaurus/Translate';
import type { Props } from '@theme/NotFound/Content';
import Heading from '@theme/Heading';
import Link from '@docusaurus/Link';

export default function NotFoundContent({ className }: Props): JSX.Element {
  return (
    <main className={clsx('container margin-vert--xl', className)}>
      <div className="row">
        <div className="col col--6 col--offset-3">
          <Heading as="h1" className="hero__title">
            <Translate id="theme.NotFound.title" description="The title of the 404 page">
              Page Not Found
            </Translate>
          </Heading>
          <p>
            <Translate id="theme.NotFound.p1" description="The first paragraph of the 404 page">
              We could not find what you were looking for.
            </Translate>
          </p>
          <p>
            <Translate>{"If you believe it's a bug, please report it at "}</Translate>
            <a href="https://github.com/neherlab/pangraph/issues" target="_blank" rel="noopener noreferrer">
              github.com/neherlab/pangraph/issues
            </a>
          </p>

          <Link to="/" className="button button--primary button--lg">
            <Translate>{'Go to main page'}</Translate>
          </Link>
        </div>
      </div>
    </main>
  );
}
