import { themes as prismThemes } from 'prism-react-renderer';
import type { Config } from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

const config: Config = {
  title: 'Pangraph docs',
  tagline: 'Bioinformatic toolkit to align genome assemblies into pangenome graphs',
  favicon: 'img/favicon.ico',

  url: 'https://pangraph.github.io', // FIXME
  baseUrl: '/',

  organizationName: 'neherlab',
  projectName: 'pangraph',

  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      {
        docs: {
          routeBasePath: '/',
          sidebarPath: './sidebars.ts',
          editUrl: 'https://github.com/neherlab/pangraph/tree/master/docs',
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
        },
        blog: false,
        theme: {
          customCss: './src/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    image: 'img/social.png',
    navbar: {
      title: 'Pangraph',
      logo: {
        alt: 'Pangraph logo',
        src: 'img/logo.svg',
      },
      items: [
        {
          href: 'https://github.com/neherlab/pangraph',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    prism: {
      theme: prismThemes.vsLight,
      darkTheme: prismThemes.nightOwl,
      additionalLanguages: ['bash', 'powershell'],
    },
  } satisfies Preset.ThemeConfig,
} satisfies Config;

export default config;
