import { describe, expect, it } from 'vitest';
import handler from './serve-index-html-for-directories.js';

const testCases = [
  { uri: '', expected: '/index.html' },
  { uri: '/', expected: '/index.html' },
  { uri: '/foo', expected: '/foo/index.html' },
  { uri: '/foo/', expected: '/foo/index.html' },
  { uri: '/foo.html', expected: '/foo.html' },
  { uri: '/foo/bar', expected: '/foo/bar/index.html' },
  { uri: '/foo/bar/', expected: '/foo/bar/index.html' },
  { uri: '/foo/bar/baz.html', expected: '/foo/bar/baz.html' },
];

describe('CloudFront URL Rewriter Tests', () => {
  testCases.forEach(({ uri, expected }) => {
    it(`handles URI '${uri}' correctly`, async () => {
      const event = { request: { uri } };
      const result = await handler(event);
      expect(result.uri).toBe(expected);
    });
  });
});
