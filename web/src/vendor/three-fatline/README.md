three-fatline
=============

[![NPM package][npm-img]][npm-url]
[![Build Size][build-size-img]][build-size-url]
[![NPM Downloads][npm-downloads-img]][npm-downloads-url]

A modularized version of [https://github.com/mrdoob/three.js/blob/master/examples/js/lines/](https://github.com/mrdoob/three.js/blob/master/examples/js/lines/).

See the [example](https://vasturiano.github.io/three-fatline/example/fat-lines/) ([source](https://github.com/vasturiano/three-fatline/blob/master/example/fat-lines/index.html))

#### Usage Example:
```js
import { Line2, LineGeometry, LineMaterial } from 'three-fatline';

const geometry = new LineGeometry();
geometry.setPositions([-160, 200, 0, 0, -200, 0, 160, 200, 0]); // [ x1, y1, z1,  x2, y2, z2, ... ] format

const material = new LineMaterial({
  color: 'red',
  linewidth: 10, // px
  resolution: new THREE.Vector2(640, 480) // resolution of the viewport
  // dashed, dashScale, dashSize, gapSize
});

const myLine = new Line2(geometry, matLine);

myLine.computeLineDistances();
```

[npm-img]: https://img.shields.io/npm/v/three-fatline
[npm-url]: https://npmjs.org/package/three-fatline
[build-size-img]: https://img.shields.io/bundlephobia/minzip/three-fatline
[build-size-url]: https://bundlephobia.com/result?p=three-fatline
[npm-downloads-img]: https://img.shields.io/npm/dt/three-fatline
[npm-downloads-url]: https://www.npmtrends.com/three-fatline
