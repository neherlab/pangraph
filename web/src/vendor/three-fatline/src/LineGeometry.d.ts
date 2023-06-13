import {
	Line
} from 'three';

import { LineSegmentsGeometry } from './LineSegmentsGeometry';

export declare class LineGeometry extends LineSegmentsGeometry {
	constructor();
	readonly isLineGeometry: true;

	fromLine(line: Line): this;
}