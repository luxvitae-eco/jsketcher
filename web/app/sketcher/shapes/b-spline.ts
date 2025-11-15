import { EndPoint } from "./point";
import { Segment } from "./segment";
import Vector from "math/vector";
import { SketchObject } from "./sketch-object";
import { Layer, Viewer } from "../viewer2d";
import { TOLERANCE, areEqual, arePointsEqual } from "math/equality";
import { lu_solve } from "math/optim/dogleg";
import { isPointInsidePolygon, polygonOffset, ConvexHull2D } from "geom/euclidean";

type IPolynomialFunc = (t: number) => number;
type IPoint = { x: number; y: number; z?: number };
export const getDividedValue = (numerator: number, denominator: number) => {
  if (denominator === 0) {
    return 0;
  } else {
    return numerator / denominator;
  }
};

export class BSplinePolynomial {
  /** * B-spline polynomial with variable coefficients */
  readonly order: number;
  kValues: number[];
  maxIndex: number;
  private cache: Map<number, IPolynomialFunc> = new Map();
  private polynomialArray: BSplinePolynomial[] = [];
  /**
   * @param kValues Polynomial value boundary points Node vector
   * @param order The degree of the polynomial (for example, 3rd order is 2nd order, 5th order is 4th order) degree = order - 1 degree is the degree of the B-spline curve
   */
  constructor(kValues: number[], order: number) {
    this.order = order;
    this.kValues = kValues;
    this.maxIndex = kValues.length - this.order;
    this.polynomialArray = [];
    for (let i = 0; i < this.order; i += 1) {
      this.polynomialArray.push(new BSplinePolynomial(this.kValues, i));
    }
  }

  updateKValues(kValues: number[]) {
    this.kValues = kValues;
    this.maxIndex = kValues.length - this.order;

    this.cache.clear();

    if (this.polynomialArray.length === 0) {
      for (let i = 0; i < this.order; i += 1) {
        this.polynomialArray.push(new BSplinePolynomial(this.kValues, i));
      }
    } else {
      for (let i = 0; i < this.order; i += 1) {
        this.polynomialArray[i].updateKValues(this.kValues);
      }
    }
  }

  get(index: number): IPolynomialFunc {
    if (index > this.maxIndex) {
      return () => 0;
    }
    const cacheFunc = this.cache.get(index);
    if (cacheFunc) {
      return cacheFunc;
    }
    return this.getPolynomialFunc(index);
  }

  /**
   * Get the polynomial evaluation function, that is, the vector polynomial variable coefficient,
   * and return a function that accepts the t parameter
   */
  private getPolynomialFunc(index: number): IPolynomialFunc {
    const tList = this.kValues;
    const { order } = this;
    const polynomialIndexSubtractOne = this.polynomialArray[order - 1];
    let func: IPolynomialFunc;
    if (order === 1) {
      func = (t: number) => {
        if (t >= this.kValues[index] && t < this.kValues[index + 1]) {
          return 1;
        } else {
          return 0;
        }
      };
    } else {
      const k1 = tList[index + order - 1] - tList[index];
      const k2 = tList[index + order] - tList[index + 1];
      func = (t: number) =>
        getDividedValue(t - tList[index], k1) * polynomialIndexSubtractOne.get(index)(t) +
        getDividedValue(tList[index + order] - t, k2) * polynomialIndexSubtractOne.get(index + 1)(t);
    }
    this.cache.set(index, func);
    return func;
  }
}

// B-spline curve interpolation drawing method and control point drawing method
class BSplineInterpolation {
  // Interpolation only supports cubic B-spline
  interpolation: boolean;
  degree: number;
  kSolver: CentripetalParameterMethod;
  cSolver: CPointsCalculator;
  cPoints: EndPoint[];
  kValues: number[];
  fPoints: EndPoint[];
  constructor(degree: number, interpolation: boolean) {
    this.degree = degree;
    this.interpolation = interpolation;
    this.init();
  }

  init() {
    this.fPoints = [];
    this.cPoints = [];
    this.kValues = [];
    if (this.interpolation) {
      this.kSolver = new CentripetalParameterMethod();
      this.cSolver = new CPointsCalculator();
    } else {
      this.kSolver = null;
      this.cSolver = null;
    }
  }

  update(fPoints: EndPoint[]) {
    if (!this.interpolation) {
      return;
    }
    this.fPoints = fPoints;
    if (this.fPoints.length < 3) {
      this.interpolation = false;
    }
  }

  solve() {
    // Solve the nodes and control points according to the interpolation points fPoints
    if (!this.interpolation) {
      return;
    }
    this.kValues = this.kSolver.calculate(this.fPoints, this.degree);
    this.cSolver.setup(this.kValues, this.fPoints, this.degree);
    const cPointsCoordinates = this.cSolver.calculate();
    this.cPoints.length = cPointsCoordinates.length;
    cPointsCoordinates.forEach((item, index) => {
      this.cPoints[index] = new EndPoint(item.x, item.y);
    });
  }
}

class BSplineControlVertices {
  // B-spline curve control point drawing method, support drawing spline curves of different degrees
  CVModel: boolean;
  maxDegree: number;
  degree: number;
  cPoints: EndPoint[];
  kValues: number[];
  fPoints: EndPoint[];
  constructor(degree: number, CVModel: boolean) {
    this.maxDegree = degree;
    this.CVModel = CVModel;
    this.fPoints = [];
    this.cPoints = [];
    this.kValues = [];
  }

  update(cPoints: EndPoint[], degree: number) {
    if (!this.CVModel) {
      return;
    }
    this.cPoints = cPoints;
    if (this.cPoints.length < 3) {
      this.CVModel = false;
    }
    this.maxDegree = degree;
    if (this.cPoints.length < this.maxDegree + 1) {
      this.degree = this.cPoints.length - 1;
    } else {
      this.degree = this.maxDegree;
    }
  }

  solve() {
    // Solve the nodes, the number of which meets the control point and order requirements,
    // and fill 0 and 1 at both ends to make the spline curve clamped and evenly segmented in the middle.
    this.fPoints = [this.cPoints[0], this.cPoints[this.cPoints.length - 1]];
    this.kValues = [...new Array(this.degree + 1).fill(0.0)];
    for (let i = 0; i < this.cPoints.length - this.degree - 1; ++i) {
      this.kValues.push((i + 1) / (this.cPoints.length - this.degree));
    }
    this.kValues.push(...new Array(this.degree + 1).fill(1.0));
  }
}

export interface IBSplineOpts {
  degree: number;
  cPoints: IPoint[];
  fPoints: IPoint[];
  kValues: number[];
}

export class BSpline extends SketchObject {
  ctx: CanvasRenderingContext2D | undefined;

  scale: number;

  degree: number;

  closed: boolean;

  order: number;

  cPoints: EndPoint[]; // Spline control points

  kValues: number[]; // Spline Nodes

  knots: number[]; // Curve nodes after deduplication

  fPoints: EndPoint[]; // Fitting points for easy curve adjustment

  a: EndPoint; // Start point of the curve

  b: EndPoint; // End point of the curve

  numberOfKnots: number;

  numberOfControlPoints: number;

  numberOfFitPoints: number;

  bSplinePolynomial: BSplinePolynomial;

  derivativePolynomial: BSplinePolynomial;

  bSplineInterpolation: BSplineInterpolation;

  bSplineControlVertices: BSplineControlVertices;

  hull: Vector[]; // Curved polygonal bounding box

  // discretePoints: EndPoint[]; // Fixed curve discrete points for distance calculation
  discretePointsWithScale: { [key: number]: EndPoint[] }; // Record discrete points at different ratios to save computing resources

  step: number; // 曲线离散步长

  dragging: boolean = false;

  constructor(
    opts: IBSplineOpts,
    interpolation: boolean = false, // If true, the interpolation method is manually drawn
    CVModel: boolean = false, // If true, the CV method is used for manual drawing. If both are false, the data is read and drawn.
    id?: string,
    ctx?: CanvasRenderingContext2D,
    scale?: number,
  ) {
    super(id);
    this.ctx = ctx;
    this.scale = scale || 1;
    this.degree = opts.degree;
    this.order = this.degree + 1;
    this.closed = false;
    this.numberOfControlPoints = opts.cPoints.length;
    this.numberOfKnots = opts.kValues.length;
    this.numberOfFitPoints = opts.fPoints.length;
    if (arePointsEqual(opts.cPoints[0], opts.cPoints[this.numberOfControlPoints - 1], TOLERANCE)) {
      this.closed = true;
    }
    this.cPoints = [];
    for (const [i, point] of opts.cPoints.entries()) {
      const cPointId = "spline${this.id}_cPoint${i}";
      const cPoint = new EndPoint(point.x, point.y, cPointId);
      this.addChild(cPoint);
      this.cPoints.push(cPoint);
      cPoint.visible = false;
    }
    this.kValues = opts.kValues;
    this.updateKnots();
    this.fPoints = [];
    if (opts.fPoints.length) {
      for (const [i, point] of opts.fPoints.entries()) {
        const fPointId = "spline${this.id}_fPoint${i}";
        const fPoint = new EndPoint(point.x, point.y, fPointId);
        this.addChild(fPoint);
        this.fPoints.push(fPoint);
      }
      this.a = this.fPoints[0];
      this.b = this.fPoints[this.numberOfFitPoints - 1];
    } else {
      this.a = this.cPoints[0];
      this.b = this.cPoints[this.numberOfControlPoints - 1];
      this.a.visible = true;
      this.b.visible = true;
    }
    if (this.degree >= this.numberOfControlPoints) {
      throw new Error(
        "the degree(${this.degree}) should be smaller than the length of control point(${this.numberOfControlPoints}).",
      );
    }
    if (this.degree < 1) {
      throw new Error("degree cannot be less than 1.");
    }
    if (this.numberOfKnots !== this.numberOfControlPoints + this.order) {
      throw new Error(
        "the array length of parameter t (${this.numberOfKnots}) must be equal to the sum of the length of cPoints (${this.numberOfControlPoints}) and the degree (${this.degree}). and 1",
      );
    }
    this.bSplinePolynomial = new BSplinePolynomial(this.kValues, this.order);
    const newKValues = this.kValues.slice(1, this.kValues.length - 1);
    this.derivativePolynomial = new BSplinePolynomial(newKValues, this.order - 1);
    this.bSplineInterpolation = new BSplineInterpolation(this.degree, interpolation);
    this.bSplineControlVertices = new BSplineControlVertices(this.degree, CVModel);
    if (interpolation) {
      this.bSplineInterpolation.update(this.fPoints);
    }
    if (CVModel) {
      this.bSplineControlVertices.update(this.cPoints, this.degree);
    }
    this.step = 0.1;
    this.discretePointsWithScale = {
      1: this.transToEndPoints(this.getDiscretePoints(1)),
    };
  }

  updateKnots() {
    this.knots = Array.from(new Set(this.kValues)).sort((a, b) => a - b);
  }

  getPoint(t: number) {
    let x = 0;
    let y = 0;
    for (let index = 0; index < this.numberOfControlPoints; ++index) {
      const ratio = this.bSplinePolynomial.get(index)(t);
      x += ratio * this.cPoints[index].x;
      y += ratio * this.cPoints[index].y;
    }
    return { x, y };
  }

  basisFunction(i, p, u, knots) {
    if (p === 0) {
      return knots[i] <= u && u < knots[i + 1] ? 1.0 : 0.0;
    }
    const left = (u - knots[i]) / (knots[i + p] - knots[i]) || 0;
    const right = (knots[i + p + 1] - u) / (knots[i + p + 1] - knots[i + 1]) || 0;
    return left * this.basisFunction(i, p - 1, u, knots) + right * this.basisFunction(i + 1, p - 1, u, knots);
  }

  derivativeBSpline(t: number) {
    const n = this.cPoints.length - 1;
    let dx = 0,
      dy = 0;
    for (let i = 0; i < n; i++) {
      const denom = this.kValues[i + this.degree + 1] - this.kValues[i + 1];
      if (denom === 0) continue;
      const coeff = this.degree / denom;
      const diffX = this.cPoints[i + 1].x - this.cPoints[i].x;
      const diffY = this.cPoints[i + 1].y - this.cPoints[i].y;
      const Ni = this.derivativePolynomial.get(i)(t);
      dx += coeff * diffX * Ni;
      dy += coeff * diffY * Ni;
    }
    return { x: dx, y: dy };
  }

  addChildPoint(point: EndPoint): void {
    point.id = this.id;
    this.addChild(point);
  }

  removeChildPoint(point: EndPoint) {
    this.children.forEach((item, index) => {
      if (item === point) {
        this.children.splice(index, 1);
      }
    });
  }

  setChildPoint(points: EndPoint[]) {
    this.children = points;
    points.forEach((item) => {
      item.parent = this;
    });
  }

  updatePoint(index: number, point: EndPoint) {
    if (index < 0 || index > this.cPoints.length - 1) {
      throw new Error("parameter index error.");
    }
    if (typeof this.cPoints[index] !== "undefined") {
      this.removeChildPoint(this.cPoints[index]);
    }
    this.addChildPoint(point);
    this.cPoints[index] = point;
  }

  addCPoint(point: EndPoint) {
    this.addChildPoint(point);
    this.cPoints.push(point);
    point.visible = false;
    this.numberOfControlPoints += 1;
  }

  setCPoints(points: EndPoint[]) {
    this.cPoints.length = points.length;
    this.numberOfControlPoints = this.cPoints.length;
    points.forEach((item, index) => {
      this.updatePoint(index, item);
      item.visible = false;
    });
  }

  resetCPoints(points: EndPoint[], visible: boolean) {
    this.cPoints = points;
    this.numberOfControlPoints = this.cPoints.length;
    this.cPoints.forEach((item) => {
      item.visible = visible;
    });
  }

  addFPoint(point: EndPoint) {
    if (point.id !== this.fPoints[this.fPoints.length - 1].id) {
      this.addChild(point);
      this.fPoints.push(point);
      this.numberOfFitPoints = this.fPoints.length;
    }
  }

  removeFPoint() {
    const point = this.fPoints.pop();
    this.removeChildPoint(point);
  }

  setFPoint(points: EndPoint[]) {
    this.fPoints.length = points.length;
    this.numberOfFitPoints = this.fPoints.length;
    for (let i = 0; i < points.length; ++i) {
      if (this.fPoints[i] !== points[i]) {
        this.removeChildPoint(this.fPoints[i]);
        this.fPoints[i] = points[i];
        this.addChild(this.fPoints[i]);
        this.fPoints[i].visible = true;
      }
    }
  }

  resetFPoints(points: EndPoint[], visible: boolean) {
    this.fPoints = points;
    this.numberOfFitPoints = this.fPoints.length;
    this.fPoints.forEach((item) => {
      item.visible = visible;
    });
  }

  setPointA(point: EndPoint) {
    this.a = point;
  }

  setPointB(point: EndPoint) {
    this.b = point;
  }

  setKValues(kValues: number[]) {
    this.kValues = kValues;
    this.updateKnots();
    this.bSplinePolynomial.updateKValues(this.kValues);
    const newKValues = this.kValues.slice(1, this.kValues.length - 1);
    this.derivativePolynomial.updateKValues(newKValues);
    this.numberOfKnots = this.kValues.length;
  }

  interpolate(fPoints: EndPoint[]) {
    if (fPoints.length > 2) {
      this.bSplineInterpolation.update(fPoints);
    }
  }

  update() {
    this.bSplineInterpolation.solve();
    this.resetFPoints(this.bSplineInterpolation.fPoints, true);
    this.setKValues(this.bSplineInterpolation.kValues);
    this.resetCPoints(this.bSplineInterpolation.cPoints, false);
    this.setPointA(this.fPoints[0]);
    this.setPointB(this.fPoints[this.fPoints.length - 1]);
    this.setChildPoint([...this.fPoints, ...this.cPoints]);
  }

  cvReset(cPoints: EndPoint[], degree: number) {
    if (cPoints.length > 2) {
      this.bSplineControlVertices.update(cPoints, degree);
      this.degree = degree;
      this.order = this.degree + 1;
    }
  }

  cvUpdate() {
    this.bSplineControlVertices.solve();
    this.resetFPoints(this.bSplineControlVertices.fPoints, false);
    this.setKValues(this.bSplineControlVertices.kValues);
    this.resetCPoints(this.bSplineControlVertices.cPoints, true);
    this.setPointA(this.cPoints[0]);
    this.setPointB(this.cPoints[this.cPoints.length - 1]);
    this.setChildPoint([...this.cPoints]);
  }

  getDiscretePoints(scale: number) {
    const ratio = 1 / scale;
    const discretePoints: IPoint[] = [];
    for (let index = 0; index < this.knots.length - 1; ++index) {
      if (index === 0) {
        discretePoints.push(this.a);
      } else {
        const fPoint = this.getPoint(this.knots[index]);
        if (this.bSplineInterpolation.interpolation && this.knots.length === this.fPoints.length) {
          discretePoints.push(...this.transToIPoints([this.fPoints[index]]));
        } else {
          discretePoints.push(fPoint);
        }
      }

      if (this.knots[index + 1] - this.knots[index] < this.step * ratio) {
        continue;
      }
      for (let k = this.knots[index] + this.step * ratio; k <= this.knots[index + 1]; k += this.step * ratio) {
        const p = this.getPoint(k);
        discretePoints.push(p);
      }
    }
    discretePoints.push(this.b);
    return discretePoints;
  }

  visitParams(callback) {
    for (const point of this.cPoints) {
      point.visitParams(callback);
    }
  }

  normalDistance(aim: Vector, scale: number) {
    // Get the vertices of the convex polygon surrounded by control points in sequence
    const boundaryPoints = [...this.cPoints]; // Deep copy avoids ConvexHull2D function sorting affecting this.cPoints
    const hullPoints = ConvexHull2D(boundaryPoints);

    // Get the point vector after the convex polygon is expanded
    // (the center point position of the convex polygon quadrilateral bounding box remains unchanged)
    this.hull = polygonOffset(hullPoints, 1 + 0.3 / scale);
    if (isPointInsidePolygon(aim, this.hull)) {
      const discreteScale = this.getDiscreteScale(scale);
      return this.closestNormalDistance(aim, this.discretePointsWithScale[discreteScale]);
    }
    return -1;
  }

  closestNormalDistance(aim: Vector, segments: EndPoint[]) {
    let hero = -1;
    for (let p = segments.length - 1, q = 0; q < segments.length; p = q++) {
      const dist = Math.min(Segment.calcNormalDistance(aim, segments[p], segments[q]));
      if (dist !== -1) {
        hero = hero === -1 ? dist : Math.min(dist, hero);
      }
    }
    return hero;
  }

  transToEndPoints(points: IPoint[]) {
    const endPoints = [];
    for (const point of points) {
      endPoints.push(new EndPoint(point.x, point.y));
    }
    return endPoints;
  }

  transToIPoints(points: EndPoint[]) {
    const IPoints = [];
    for (const point of points) {
      IPoints.push({ x: point.x, y: point.y, z: 0.0 });
    }
    return IPoints;
  }

  bsplineToBezierSegments() {
    // Each interval is a Bézier
    const bezierSegments = [];
    if (!this.closed && this.bSplineInterpolation) {
      for (let i = 0; i < this.fPoints.length - 1; i++) {
        const a = { x: this.fPoints[i].x, y: this.fPoints[i].y };
        const b = { x: this.fPoints[i + 1].x, y: this.fPoints[i + 1].y };
        const derivative1 = this.derivativeBSpline(this.knots[i]);
        const derivative2 = this.derivativeBSpline(this.knots[i + 1]);
        let cp1X;
        let cp1Y;
        let cp2X;
        let cp2Y;
        if (areEqual(derivative1.x, 0, TOLERANCE)) {
          cp1X = this.fPoints[i].x;
          if (areEqual(this.cPoints[i + 2].x, this.cPoints[i + 1].x, TOLERANCE)) {
            cp1Y = this.fPoints[i].y;
          } else {
            cp1Y =
              ((this.cPoints[i + 2].y - this.cPoints[i + 1].y) / (this.cPoints[i + 2].x - this.cPoints[i + 1].x)) *
                (cp1X - this.cPoints[i + 1].x) +
              this.cPoints[i + 1].y;
          }
        } else {
          const k = derivative1.y / derivative1.x;
          if (areEqual(this.cPoints[i + 2].x, this.cPoints[i + 1].x, TOLERANCE)) {
            cp1X = this.cPoints[i + 1].x;
            cp1Y = this.fPoints[i].y + k * (cp1X - this.fPoints[i].x);
          } else {
            const k1 =
              (this.cPoints[i + 2].y - this.cPoints[i + 1].y) / (this.cPoints[i + 2].x - this.cPoints[i + 1].x);
            cp1X = areEqual(k, k1, TOLERANCE)
              ? this.fPoints[i].x
              : (this.fPoints[i].y - this.cPoints[i + 1].y + k1 * this.cPoints[i + 1].x - k * this.fPoints[i].x) /
                (k1 - k);
            cp1Y =
              ((this.cPoints[i + 2].y - this.cPoints[i + 1].y) / (this.cPoints[i + 2].x - this.cPoints[i + 1].x)) *
                (cp1X - this.cPoints[i + 1].x) +
              this.cPoints[i + 1].y;
          }
        }
        if (areEqual(derivative2.x, 0, TOLERANCE)) {
          cp2X = this.fPoints[i + 1].x;
          if (areEqual(this.cPoints[i + 2].x, this.cPoints[i + 1].x, TOLERANCE)) {
            cp2Y = this.fPoints[i + 1].y;
          } else {
            cp2Y =
              ((this.cPoints[i + 2].y - this.cPoints[i + 1].y) / (this.cPoints[i + 2].x - this.cPoints[i + 1].x)) *
                (cp2X - this.cPoints[i + 1].x) +
              this.cPoints[i + 1].y;
          }
        } else {
          const k = derivative2.y / derivative2.x;
          if (areEqual(this.cPoints[i + 2].x, this.cPoints[i + 1].x, TOLERANCE)) {
            cp2X = this.cPoints[i + 1].x;
            cp2Y = this.fPoints[i + 1].y + k * (cp2X - this.fPoints[i + 1].x);
          } else {
            const k1 =
              (this.cPoints[i + 2].y - this.cPoints[i + 1].y) / (this.cPoints[i + 2].x - this.cPoints[i + 1].x);
            cp2X = areEqual(k, k1, TOLERANCE)
              ? this.fPoints[i + 1].x
              : (this.fPoints[i + 1].y -
                  this.cPoints[i + 1].y +
                  k1 * this.cPoints[i + 1].x -
                  k * this.fPoints[i + 1].x) /
                (k1 - k);
            cp2Y =
              ((this.cPoints[i + 2].y - this.cPoints[i + 1].y) / (this.cPoints[i + 2].x - this.cPoints[i + 1].x)) *
                (cp2X - this.cPoints[i + 1].x) +
              this.cPoints[i + 1].y;
          }
        }
        const cp1 = { x: cp1X, y: cp1Y };
        const cp2 = { x: cp2X, y: cp2Y };
        bezierSegments.push({ a, b, cp1, cp2 });
      }
    }
    return bezierSegments;
  }

  /**
   * Draw B-spline curves (converted to Bézier segments)
   */
  drawBSplineBezier(ctx: CanvasRenderingContext2D) {
    const segments = this.bsplineToBezierSegments();

    ctx.beginPath();

    for (const seg of segments) {
      ctx.moveTo(seg.a.x, seg.a.y);
      ctx.bezierCurveTo(seg.cp1.x, seg.cp1.y, seg.cp2.x, seg.cp2.y, seg.b.x, seg.b.y);
    }
    ctx.stroke();
  }

  drawBSplineLine(ctx: CanvasRenderingContext2D, scale: number) {
    const discretePoints = this.getDiscretePointsWithScale(scale);
    const len = discretePoints.length;
    if (len === 0) {
      return;
    }
    ctx.beginPath();
    ctx.moveTo(discretePoints[0].x, discretePoints[0].y);
    for (const point of discretePoints) {
      ctx.lineTo(point.x, point.y);
    }
    ctx.stroke();
  }

  drawImpl(ctx: CanvasRenderingContext2D, scale: number, viewer: Viewer) {
    // This function will be called multiple times to draw the image.
    const discreteScale = this.getDiscreteScale(scale);
    if (this.bSplineInterpolation.interpolation && this.dragging !== true) {
      this.update();
    } else if (this.bSplineControlVertices.CVModel && this.dragging !== true) {
      this.cvUpdate();
    } else {
      if (!this.discretePointsWithScale[discreteScale]) {
        this.discretePointsWithScale[discreteScale] = this.transToEndPoints(this.getDiscretePoints(discreteScale));
      }
    }
    // this.drawBSplineLine(ctx, discreteScale);
    this.drawBSplineBezier(ctx);
  }

  getDiscreteScale(scale: number) {
    let discreteScale = Math.ceil(Math.log2(scale));
    if (discreteScale < -3) {
      discreteScale = 0.125;
    } else if (discreteScale <= 1) {
      discreteScale = 2 ** discreteScale;
    } else {
      discreteScale = 2;
    }
    return discreteScale;
  }

  getDiscretePointsWithScale(scale: number) {
    const discreteScale = this.getDiscreteScale(scale);
    this.discretePointsWithScale[discreteScale] = this.transToEndPoints(this.getDiscretePoints(discreteScale));
    return this.discretePointsWithScale[discreteScale];
  }

  write() {
    return {
      degree: this.degree,
      cPoints: this.transToIPoints(this.cPoints),
      fPoints: this.transToIPoints(this.fPoints),
      kValues: this.kValues,
    };
  }

  static read(id: string, bSplineData: IBSplineOpts) {
    return new BSpline(bSplineData, false, false, id);
  }

  drag(x, y, dx, dy) {
    this.dragging = true;
    this.translate(dx, dy);
  }

  stabilize(viewer: Viewer) {
    this.children.forEach((c) => c.stabilize(viewer));
  }
}

interface KnotsCalculator {
  calculate(modelPoints: EndPoint[], degree: number): number[];
}

export class CentripetalParameterMethod implements KnotsCalculator {
  calculate(modelPoints: EndPoint[], degree: number) {
    const n = modelPoints.length;
    const accumulatedLengths = [0.0];
    const knotValues = new Array(degree).fill(0.0);
    for (let i = 0; i < n - 1; ++i) {
      const lineLength = Math.sqrt(
        (modelPoints[i + 1].x - modelPoints[i].x) ** 2 + (modelPoints[i + 1].y - modelPoints[i].y) ** 2,
      );
      accumulatedLengths.push(accumulatedLengths[accumulatedLengths.length - 1] + Math.sqrt(lineLength));
    }
    for (let i = 0; i < n; ++i) {
      knotValues.push(accumulatedLengths[i] / accumulatedLengths[n - 1]);
    }
    knotValues.push(...new Array(degree).fill(1.0));
    return knotValues;
  }
}

export class CPointsCalculator {
  cPoints: Array<{ x: number; y: number; z: 0 }>;
  knotValues: number[];
  degree: number;
  modelPoints: EndPoint[];

  constructor() {
    this.cPoints = [];
    this.knotValues = [];
    this.degree = 3;
    this.modelPoints = [];
  }

  setup(knotValues: number[], modelPoints: EndPoint[], degree: number) {
    this.knotValues = knotValues;
    this.degree = degree;
    this.modelPoints = modelPoints;
    const x = this.modelPoints.length;
    if (x < this.degree) {
      throw new Error("too less points !");
    }
  }

  calculate() {
    const n = this.modelPoints.length + this.degree - 1;
    const matrixN = new Array(n);
    const polynomial = new BSplinePolynomial(this.knotValues, this.degree + 1);
    const start = new BesselTangentMethod();
    start.calculate(this.modelPoints[0], this.modelPoints[1], this.modelPoints[2]);
    const end = new BesselTangentMethod();
    end.calculate(
      this.modelPoints[this.modelPoints.length - 3],
      this.modelPoints[this.modelPoints.length - 2],
      this.modelPoints[this.modelPoints.length - 1],
    );
    const { startTangent } = start;
    const { endTangent } = end;
    const matrixP = new Array(n);
    const matrixFX = [(startTangent.x * (this.knotValues[this.degree + 1] - this.knotValues[1])) / this.degree];
    const matrixFY = [(startTangent.y * (this.knotValues[this.degree + 1] - this.knotValues[1])) / this.degree];
    matrixN[0] = [-1, 1, ...new Array(n - 2).fill(0)];
    for (let i = 1; i < this.modelPoints.length; ++i) {
      matrixN[i] = new Array(n);
      matrixFX.push(this.modelPoints[i - 1].x);
      matrixFY.push(this.modelPoints[i - 1].y);
      for (let j = 0; j < n; ++j) {
        matrixN[i][j] = polynomial.get(j)(this.knotValues[i - 1 + this.degree]);
      }
    }
    matrixN[this.modelPoints.length] = [...new Array(n - 1).fill(0), 1];
    matrixFX.push(this.modelPoints[this.modelPoints.length - 1].x);
    matrixFY.push(this.modelPoints[this.modelPoints.length - 1].y);
    matrixFX.push((endTangent.x * (this.knotValues[this.degree + n - 1] - this.knotValues[n - 1])) / this.degree);
    matrixFY.push((endTangent.y * (this.knotValues[this.degree + n - 1] - this.knotValues[n - 1])) / this.degree);
    matrixN[n - 1] = [...new Array(n - 2).fill(0), -1, 1];
    const matrixPX = lu_solve(matrixN, matrixFX, false);
    const matrixPY = lu_solve(matrixN, matrixFY, false);
    this.cPoints = [];
    for (let i = 0; i < n; ++i) {
      this.cPoints[i] = { x: matrixPX[i], y: matrixPY[i], z: 0 };
    }
    return this.cPoints;
  }
}

class BesselTangentMethod {
  startTangent: EndPoint;
  middleTangent: EndPoint;
  endTangent: EndPoint;

  calculate(pointA: EndPoint, pointB: EndPoint, pointC: EndPoint) {
    const distanceAB = Math.sqrt((pointB.x - pointA.x) ** 2 + (pointB.y - pointA.y) ** 2);
    const distanceBC = Math.sqrt((pointC.x - pointB.x) ** 2 + (pointC.y - pointB.y) ** 2);
    const sum = distanceAB + distanceBC;
    const deltaAB = new EndPoint((pointB.x - pointA.x) / distanceAB, (pointB.y - pointA.y) / distanceAB);
    const deltaBC = new EndPoint((pointC.x - pointB.x) / distanceBC, (pointC.y - pointB.y) / distanceBC);
    this.middleTangent = new EndPoint(
      (distanceAB / sum) * deltaAB.x + (distanceBC / sum) * deltaBC.x,
      (distanceAB / sum) * deltaAB.y + (distanceBC / sum) * deltaBC.y,
    );
    this.startTangent = new EndPoint(2 * deltaAB.x - this.middleTangent.x, 2 * deltaAB.y - this.middleTangent.y);
    this.endTangent = new EndPoint(2 * deltaBC.x - this.middleTangent.x, 2 * deltaBC.y - this.middleTangent.y);
  }
}
